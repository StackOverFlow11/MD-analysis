"""Surface charge density and per-atom charge analysis from Bader output."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms

logger = logging.getLogger(__name__)

from ...utils._io_helpers import _cumulative_average, _write_csv
from ...utils.BaderParser import load_bader_atoms
from ...utils.config import (
    AREA_VECTOR_INDICES,
    AXIS_MAP,
    CHARGE_METHOD_COUNTERION,
    CHARGE_METHOD_LAYER,
    DEFAULT_LAYER_TOL_A,
)
from ...utils.StructureParser.LayerParser import (
    circular_mean_fractional,
    detect_interface_layers,
    mic_delta_fractional,
)
from ...utils.StructureParser.WaterParser import detect_water_molecule_indices
from .config import (
    DEFAULT_ACF_FILENAME,
    DEFAULT_DIR_PATTERN,
    DEFAULT_N_SURFACE_LAYERS,
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
    DEFAULT_SURFACE_CHARGE_CSV_NAME,
    DEFAULT_SURFACE_CHARGE_PNG_NAME,
    E_PER_A2_TO_UC_PER_CM2,
)

_T_VALUE_RE = re.compile(r"_t(\d+)")


def _extract_t_value(dirname: str) -> int:
    """Extract numeric t value from a directory name like ``bader_t50_i0``."""
    m = _T_VALUE_RE.search(dirname)
    return int(m.group(1)) if m else 0


def _sorted_frame_dirs(root: Path, dir_pattern: str) -> list[Path]:
    """Discover and numerically sort frame subdirectories by t value."""
    frame_dirs = sorted(root.glob(dir_pattern), key=lambda p: _extract_t_value(p.name))
    if not frame_dirs:
        raise FileNotFoundError(
            f"No subdirectories matching '{dir_pattern}' in {root}"
        )
    return frame_dirs


# ---------------------------------------------------------------------------
# Single-frame analysis
# ---------------------------------------------------------------------------

_VALID_METHODS = (CHARGE_METHOD_COUNTERION, CHARGE_METHOD_LAYER)


def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = "counterion",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
) -> Atoms:
    """Compute surface charge density for a single frame.

    Parameters
    ----------
    atoms
        ASE Atoms with ``"bader_net_charge"`` per-atom array.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).
    method
        ``"counterion"`` — only non-water, non-metal species contribute to σ.
        ``"layer"`` — sum net charges of the interface-layer metal atoms.
    layer_tol_A
        Layer clustering tolerance in Ångströms for ``detect_interface_layers``.
    n_surface_layers
        Number of metal layers (counted inward from each interface) whose net
        charge is summed to compute σ. Only used when ``method="layer"``.

    Returns
    -------
    Atoms
        The same object, mutated with results in ``atoms.info``.
    """
    if method not in _VALID_METHODS:
        raise ValueError(
            f"method must be one of {_VALID_METHODS}, got {method!r}"
        )
    if method == CHARGE_METHOD_LAYER:
        return _compute_surface_charge_layer(
            atoms, metal_symbols=metal_symbols, normal=normal,
            layer_tol_A=layer_tol_A, n_surface_layers=n_surface_layers,
        )
    return _compute_surface_charge_counterion(
        atoms, metal_symbols=metal_symbols, normal=normal,
        layer_tol_A=layer_tol_A,
    )


def _compute_surface_charge_layer(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
) -> Atoms:
    """Surface charge = Σ(net charge of surface-layer atoms) / area."""
    if normal not in AREA_VECTOR_INDICES:
        raise ValueError(
            f"normal must be one of {set(AREA_VECTOR_INDICES)}, got {normal!r}"
        )

    if "bader_net_charge" not in atoms.arrays:
        raise ValueError(
            "atoms.arrays must contain 'bader_net_charge'. "
            "Use load_bader_atoms() first."
        )

    if n_surface_layers < 1:
        raise ValueError(
            f"n_surface_layers must be >= 1, got {n_surface_layers}"
        )

    net_charge = atoms.arrays["bader_net_charge"]

    det = detect_interface_layers(
        atoms, metal_symbols=metal_symbols, normal=normal,
        layer_tol_A=layer_tol_A,
    )

    layers = det.metal_layers_sorted
    n_total = len(layers)
    if n_surface_layers > n_total:
        raise ValueError(
            f"n_surface_layers={n_surface_layers} exceeds total "
            f"metal layers={n_total}"
        )

    aligned_layers = layers[:n_surface_layers]
    opposed_layers = layers[-n_surface_layers:]

    i0, i1 = AREA_VECTOR_INDICES[normal]
    cell = np.asarray(atoms.cell.array, dtype=float)
    area_A2 = float(np.linalg.norm(np.cross(cell[i0], cell[i1])))

    sigma_e_A2 = np.empty(2)
    n_atoms_per_surface = []
    q_per_surface = []
    for i, side_layers in enumerate([aligned_layers, opposed_layers]):
        idx = np.concatenate([np.array(L.atom_indices) for L in side_layers])
        q = float(net_charge[idx].sum())
        sigma_e_A2[i] = q / area_A2
        n_atoms_per_surface.append(len(idx))
        q_per_surface.append(q)

    sigma_uC_cm2 = sigma_e_A2 * E_PER_A2_TO_UC_PER_CM2

    atoms.info["surface_charge_density_e_A2"] = sigma_e_A2.tolist()
    atoms.info["surface_charge_density_uC_cm2"] = sigma_uC_cm2.tolist()
    atoms.info["n_charged_atoms_per_surface"] = n_atoms_per_surface
    atoms.info["charge_per_surface_e"] = q_per_surface

    return atoms


def _compute_surface_charge_counterion(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
) -> Atoms:
    """Surface charge from non-water, non-metal species near each surface."""
    if normal not in AREA_VECTOR_INDICES:
        raise ValueError(
            f"normal must be one of {set(AREA_VECTOR_INDICES)}, got {normal!r}"
        )

    if "bader_net_charge" not in atoms.arrays:
        raise ValueError(
            "atoms.arrays must contain 'bader_net_charge'. "
            "Use load_bader_atoms() first."
        )

    net_charge = atoms.arrays["bader_net_charge"]

    det = detect_interface_layers(
        atoms, metal_symbols=metal_symbols, normal=normal,
        layer_tol_A=layer_tol_A,
    )
    iface_aligned = det.interface_normal_aligned()
    iface_opposed = det.interface_normal_opposed()

    i0, i1 = AREA_VECTOR_INDICES[normal]
    cell = np.asarray(atoms.cell.array, dtype=float)
    area_A2 = float(np.linalg.norm(np.cross(cell[i0], cell[i1])))

    # Exclude water and metal atoms — only solute/counterion species contribute
    water_mol = detect_water_molecule_indices(atoms)
    exclude = set(water_mol.ravel().tolist()) | set(det.metal_indices)

    charged_idx = np.array(
        [i for i in np.where(net_charge != 0.0)[0] if i not in exclude],
        dtype=int,
    )

    # No charged atoms → σ = 0
    if charged_idx.size == 0:
        atoms.info["surface_charge_density_e_A2"] = [0.0, 0.0]
        atoms.info["surface_charge_density_uC_cm2"] = [0.0, 0.0]
        atoms.info["n_charged_atoms_per_surface"] = [0, 0]
        atoms.info["charge_per_surface_e"] = [0.0, 0.0]
        return atoms

    axis_idx = AXIS_MAP[normal]
    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)

    frac_aligned = iface_aligned.center_frac
    frac_opposed = iface_opposed.center_frac

    metal_frac = scaled[list(det.metal_indices), axis_idx]
    midplane = circular_mean_fractional(metal_frac)
    gap_frac = (midplane + 0.5) % 1.0

    frac_charged = scaled[charged_idx, axis_idx]

    gap_dir_aligned = mic_delta_fractional(gap_frac - frac_aligned)
    delta_aligned = mic_delta_fractional(frac_charged - frac_aligned)
    assigned_aligned = (delta_aligned * gap_dir_aligned > 0) & (
        np.abs(delta_aligned) < np.abs(gap_dir_aligned)
    )

    gap_dir_opposed = mic_delta_fractional(gap_frac - frac_opposed)
    delta_opposed = mic_delta_fractional(frac_charged - frac_opposed)
    assigned_opposed = (delta_opposed * gap_dir_opposed > 0) & (
        np.abs(delta_opposed) < np.abs(gap_dir_opposed)
    )

    q_aligned = float(net_charge[charged_idx[assigned_aligned]].sum()) if assigned_aligned.any() else 0.0
    q_opposed = float(net_charge[charged_idx[assigned_opposed]].sum()) if assigned_opposed.any() else 0.0

    sigma_e_A2 = np.array([q_aligned / area_A2, q_opposed / area_A2])
    sigma_uC_cm2 = sigma_e_A2 * E_PER_A2_TO_UC_PER_CM2

    atoms.info["surface_charge_density_e_A2"] = sigma_e_A2.tolist()
    atoms.info["surface_charge_density_uC_cm2"] = sigma_uC_cm2.tolist()
    atoms.info["n_charged_atoms_per_surface"] = [
        int(assigned_aligned.sum()),
        int(assigned_opposed.sum()),
    ]
    atoms.info["charge_per_surface_e"] = [q_aligned, q_opposed]

    return atoms


# ---------------------------------------------------------------------------
# Single-frame indexed atom charges
# ---------------------------------------------------------------------------

def frame_indexed_atom_charges(
    atoms: Atoms,
    atom_indices: np.ndarray,
) -> np.ndarray:
    """Extract net charges for specified atom indices from a single frame.

    Parameters
    ----------
    atoms
        ASE Atoms with ``"bader_net_charge"`` per-atom array.
    atom_indices
        1-D integer array of 0-based atom indices.

    Returns
    -------
    np.ndarray
        Array of shape ``(N, 2)`` where ``[:, 0]`` contains the
        (echoed-back) atom indices and ``[:, 1]`` contains the
        corresponding Bader net charges.
    """
    idx = np.asarray(atom_indices)
    if idx.ndim != 1:
        raise ValueError(
            f"atom_indices must be 1-D, got {idx.ndim}-D"
        )
    if not np.issubdtype(idx.dtype, np.integer):
        raise ValueError(
            f"atom_indices must have integer dtype, got {idx.dtype}"
        )
    if np.any(idx < 0):
        raise ValueError("atom_indices contains negative indices")

    if "bader_net_charge" not in atoms.arrays:
        raise ValueError(
            "atoms.arrays must contain 'bader_net_charge'. "
            "Use load_bader_atoms() first."
        )

    net_charge = atoms.arrays["bader_net_charge"]
    n_atoms = len(atoms)

    oob = idx[idx >= n_atoms]
    if oob.size > 0:
        raise IndexError(
            f"Atom index {int(oob[0])} out of bounds for frame "
            f"with {n_atoms} atoms"
        )

    result = np.empty((len(idx), 2), dtype=float)
    result[:, 0] = idx.astype(float)
    result[:, 1] = net_charge[idx]
    return result


# ---------------------------------------------------------------------------
# Trajectory indexed atom charges
# ---------------------------------------------------------------------------

def trajectory_indexed_atom_charges(
    root_dir: str | Path,
    atom_index_matrix: np.ndarray,
    *,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray:
    """Extract net charges for specified atom indices across trajectory frames.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    atom_index_matrix
        2-D integer array of shape ``(t, N)`` with 0-based atom indices.
        Row *i* lists the atom indices to query in frame *i*.
    dir_pattern
        Glob pattern for discovering frame subdirectories.
    structure_filename
        Name of the structure file (POSCAR) in each subdirectory.
    acf_filename
        Name of the ACF.dat file in each subdirectory.
    potcar_filename
        Name of the POTCAR file in each subdirectory.

    Returns
    -------
    np.ndarray
        Array of shape ``(t, N, 2)`` where ``[:, :, 0]`` contains the
        (echoed-back) atom indices and ``[:, :, 1]`` contains the
        corresponding Bader net charges.
    """
    # --- Validate atom_index_matrix ---
    arr = np.asarray(atom_index_matrix)
    if arr.ndim != 2:
        raise ValueError(
            f"atom_index_matrix must be 2-D, got {arr.ndim}-D"
        )
    if not np.issubdtype(arr.dtype, np.integer):
        raise ValueError(
            f"atom_index_matrix must have integer dtype, got {arr.dtype}"
        )
    if np.any(arr < 0):
        raise ValueError("atom_index_matrix contains negative indices")

    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = _sorted_frame_dirs(root, dir_pattern)

    t_expected = arr.shape[0]
    if t_expected != len(frame_dirs):
        raise ValueError(
            f"atom_index_matrix has {t_expected} rows but found "
            f"{len(frame_dirs)} frame directories"
        )

    t, n = arr.shape
    result = np.empty((t, n, 2), dtype=float)

    for i, frame_dir in enumerate(frame_dirs):
        fname = frame_dir.name
        poscar = frame_dir / structure_filename
        acf = frame_dir / acf_filename
        potcar = frame_dir / potcar_filename

        for path, label in [(poscar, structure_filename),
                            (acf, acf_filename),
                            (potcar, potcar_filename)]:
            if not path.exists():
                raise FileNotFoundError(
                    f"{label} not found in frame {fname}: {path}"
                )

        atoms = load_bader_atoms(poscar, acf, potcar)

        try:
            frame_result = frame_indexed_atom_charges(atoms, arr[i])
        except IndexError as exc:
            # Re-raise with frame name for context
            raise IndexError(
                f"Atom index {int(arr[i][arr[i] >= len(atoms)][0])} "
                f"out of bounds for frame {fname} "
                f"with {len(atoms)} atoms"
            ) from exc

        result[i] = frame_result

    return result


# ---------------------------------------------------------------------------
# Trajectory surface charge density
# ---------------------------------------------------------------------------

def trajectory_surface_charge(
    root_dir: str | Path,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = "counterion",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray:
    """Compute surface charge density across trajectory frames.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).
    dir_pattern
        Glob pattern for discovering frame subdirectories.
    structure_filename
        Name of the structure file (POSCAR) in each subdirectory.
    acf_filename
        Name of the ACF.dat file in each subdirectory.
    potcar_filename
        Name of the POTCAR file in each subdirectory.

    Returns
    -------
    np.ndarray
        Array of shape ``(t, 2)`` where ``[:, 0]`` is σ_aligned
        and ``[:, 1]`` is σ_opposed, in μC/cm².
    """
    if normal not in AREA_VECTOR_INDICES:
        raise ValueError(
            f"normal must be one of {set(AREA_VECTOR_INDICES)}, got {normal!r}"
        )

    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = _sorted_frame_dirs(root, dir_pattern)

    rows = []
    for frame_dir in frame_dirs:
        fname = frame_dir.name
        poscar = frame_dir / structure_filename
        acf = frame_dir / acf_filename
        potcar = frame_dir / potcar_filename

        for path, label in [(poscar, structure_filename),
                            (acf, acf_filename),
                            (potcar, potcar_filename)]:
            if not path.exists():
                raise FileNotFoundError(
                    f"{label} not found in frame {fname}: {path}"
                )

        atoms = load_bader_atoms(poscar, acf, potcar)
        compute_frame_surface_charge(
            atoms, metal_symbols=metal_symbols, normal=normal, method=method,
            layer_tol_A=layer_tol_A, n_surface_layers=n_surface_layers,
        )
        sigma = atoms.info["surface_charge_density_uC_cm2"]
        rows.append(sigma)

    return np.array(rows, dtype=float)


# ---------------------------------------------------------------------------
# Private helpers for analysis output
# ---------------------------------------------------------------------------

def _plot_surface_charge(
    png_path: Path,
    steps: np.ndarray,
    sigma_aligned: np.ndarray,
    sigma_opposed: np.ndarray,
    sigma_aligned_cum: np.ndarray,
    sigma_opposed_cum: np.ndarray,
) -> None:
    """Plot surface charge density with instantaneous and cumulative average."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(steps, sigma_aligned, lw=1.0, alpha=0.65, color="tab:blue", label="aligned (inst.)")
    ax.plot(steps, sigma_aligned_cum, lw=2.0, color="tab:blue", ls="--", label="aligned (cum. avg)")
    ax.plot(steps, sigma_opposed, lw=1.0, alpha=0.65, color="tab:orange", label="opposed (inst.)")
    ax.plot(steps, sigma_opposed_cum, lw=2.0, color="tab:orange", ls="--", label="opposed (cum. avg)")
    ax.set_xlabel("MD step")
    ax.set_ylabel(r"$\sigma$ ($\mu$C/cm$^2$)")
    ax.set_title("Surface charge density")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# End-to-end surface charge analysis
# ---------------------------------------------------------------------------

def surface_charge_analysis(
    root_dir: str | Path = ".",
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = "counterion",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    output_dir: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """End-to-end surface charge density analysis with CSV + PNG output.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).
    dir_pattern
        Glob pattern for discovering frame subdirectories.
    structure_filename, acf_filename, potcar_filename
        Per-frame file names.
    output_dir
        Where to write CSV and PNG. Defaults to *root_dir*.
    frame_start, frame_end, frame_step
        Slice parameters applied to the sorted frame list.
    verbose
        Print progress.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    if normal not in AREA_VECTOR_INDICES:
        raise ValueError(
            f"normal must be one of {set(AREA_VECTOR_INDICES)}, got {normal!r}"
        )

    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = _sorted_frame_dirs(root, dir_pattern)
    frame_dirs = frame_dirs[frame_start:frame_end:frame_step]
    if not frame_dirs:
        raise FileNotFoundError("No frame directories after slicing")
    logger.info("Surface charge analysis: %d frames, method=%s", len(frame_dirs), method)

    if output_dir is None:
        output_dir = root
    output_dir = Path(output_dir)

    steps_list: list[int] = []
    sigma_aligned_list: list[float] = []
    sigma_opposed_list: list[float] = []

    iterator = enumerate(frame_dirs)
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frame_dirs, desc="Surface charge", unit="frame", ascii=" =")
        iterator = enumerate(iterator)

    for idx, frame_dir in iterator:
        fname = frame_dir.name
        poscar = frame_dir / structure_filename
        acf = frame_dir / acf_filename
        potcar = frame_dir / potcar_filename

        for path, label in [(poscar, structure_filename),
                            (acf, acf_filename),
                            (potcar, potcar_filename)]:
            if not path.exists():
                raise FileNotFoundError(
                    f"{label} not found in frame {fname}: {path}"
                )

        atoms = load_bader_atoms(poscar, acf, potcar)
        compute_frame_surface_charge(
            atoms, metal_symbols=metal_symbols, normal=normal, method=method,
            layer_tol_A=layer_tol_A, n_surface_layers=n_surface_layers,
        )
        sigma = atoms.info["surface_charge_density_uC_cm2"]
        step = _extract_t_value(fname)
        steps_list.append(step)
        sigma_aligned_list.append(sigma[0])
        sigma_opposed_list.append(sigma[1])

    steps = np.array(steps_list)
    sigma_aligned = np.array(sigma_aligned_list)
    sigma_opposed = np.array(sigma_opposed_list)
    sigma_aligned_cum = _cumulative_average(sigma_aligned)
    sigma_opposed_cum = _cumulative_average(sigma_opposed)

    # CSV
    fieldnames = [
        "step",
        "sigma_aligned_uC_cm2",
        "sigma_opposed_uC_cm2",
        "sigma_aligned_cumavg_uC_cm2",
        "sigma_opposed_cumavg_uC_cm2",
    ]
    csv_rows: list[dict] = []
    for i in range(len(steps)):
        csv_rows.append({
            "step": int(steps[i]),
            "sigma_aligned_uC_cm2": float(sigma_aligned[i]),
            "sigma_opposed_uC_cm2": float(sigma_opposed[i]),
            "sigma_aligned_cumavg_uC_cm2": float(sigma_aligned_cum[i]),
            "sigma_opposed_cumavg_uC_cm2": float(sigma_opposed_cum[i]),
        })
    csv_path = output_dir / DEFAULT_SURFACE_CHARGE_CSV_NAME
    _write_csv(csv_path, csv_rows, fieldnames)

    # PNG
    png_path = output_dir / DEFAULT_SURFACE_CHARGE_PNG_NAME
    _plot_surface_charge(png_path, steps, sigma_aligned, sigma_opposed,
                         sigma_aligned_cum, sigma_opposed_cum)

    return csv_path
