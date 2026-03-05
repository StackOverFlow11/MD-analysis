"""Surface charge density and per-atom charge analysis from Bader output."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms

from ..utils.BaderParser import load_bader_atoms
from ..utils.LayerParser import detect_interface_layers
from .config import (
    DEFAULT_ACF_FILENAME,
    DEFAULT_DIR_PATTERN,
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
    DEFAULT_SURFACE_CHARGE_CSV_NAME,
    DEFAULT_SURFACE_CHARGE_PNG_NAME,
    E_PER_A2_TO_UC_PER_CM2,
)

# Map normal axis to the two cell-vector indices spanning the surface plane
_AREA_VECTORS = {"a": (1, 2), "b": (0, 2), "c": (0, 1)}

_T_VALUE_RE = re.compile(r"_t(\d+)")


def _extract_t_value(dirname: str) -> int:
    """Extract numeric t value from a directory name like ``calc_t50_i0``."""
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

def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
) -> Atoms:
    """Compute surface charge density for a single frame.

    The input ``atoms`` must already carry ``"bader_net_charge"`` in
    ``atoms.arrays`` (as produced by :func:`load_bader_atoms`).

    Results are stored in-place on the returned Atoms:
    - ``atoms.info["surface_charge_density_e_A2"]`` — [σ_bottom, σ_top] in e/Å²
    - ``atoms.info["surface_charge_density_uC_cm2"]`` — same in μC/cm²

    Parameters
    ----------
    atoms
        ASE Atoms with ``"bader_net_charge"`` per-atom array.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).

    Returns
    -------
    Atoms
        The same object, mutated with results in ``atoms.info``.
    """
    if normal not in _AREA_VECTORS:
        raise ValueError(
            f"normal must be one of {set(_AREA_VECTORS)}, got {normal!r}"
        )

    if "bader_net_charge" not in atoms.arrays:
        raise ValueError(
            "atoms.arrays must contain 'bader_net_charge'. "
            "Use load_bader_atoms() first."
        )

    net_charge = atoms.arrays["bader_net_charge"]

    # Detect interface layers
    det = detect_interface_layers(atoms, metal_symbols=metal_symbols, normal=normal)
    iface_layers = det.interface_layers()
    if len(iface_layers) != 2:
        raise ValueError(
            f"Expected exactly 2 interface layers, got {len(iface_layers)}"
        )

    # Surface area from the two cell vectors spanning the surface plane
    i0, i1 = _AREA_VECTORS[normal]
    cell = np.asarray(atoms.cell.array, dtype=float)
    area_A2 = float(np.linalg.norm(np.cross(cell[i0], cell[i1])))

    # Sort interface layers by center_s (bottom first)
    sorted_iface = sorted(iface_layers, key=lambda L: L.center_s)

    sigma_e_A2 = np.empty(2)
    for i, layer in enumerate(sorted_iface):
        idx = np.array(layer.atom_indices)
        sigma_e_A2[i] = net_charge[idx].sum() / area_A2

    sigma_uC_cm2 = sigma_e_A2 * E_PER_A2_TO_UC_PER_CM2

    atoms.info["surface_charge_density_e_A2"] = sigma_e_A2.tolist()
    atoms.info["surface_charge_density_uC_cm2"] = sigma_uC_cm2.tolist()

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
        Array of shape ``(t, 2)`` where ``[:, 0]`` is σ_bottom and
        ``[:, 1]`` is σ_top, in μC/cm².
    """
    if normal not in _AREA_VECTORS:
        raise ValueError(
            f"normal must be one of {set(_AREA_VECTORS)}, got {normal!r}"
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
            atoms, metal_symbols=metal_symbols, normal=normal
        )
        sigma = atoms.info["surface_charge_density_uC_cm2"]
        rows.append(sigma)

    return np.array(rows, dtype=float)


# ---------------------------------------------------------------------------
# Private helpers for analysis output
# ---------------------------------------------------------------------------

def _cumulative_average(values: np.ndarray) -> np.ndarray:
    csum = np.cumsum(values, dtype=float)
    return csum / np.arange(1, values.size + 1, dtype=float)


def _write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def _plot_surface_charge(
    png_path: Path,
    steps: np.ndarray,
    sigma_bottom: np.ndarray,
    sigma_top: np.ndarray,
    sigma_bottom_cum: np.ndarray,
    sigma_top_cum: np.ndarray,
) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(steps, sigma_bottom, lw=1.0, alpha=0.65, color="tab:blue", label="bottom (inst.)")
    ax.plot(steps, sigma_bottom_cum, lw=2.0, color="tab:blue", ls="--", label="bottom (cum. avg)")
    ax.plot(steps, sigma_top, lw=1.0, alpha=0.65, color="tab:orange", label="top (inst.)")
    ax.plot(steps, sigma_top_cum, lw=2.0, color="tab:orange", ls="--", label="top (cum. avg)")
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
    if normal not in _AREA_VECTORS:
        raise ValueError(
            f"normal must be one of {set(_AREA_VECTORS)}, got {normal!r}"
        )

    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = _sorted_frame_dirs(root, dir_pattern)
    frame_dirs = frame_dirs[frame_start:frame_end:frame_step]
    if not frame_dirs:
        raise FileNotFoundError("No frame directories after slicing")

    if output_dir is None:
        output_dir = root
    output_dir = Path(output_dir)

    steps_list: list[int] = []
    sigma_bottom_list: list[float] = []
    sigma_top_list: list[float] = []

    for idx, frame_dir in enumerate(frame_dirs):
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
            atoms, metal_symbols=metal_symbols, normal=normal,
        )
        sigma = atoms.info["surface_charge_density_uC_cm2"]
        step = _extract_t_value(fname)
        steps_list.append(step)
        sigma_bottom_list.append(sigma[0])
        sigma_top_list.append(sigma[1])

        if verbose:
            print(f"  [{idx + 1}/{len(frame_dirs)}] {fname}: "
                  f"σ_bottom={sigma[0]:.4f}, σ_top={sigma[1]:.4f} μC/cm²")

    steps = np.array(steps_list)
    sigma_bottom = np.array(sigma_bottom_list)
    sigma_top = np.array(sigma_top_list)
    sigma_bottom_cum = _cumulative_average(sigma_bottom)
    sigma_top_cum = _cumulative_average(sigma_top)

    # CSV
    fieldnames = [
        "step",
        "sigma_bottom_uC_cm2",
        "sigma_top_uC_cm2",
        "sigma_bottom_cumavg_uC_cm2",
        "sigma_top_cumavg_uC_cm2",
    ]
    csv_rows: list[dict] = []
    for i in range(len(steps)):
        csv_rows.append({
            "step": int(steps[i]),
            "sigma_bottom_uC_cm2": float(sigma_bottom[i]),
            "sigma_top_uC_cm2": float(sigma_top[i]),
            "sigma_bottom_cumavg_uC_cm2": float(sigma_bottom_cum[i]),
            "sigma_top_cumavg_uC_cm2": float(sigma_top_cum[i]),
        })
    csv_path = output_dir / DEFAULT_SURFACE_CHARGE_CSV_NAME
    _write_csv(csv_path, csv_rows, fieldnames)

    # PNG
    png_path = output_dir / DEFAULT_SURFACE_CHARGE_PNG_NAME
    _plot_surface_charge(png_path, steps, sigma_bottom, sigma_top,
                         sigma_bottom_cum, sigma_top_cum)

    return csv_path
