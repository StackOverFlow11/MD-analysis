"""Surface charge density computation from Bader charges."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms

from ....utils._io_helpers import _cumulative_average, _write_csv
from ....utils.BaderParser import load_bader_atoms
from ....utils.config import (
    AREA_VECTOR_INDICES,
    AXIS_MAP,
    CHARGE_METHOD_COUNTERION,
    CHARGE_METHOD_LAYER,
    DEFAULT_LAYER_TOL_A,
)
from ....utils.StructureParser.LayerParser import (
    circular_mean_fractional,
    detect_interface_layers,
    mic_delta_fractional,
)
from ....utils.StructureParser.WaterParser import detect_water_molecule_indices
from ..config import (
    DEFAULT_ACF_FILENAME,
    DEFAULT_DIR_PATTERN,
    DEFAULT_N_SURFACE_LAYERS,
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
    DEFAULT_SURFACE_CHARGE_CSV_NAME,
    DEFAULT_SURFACE_CHARGE_PNG_NAME,
    E_PER_A2_TO_UC_PER_CM2,
)
from ._frame_utils import _extract_t_value, _sorted_frame_dirs

logger = logging.getLogger(__name__)


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

    q_ci_aligned = float(net_charge[charged_idx[assigned_aligned]].sum()) if assigned_aligned.any() else 0.0
    q_ci_opposed = float(net_charge[charged_idx[assigned_opposed]].sum()) if assigned_opposed.any() else 0.0

    # Surface charge = negative of counterion charge (charge neutrality)
    q_aligned = -q_ci_aligned
    q_opposed = -q_ci_opposed

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
    *,
    phi_aligned: np.ndarray | None = None,
    phi_opposed: np.ndarray | None = None,
    phi_aligned_cum: np.ndarray | None = None,
    phi_opposed_cum: np.ndarray | None = None,
    fit_rmse: float | None = None,
    potential_reference: str = "SHE",
) -> None:
    """Plot surface charge density with instantaneous and cumulative average.

    If *phi_aligned* / *phi_opposed* are provided, a secondary y-axis shows
    the extrapolated electrode potential.
    """
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_phi = phi_aligned is not None

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(steps, sigma_aligned, lw=1.0, alpha=0.65, color="tab:blue", label="aligned σ (inst.)")
    ax.plot(steps, sigma_aligned_cum, lw=2.0, color="tab:blue", ls="--", label="aligned σ (cum. avg)")
    ax.plot(steps, sigma_opposed, lw=1.0, alpha=0.65, color="tab:orange", label="opposed σ (inst.)")
    ax.plot(steps, sigma_opposed_cum, lw=2.0, color="tab:orange", ls="--", label="opposed σ (cum. avg)")
    ax.set_xlabel("MD step")
    ax.set_ylabel(r"$\sigma$ ($\mu$C/cm$^2$)")
    ax.set_title("Surface charge density" + (" + extrapolated potential" if has_phi else ""))
    ax.grid(True, alpha=0.25)

    if has_phi:
        ax2 = ax.twinx()
        ax2.plot(steps, phi_aligned, lw=1.0, alpha=0.45, color="tab:green",
                 label="aligned φ (inst.)")
        ax2.plot(steps, phi_aligned_cum, lw=2.0, color="tab:green", ls="--",
                 label="aligned φ (cum. avg)")
        ax2.plot(steps, phi_opposed, lw=1.0, alpha=0.45, color="tab:red",
                 label="opposed φ (inst.)")
        ax2.plot(steps, phi_opposed_cum, lw=2.0, color="tab:red", ls="--",
                 label="opposed φ (cum. avg)")
        ax2.set_ylabel(f"φ (V vs {potential_reference})")

        # Combine legends from both axes
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)

        # Annotate fit RMSE
        if fit_rmse is not None:
            ax2.annotate(
                f"fit RMSE = {fit_rmse:.4e} V",
                xy=(0.98, 0.02), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=8, fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.2", fc="wheat", alpha=0.5),
            )
    else:
        ax.legend()

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def _plot_single_side_charge(
    png_path: Path,
    steps: np.ndarray,
    sigma: np.ndarray,
    sigma_cum: np.ndarray,
    *,
    side_label: str = "aligned",
    phi: np.ndarray | None = None,
    phi_cum: np.ndarray | None = None,
    fit_rmse: float | None = None,
    potential_reference: str = "SHE",
) -> None:
    """Plot single-side surface charge density (± extrapolated potential)."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_phi = phi is not None

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(steps, sigma, lw=1.0, alpha=0.65, color="tab:blue",
            label=f"{side_label} σ (inst.)")
    ax.plot(steps, sigma_cum, lw=2.0, color="tab:blue", ls="--",
            label=f"{side_label} σ (cum. avg)")
    ax.set_xlabel("MD step")
    ax.set_ylabel(r"$\sigma$ ($\mu$C/cm$^2$)")
    ax.set_title(
        f"Surface charge density ({side_label})"
        + (" + extrapolated potential" if has_phi else "")
    )
    ax.grid(True, alpha=0.25)

    if has_phi:
        ax2 = ax.twinx()
        ax2.plot(steps, phi, lw=1.0, alpha=0.45, color="tab:green",
                 label=f"{side_label} φ (inst.)")
        ax2.plot(steps, phi_cum, lw=2.0, color="tab:green", ls="--",
                 label=f"{side_label} φ (cum. avg)")
        ax2.set_ylabel(f"φ (V vs {potential_reference})")

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)

        if fit_rmse is not None:
            ax2.annotate(
                f"fit RMSE = {fit_rmse:.4e} V",
                xy=(0.98, 0.02), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=8, fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.2", fc="wheat", alpha=0.5),
            )
    else:
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
    potential_reference: str = "SHE",
    potential_pH: float = 0.0,
    potential_temperature_K: float = 298.15,
    potential_phi_pzc: float | None = None,
    target_side: str | None = None,
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
    potential_reference
        Output potential reference scale (``"SHE"``, ``"RHE"``, ``"PZC"``).
    potential_pH
        Solution pH for RHE conversion.
    potential_temperature_K
        Temperature in K for RHE conversion.
    potential_phi_pzc
        Potential of zero charge in V vs SHE for PZC conversion.
    target_side
        ``None`` — output both sides (default).
        ``"aligned"`` or ``"opposed"`` — output only that side.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    if target_side is not None and target_side not in ("aligned", "opposed"):
        raise ValueError(
            f"target_side must be 'aligned', 'opposed', or None, got {target_side!r}"
        )
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

    # --- Optional: extrapolate potential from calibration ---
    phi_aligned: np.ndarray | None = None
    phi_opposed: np.ndarray | None = None
    phi_aligned_cum: np.ndarray | None = None
    phi_opposed_cum: np.ndarray | None = None
    fit_rmse: float | None = None

    ref_tag = potential_reference.upper()

    try:
        from ...calibration._data import load_calibration_json
        from ...calibration._mapper import mapper_from_dict

        cal_data, fit_params = load_calibration_json()
        mapper = mapper_from_dict(fit_params)
        phi_aligned = mapper.predict(sigma_aligned)
        phi_opposed = mapper.predict(sigma_opposed)

        # Convert reference scale if needed
        if ref_tag != "SHE":
            from ...calibration.CalibrationWorkflow import convert_reference

            phi_aligned = convert_reference(
                phi_aligned, from_ref="SHE", to_ref=ref_tag,
                temperature_K=potential_temperature_K, pH=potential_pH,
                phi_pzc=potential_phi_pzc,
            )
            phi_opposed = convert_reference(
                phi_opposed, from_ref="SHE", to_ref=ref_tag,
                temperature_K=potential_temperature_K, pH=potential_pH,
                phi_pzc=potential_phi_pzc,
            )

        phi_aligned_cum = _cumulative_average(phi_aligned)
        phi_opposed_cum = _cumulative_average(phi_opposed)
        fit_rmse = fit_params.get("rmse")
        logger.info(
            "Calibration loaded — extrapolated potential appended "
            "(reference=%s, fit RMSE=%.4e V)", ref_tag, fit_rmse or 0.0,
        )
    except FileNotFoundError:
        logger.info(
            "No calibration file found at default location. "
            "Run 'Charge-Potential Calibration' (menu 23) to enable "
            "automatic potential extrapolation from surface charge density."
        )
    except Exception as exc:
        logger.warning("Failed to load calibration: %s", exc)

    # --- Select side(s) for output ---
    if target_side is not None:
        side_idx = 0 if target_side == "aligned" else 1
        sigma_side = np.array([sigma_aligned, sigma_opposed][side_idx])
        sigma_side_cum = _cumulative_average(sigma_side)
        phi_side = [phi_aligned, phi_opposed][side_idx] if phi_aligned is not None else None
        phi_side_cum = _cumulative_average(phi_side) if phi_side is not None else None

        # CSV (single-side)
        fieldnames = ["step", "sigma_uC_cm2", "sigma_cumavg_uC_cm2"]
        if phi_side is not None:
            fieldnames.extend([f"phi_V_vs_{ref_tag}", f"phi_cumavg_V_vs_{ref_tag}"])

        csv_rows: list[dict] = []
        for i in range(len(steps)):
            row: dict = {
                "step": int(steps[i]),
                "sigma_uC_cm2": float(sigma_side[i]),
                "sigma_cumavg_uC_cm2": float(sigma_side_cum[i]),
            }
            if phi_side is not None:
                row[f"phi_V_vs_{ref_tag}"] = float(phi_side[i])
                row[f"phi_cumavg_V_vs_{ref_tag}"] = float(phi_side_cum[i])
            csv_rows.append(row)
        csv_path = output_dir / DEFAULT_SURFACE_CHARGE_CSV_NAME
        _write_csv(csv_path, csv_rows, fieldnames)

        # PNG (single-side)
        png_path = output_dir / DEFAULT_SURFACE_CHARGE_PNG_NAME
        _plot_single_side_charge(
            png_path, steps, sigma_side, sigma_side_cum,
            side_label=target_side,
            phi=phi_side, phi_cum=phi_side_cum,
            fit_rmse=fit_rmse, potential_reference=ref_tag,
        )

        # Summary
        if phi_side is not None:
            logger.info(
                "Extrapolated potential — %s side (cum. avg over %d frames):\n"
                "  %.6f ± %.6f V vs %s\n"
                "  fit RMSE: %.4e V",
                target_side, len(phi_side),
                phi_side_cum[-1], float(np.std(phi_side)), ref_tag,
                fit_rmse or 0.0,
            )

        return csv_path

    # CSV (both sides)
    fieldnames = [
        "step",
        "sigma_aligned_uC_cm2",
        "sigma_opposed_uC_cm2",
        "sigma_aligned_cumavg_uC_cm2",
        "sigma_opposed_cumavg_uC_cm2",
    ]
    if phi_aligned is not None:
        fieldnames.extend([
            f"phi_aligned_V_vs_{ref_tag}",
            f"phi_opposed_V_vs_{ref_tag}",
            f"phi_aligned_cumavg_V_vs_{ref_tag}",
            f"phi_opposed_cumavg_V_vs_{ref_tag}",
        ])

    csv_rows: list[dict] = []
    for i in range(len(steps)):
        row: dict = {
            "step": int(steps[i]),
            "sigma_aligned_uC_cm2": float(sigma_aligned[i]),
            "sigma_opposed_uC_cm2": float(sigma_opposed[i]),
            "sigma_aligned_cumavg_uC_cm2": float(sigma_aligned_cum[i]),
            "sigma_opposed_cumavg_uC_cm2": float(sigma_opposed_cum[i]),
        }
        if phi_aligned is not None:
            row[f"phi_aligned_V_vs_{ref_tag}"] = float(phi_aligned[i])
            row[f"phi_opposed_V_vs_{ref_tag}"] = float(phi_opposed[i])
            row[f"phi_aligned_cumavg_V_vs_{ref_tag}"] = float(phi_aligned_cum[i])
            row[f"phi_opposed_cumavg_V_vs_{ref_tag}"] = float(phi_opposed_cum[i])
        csv_rows.append(row)
    csv_path = output_dir / DEFAULT_SURFACE_CHARGE_CSV_NAME
    _write_csv(csv_path, csv_rows, fieldnames)

    # PNG (both sides)
    png_path = output_dir / DEFAULT_SURFACE_CHARGE_PNG_NAME
    _plot_surface_charge(
        png_path, steps, sigma_aligned, sigma_opposed,
        sigma_aligned_cum, sigma_opposed_cum,
        phi_aligned=phi_aligned,
        phi_opposed=phi_opposed,
        phi_aligned_cum=phi_aligned_cum,
        phi_opposed_cum=phi_opposed_cum,
        fit_rmse=fit_rmse,
        potential_reference=ref_tag,
    )

    # Print potential summary if available
    if phi_aligned is not None:
        n = len(phi_aligned)
        logger.info(
            "Extrapolated potential (cum. avg over %d frames):\n"
            "  aligned: %.6f ± %.6f V vs %s (σ propagation)\n"
            "  opposed: %.6f ± %.6f V vs %s (σ propagation)\n"
            "  fit RMSE: %.4e V",
            n,
            phi_aligned_cum[-1], float(np.std(phi_aligned)), ref_tag,
            phi_opposed_cum[-1], float(np.std(phi_opposed)), ref_tag,
            fit_rmse or 0.0,
        )

    return csv_path
