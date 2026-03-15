"""Per-atom Bader charge extraction and tracking across trajectory frames."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms

from ....scripts.utils.IndexMapper import read_index_map_from_poscar, remap_array
from ....utils._io_helpers import _cumulative_average, _write_csv
from ....utils.BaderParser import load_bader_atoms
from ....utils.config import (
    AREA_VECTOR_INDICES,
    AXIS_MAP,
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
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
)
from ._frame_utils import _extract_step_and_time, _sorted_frame_dirs
from .BaderData import BaderTrajectoryData, load_bader_trajectory

logger = logging.getLogger(__name__)

# Default output filenames
DEFAULT_TRACKED_CHARGE_CSV = "tracked_atom_charges.csv"
DEFAULT_TRACKED_CHARGE_PNG = "tracked_atom_charges.png"
DEFAULT_COUNTERION_CHARGE_CSV = "counterion_charges.csv"
DEFAULT_COUNTERION_SUMMARY_CSV = "counterion_summary.csv"
DEFAULT_COUNTERION_CHARGE_PNG = "counterion_charges.png"


# ---------------------------------------------------------------------------
# Single-frame indexed atom charges (existing, POSCAR order)
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
# Trajectory indexed atom charges (existing, POSCAR order)
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
# Tracked atom charge analysis (NEW — XYZ order)
# ---------------------------------------------------------------------------

def tracked_atom_charge_analysis(
    root_dir: str | Path = ".",
    *,
    atom_indices_xyz: Iterable[int],
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
    """Track Bader net charges for specified XYZ atoms: time evolution + ensemble average.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    atom_indices_xyz
        Atom indices in original XYZ ordering (0-based).
    dir_pattern
        Glob pattern for frame subdirectories.
    structure_filename, acf_filename, potcar_filename
        Per-frame file names.
    output_dir
        Where to write CSV and PNG.  Defaults to *root_dir*.
    frame_start, frame_end, frame_step
        Slice parameters applied to the sorted frame list.
    verbose
        Show tqdm progress bar.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    indices = np.asarray(list(atom_indices_xyz), dtype=int)
    if indices.ndim != 1 or indices.size == 0:
        raise ValueError("atom_indices_xyz must be a non-empty 1-D sequence of integers")
    if np.any(indices < 0):
        raise ValueError("atom_indices_xyz contains negative indices")

    data = load_bader_trajectory(
        root_dir,
        dir_pattern=dir_pattern,
        structure_filename=structure_filename,
        acf_filename=acf_filename,
        potcar_filename=potcar_filename,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    n_atoms_total = data.net_charges.shape[1]
    oob = indices[indices >= n_atoms_total]
    if oob.size > 0:
        raise IndexError(
            f"Atom index {int(oob[0])} out of bounds for trajectory "
            f"with {n_atoms_total} atoms"
        )

    # Slice charges for requested atoms: (n_frames, n_selected)
    selected_charges = data.net_charges[:, indices]
    n_frames = len(data.steps)

    # Cumulative averages per atom
    cum_avgs = np.empty_like(selected_charges)
    for j in range(len(indices)):
        cum_avgs[:, j] = _cumulative_average(selected_charges[:, j])

    # CSV
    root = Path(root_dir)
    if output_dir is None:
        output_dir = root
    output_dir = Path(output_dir)

    q_fields = [f"atom{int(idx)}_q" for idx in indices]
    cumavg_fields = [f"atom{int(idx)}_cumavg" for idx in indices]
    fieldnames = ["step", "time_fs"] + q_fields + cumavg_fields

    csv_rows: list[dict] = []
    for i in range(n_frames):
        row: dict = {
            "step": int(data.steps[i]),
            "time_fs": int(data.times[i]),
        }
        for j, idx in enumerate(indices):
            row[f"atom{int(idx)}_q"] = float(selected_charges[i, j])
            row[f"atom{int(idx)}_cumavg"] = float(cum_avgs[i, j])
        csv_rows.append(row)

    csv_path = output_dir / DEFAULT_TRACKED_CHARGE_CSV
    _write_csv(csv_path, csv_rows, fieldnames)

    # PNG
    png_path = output_dir / DEFAULT_TRACKED_CHARGE_PNG
    _plot_tracked_charges(png_path, data.steps, indices, selected_charges, cum_avgs)

    logger.info("Tracked atom charge analysis: %d frames, %d atoms → %s",
                n_frames, len(indices), csv_path)
    return csv_path


def _plot_tracked_charges(
    png_path: Path,
    steps: np.ndarray,
    indices: np.ndarray,
    charges: np.ndarray,
    cum_avgs: np.ndarray,
) -> None:
    """Plot per-atom charge time evolution with cumulative averages."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    colors = plt.cm.tab10(np.linspace(0, 1, min(len(indices), 10)))

    for j, idx in enumerate(indices):
        c = colors[j % len(colors)]
        ax.plot(steps, charges[:, j], lw=0.8, alpha=0.5, color=c,
                label=f"atom {int(idx)} (inst.)")
        ax.plot(steps, cum_avgs[:, j], lw=2.0, ls="--", color=c,
                label=f"atom {int(idx)} (cum. avg)")

    ax.set_xlabel("MD step")
    ax.set_ylabel("Bader net charge (e)")
    ax.set_title("Tracked atom charges")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Counterion charge analysis (NEW — per-frame detection, XYZ order)
# ---------------------------------------------------------------------------

def counterion_charge_analysis(
    root_dir: str | Path = ".",
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
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
    """Per-frame counterion detection with charge time evolution + ensemble average.

    Detects non-water, non-metal atoms with non-zero net charge at each frame,
    remaps their indices back to original XYZ ordering, and outputs:

    - Per-frame CSV with detected atom indices and charges.
    - Summary CSV with per-atom mean charge and detection frequency.
    - PNG with charge time evolution for each unique detected atom.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    metal_symbols
        Override default metal symbols for layer/water detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).
    layer_tol_A
        Layer clustering tolerance in Ångströms.
    dir_pattern
        Glob pattern for frame subdirectories.
    structure_filename, acf_filename, potcar_filename
        Per-frame file names.
    output_dir
        Where to write CSV and PNG.  Defaults to *root_dir*.
    frame_start, frame_end, frame_step
        Slice parameters applied to the sorted frame list.
    verbose
        Show tqdm progress bar.

    Returns
    -------
    Path
        Path to the written per-frame CSV file.
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

    if output_dir is None:
        output_dir = root
    output_dir = Path(output_dir)

    logger.info("Counterion charge analysis: %d frames", len(frame_dirs))

    axis_idx = AXIS_MAP[normal]

    # Collect per-frame results: list of (step, time, {xyz_idx: net_charge})
    frame_records: list[tuple[int, int, dict[int, float]]] = []
    all_xyz_indices: set[int] = set()

    iterator: enumerate = enumerate(frame_dirs)
    if verbose:
        from tqdm import tqdm
        iterator = enumerate(tqdm(frame_dirs, desc="Counterion detect", unit="frame", ascii=" ="))

    for _idx, frame_dir in iterator:
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
        imap = read_index_map_from_poscar(poscar)
        net_charge = atoms.arrays["bader_net_charge"]

        # Detect interface layers
        det = detect_interface_layers(
            atoms, metal_symbols=metal_symbols, normal=normal,
            layer_tol_A=layer_tol_A,
        )

        # Detect water molecules
        water_mol = detect_water_molecule_indices(atoms)
        exclude = set(water_mol.ravel().tolist()) | set(det.metal_indices)

        # Filter: non-zero net charge AND not water/metal
        charged_poscar = np.array(
            [i for i in np.where(net_charge != 0.0)[0] if i not in exclude],
            dtype=int,
        )

        step, time_fs = _extract_step_and_time(fname)

        if charged_poscar.size == 0:
            frame_records.append((step, time_fs, {}))
            continue

        # Remap POSCAR indices → XYZ indices
        charged_xyz = imap.poscar_to_xyz[charged_poscar].astype(int)
        charges = net_charge[charged_poscar]

        record: dict[int, float] = {}
        for poscar_i, xyz_i, q in zip(charged_poscar, charged_xyz, charges):
            record[int(xyz_i)] = float(q)

        all_xyz_indices.update(record.keys())
        frame_records.append((step, time_fs, record))

    # --- Per-frame CSV ---
    sorted_xyz_indices = sorted(all_xyz_indices)
    idx_q_fields = []
    for xi in sorted_xyz_indices:
        idx_q_fields.append(f"atom{xi}_q")
    fieldnames = ["step", "time_fs"] + idx_q_fields

    csv_rows: list[dict] = []
    for step, time_fs, record in frame_records:
        row: dict = {"step": step, "time_fs": time_fs}
        for xi in sorted_xyz_indices:
            row[f"atom{xi}_q"] = float(record[xi]) if xi in record else ""
        csv_rows.append(row)

    csv_path = output_dir / DEFAULT_COUNTERION_CHARGE_CSV
    _write_csv(csv_path, csv_rows, fieldnames)

    # --- Summary CSV ---
    summary_fields = ["atom_idx_xyz", "mean_charge_e", "n_frames_detected", "detection_ratio"]
    summary_rows: list[dict] = []
    n_frames = len(frame_records)
    for xi in sorted_xyz_indices:
        charges_list = [rec[xi] for _, _, rec in frame_records if xi in rec]
        n_detected = len(charges_list)
        mean_q = float(np.mean(charges_list)) if charges_list else 0.0
        summary_rows.append({
            "atom_idx_xyz": xi,
            "mean_charge_e": mean_q,
            "n_frames_detected": n_detected,
            "detection_ratio": round(n_detected / n_frames, 4) if n_frames > 0 else 0.0,
        })
    summary_path = output_dir / DEFAULT_COUNTERION_SUMMARY_CSV
    _write_csv(summary_path, summary_rows, summary_fields)

    # --- PNG ---
    png_path = output_dir / DEFAULT_COUNTERION_CHARGE_PNG
    _plot_counterion_charges(png_path, frame_records, sorted_xyz_indices)

    logger.info("Counterion charge analysis: %d frames, %d unique atoms → %s",
                n_frames, len(sorted_xyz_indices), csv_path)
    return csv_path


def _plot_counterion_charges(
    png_path: Path,
    frame_records: list[tuple[int, int, dict[int, float]]],
    sorted_xyz_indices: list[int],
) -> None:
    """Plot per-atom counterion charge time evolution."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not sorted_xyz_indices:
        # Nothing to plot
        return

    steps = np.array([rec[0] for rec in frame_records])
    n_frames = len(steps)
    n_atoms = len(sorted_xyz_indices)

    # Build charge matrix (n_frames, n_atoms), NaN where not detected
    charge_matrix = np.full((n_frames, n_atoms), np.nan)
    for i, (_, _, record) in enumerate(frame_records):
        for j, xi in enumerate(sorted_xyz_indices):
            if xi in record:
                charge_matrix[i, j] = record[xi]

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    colors = plt.cm.tab10(np.linspace(0, 1, min(n_atoms, 10)))

    for j, xi in enumerate(sorted_xyz_indices):
        c = colors[j % len(colors)]
        series = charge_matrix[:, j]
        mask = ~np.isnan(series)
        if mask.any():
            ax.plot(steps[mask], series[mask], lw=0.8, alpha=0.6, color=c,
                    label=f"atom {xi}", marker=".", markersize=2)

    ax.set_xlabel("MD step")
    ax.set_ylabel("Bader net charge (e)")
    ax.set_title("Counterion charges (per-frame detection)")
    ax.grid(True, alpha=0.25)
    if n_atoms <= 20:
        ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)
