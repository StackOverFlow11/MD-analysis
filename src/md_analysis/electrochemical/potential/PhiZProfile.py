"""Full-frame φ(z) plane-averaged potential profile analysis and visualization.

Public API
----------
- ``phi_z_planeavg_analysis``
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

from ...utils.config import BOHR_TO_ANG, DEFAULT_LAYER_TOL_A, TRANSITION_METAL_SYMBOLS
from ...utils.CubeParser import (
    CubeHeader,
    _float,
    discover_cube_files,
    extract_step_from_cube_filename,
    plane_avg_phi_z_ev,
    read_cube_header_and_values,
    z_coords_ang,
)
from ...utils.StructureParser.ClusterUtils import gap_midpoint_periodic
from ...utils.StructureParser.LayerParser import detect_interface_layers
from .config import DEFAULT_PHI_Z_PNG_NAME, DEFAULT_PHI_Z_STATS_CSV_NAME


def _write_phi_z_csv(
    path: Path,
    z_ang: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    vmin: np.ndarray,
    vmax: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["z_ang", "phi_mean_ev", "phi_std_ev", "phi_min_ev", "phi_max_ev"])
        for zz, m, s, mn, mx in zip(z_ang, mean, std, vmin, vmax, strict=True):
            w.writerow([float(zz), float(m), float(s), float(mn), float(mx)])


# ---------------------------------------------------------------------------
# Private helpers — slab centering
# ---------------------------------------------------------------------------

def _read_cube_atoms(path: Path, header: CubeHeader):
    """Parse atomic positions from a cube file header and return ``ase.Atoms``.

    Reuses the already-parsed *header* for cell vectors; only re-reads the
    atom lines from *path*.
    """
    from ase import Atoms
    from ase.data import chemical_symbols

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for _ in range(6):          # 2 comments + natoms/origin + 3 voxel lines
            f.readline()
        symbols: list[str] = []
        positions_ang: list[list[float]] = []
        for _ in range(header.natoms):
            parts = f.readline().split()
            z_num = int(parts[0])
            x, y, z = _float(parts[2]), _float(parts[3]), _float(parts[4])
            symbols.append(chemical_symbols[z_num])
            positions_ang.append([x * BOHR_TO_ANG, y * BOHR_TO_ANG, z * BOHR_TO_ANG])

    cell_ang = np.array([
        header.vx_bohr * header.nx * BOHR_TO_ANG,
        header.vy_bohr * header.ny * BOHR_TO_ANG,
        header.vz_bohr * header.nz * BOHR_TO_ANG,
    ])
    return Atoms(symbols=symbols, positions=positions_ang, cell=cell_ang, pbc=True)


def _slab_center_roll(
    atoms,
    nz: int,
    *,
    metal_symbols: set[str],
    layer_tol_A: float,
) -> int:
    """Return the ``np.roll`` shift (grid points) that centres the slab at ``nz/2``."""
    detection = detect_interface_layers(
        atoms, metal_symbols=metal_symbols, normal="c", layer_tol_A=layer_tol_A,
    )
    aligned = detection.interface_normal_aligned()
    opposed = detection.interface_normal_opposed()
    slab_center_frac = gap_midpoint_periodic(
        aligned.center_frac, opposed.center_frac, 1.0,
    )
    shift_frac = 0.5 - slab_center_frac
    # Wrap to (-0.5, 0.5] then convert to grid points
    shift_frac = (shift_frac + 0.5) % 1.0 - 0.5
    return round(shift_frac * nz)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def phi_z_planeavg_analysis(
    cube_pattern: str,
    *,
    output_dir: Path | None = None,
    max_curves: int = 0,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """Full-frame φ(z) plane-averaged potential profile analysis.

    Reads all cube files matching *cube_pattern*, computes the
    xy plane-average at each z-slice, and produces:

    - ``phi_z_planeavg_stats.csv`` (mean/std/min/max over frames)
    - ``phi_z_planeavg_all_frames.png`` (overlay plot, slab centred at cell_c/2)

    Each frame's φ(z) is cyclically shifted so that the metal-slab
    midpoint sits at cell_c / 2 (the x-axis midpoint of the plot).

    Returns the PNG path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = discover_cube_files(
        cube_pattern, workdir=workdir,
        frame_start=frame_start, frame_end=frame_end, frame_step=frame_step,
    )
    logger.info("phi(z) analysis: %d cube files", len(cube_paths))

    steps: list[int] = []
    phi_list: list[np.ndarray] = []
    z_ang_ref: Optional[np.ndarray] = None
    metal_symbols: Optional[set[str]] = (
        set(metal_elements) if metal_elements is not None else None
    )

    cube_iter = cube_paths
    if verbose:
        from tqdm import tqdm
        cube_iter = tqdm(cube_paths, desc="φ(z) plane-avg", unit="cube", ascii=" =")

    shift_n: int | None = None          # computed once from the first frame

    for cp in cube_iter:
        header, values = read_cube_header_and_values(cp)
        z = z_coords_ang(header)
        phi_z_ev = plane_avg_phi_z_ev(header, values)

        if z_ang_ref is None:
            z_ang_ref = z
        else:
            if z.shape != z_ang_ref.shape or not np.allclose(z, z_ang_ref, rtol=0, atol=1e-6):
                phi_z_ev = np.interp(z_ang_ref, z, phi_z_ev)

        # --- slab centering: roll φ(z) so slab midpoint → cell_c/2 ---
        # Use the first frame's shift for all frames to keep profiles aligned.
        if shift_n is None:
            atoms = _read_cube_atoms(cp, header)
            if metal_symbols is None:
                metal_symbols = set(atoms.get_chemical_symbols()) & set(TRANSITION_METAL_SYMBOLS)
                if metal_symbols:
                    logger.info("Auto-detected metal elements: %s", metal_symbols)
            if metal_symbols:
                shift_n = _slab_center_roll(
                    atoms, header.nz,
                    metal_symbols=metal_symbols, layer_tol_A=layer_tol_ang,
                )
            else:
                shift_n = 0
        if shift_n:
            phi_z_ev = np.roll(phi_z_ev, shift_n)

        s = extract_step_from_cube_filename(cp)
        steps.append(int(s) if s is not None else 0)
        phi_list.append(phi_z_ev)

    assert z_ang_ref is not None
    steps_arr = np.array(steps, dtype=int)
    phi_mat = np.vstack(phi_list)

    # Sort by step
    order = np.argsort(steps_arr)
    steps_arr = steps_arr[order]
    phi_mat = phi_mat[order, :]

    phi_mean = phi_mat.mean(axis=0)
    phi_std = phi_mat.std(axis=0, ddof=1) if phi_mat.shape[0] > 1 else np.zeros_like(phi_mean)
    phi_min = phi_mat.min(axis=0)
    phi_max = phi_mat.max(axis=0)

    _write_phi_z_csv(outdir / DEFAULT_PHI_Z_STATS_CSV_NAME, z_ang_ref, phi_mean, phi_std, phi_min, phi_max)

    # --- Plotting ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(11, 4.8), dpi=160)

    if max_curves > 0 and phi_mat.shape[0] > max_curves:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(phi_mat.shape[0], size=max_curves, replace=False))
        curves = phi_mat[idx]
        label_suffix = f"(random {max_curves}/{phi_mat.shape[0]})"
    else:
        curves = phi_mat
        label_suffix = f"({phi_mat.shape[0]} frames)"

    for row in curves:
        ax1.plot(z_ang_ref, row, color="#1f77b4", alpha=0.05, lw=0.8)

    ax1.fill_between(z_ang_ref, phi_mean - phi_std, phi_mean + phi_std, color="k", alpha=0.15, label="mean ± 1σ")
    ax1.plot(z_ang_ref, phi_mean, color="k", lw=2.0, label=f"mean {label_suffix}")
    ax1.plot(z_ang_ref, phi_min, color="k", lw=1.0, ls="--", alpha=0.45, label="min/max envelope")
    ax1.plot(z_ang_ref, phi_max, color="k", lw=1.0, ls="--", alpha=0.45)

    ax1.set_xlabel("z (Å)")
    ax1.set_ylabel("φ(z) (eV)")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", frameon=True)

    fig.tight_layout()
    out_png = outdir / DEFAULT_PHI_Z_PNG_NAME
    fig.savefig(out_png)
    plt.close(fig)

    return out_png
