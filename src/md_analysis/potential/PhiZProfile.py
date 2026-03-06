"""Full-frame φ(z) plane-averaged potential profile analysis and visualization.

Public API
----------
- ``phi_z_planeavg_analysis``
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import numpy as np

from ..utils.CubeParser import (
    extract_step_from_cube_filename,
    plane_avg_phi_z_ev,
    read_cube_header_and_values,
    z_coords_ang,
)
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
# Public API
# ---------------------------------------------------------------------------

def phi_z_planeavg_analysis(
    cube_pattern: str,
    *,
    output_dir: Path | None = None,
    max_curves: int = 0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """Full-frame φ(z) plane-averaged potential profile analysis.

    Reads all cube files matching *cube_pattern*, computes the
    xy plane-average at each z-slice, and produces:

    - ``phi_z_planeavg_stats.csv`` (mean/std/min/max over frames)
    - ``phi_z_planeavg_all_frames.png`` (overlay plot)

    Returns the PNG path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = [Path(p) for p in sorted(workdir.glob(cube_pattern))]
    if not cube_paths:
        raise FileNotFoundError(f"No cube files matched pattern: {cube_pattern!r} in {workdir}")
    cube_paths = cube_paths[frame_start:frame_end:frame_step]

    steps: list[int] = []
    phi_list: list[np.ndarray] = []
    z_ang_ref: Optional[np.ndarray] = None

    cube_iter = cube_paths
    if verbose:
        from tqdm import tqdm
        cube_iter = tqdm(cube_paths, desc="φ(z) plane-avg", unit="cube", ascii=" =")

    for cp in cube_iter:
        header, values = read_cube_header_and_values(cp)
        z = z_coords_ang(header)
        phi_z_ev = plane_avg_phi_z_ev(header, values)

        if z_ang_ref is None:
            z_ang_ref = z
        else:
            if z.shape != z_ang_ref.shape or not np.allclose(z, z_ang_ref, rtol=0, atol=1e-6):
                phi_z_ev = np.interp(z_ang_ref, z, phi_z_ev)

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
