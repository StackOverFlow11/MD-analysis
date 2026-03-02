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


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _moving_slab_average_z_periodic(
    phi_z: np.ndarray,
    z_ang: np.ndarray,
    window_ang: float,
) -> np.ndarray:
    """Moving average along z with periodic boundary conditions.

    For each z position, averages φ within ±window_ang/2 using periodic distance.
    Assumes a uniform grid (typical for CP2K cube files).
    """
    phi_z = np.asarray(phi_z, dtype=float)
    z_ang = np.asarray(z_ang, dtype=float)
    if phi_z.ndim != 1 or z_ang.ndim != 1 or phi_z.shape != z_ang.shape:
        raise ValueError("phi_z and z_ang must be 1D arrays of the same shape")
    if window_ang <= 0:
        return phi_z.copy()

    dz = float(np.mean(np.diff(z_ang)))
    if not np.isfinite(dz) or dz <= 0:
        raise ValueError(f"Invalid dz from z grid: {dz}")
    lz = dz * float(z_ang.size)

    z_rel = (z_ang - float(z_ang[0])) % lz
    half = 0.5 * float(window_ang)

    out = np.empty_like(phi_z, dtype=float)
    for i, zi in enumerate(z_rel):
        dist = np.abs(z_rel - zi)
        dist = np.minimum(dist, lz - dist)
        m = dist <= half
        out[i] = float(phi_z[m].mean())
    return out


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
    z_mavg_window_ang: float = 7.0,
    max_curves: int = 0,
) -> Path:
    """Full-frame φ(z) plane-averaged potential profile analysis.

    Reads all cube files matching *cube_pattern*, computes the
    xy plane-average at each z-slice, and produces:

    - ``phi_z_planeavg_stats.csv`` (mean/std/min/max over frames)
    - ``phi_z_planeavg_all_frames.png`` (heatmap + overlay plot)
    - ``phi_z_planeavg_mavg_{window}A_stats.csv`` (if window > 0)

    Returns the PNG path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = [Path(p) for p in sorted(workdir.glob(cube_pattern))]
    if not cube_paths:
        raise FileNotFoundError(f"No cube files matched pattern: {cube_pattern!r} in {workdir}")

    steps: list[int] = []
    phi_list: list[np.ndarray] = []
    z_ang_ref: Optional[np.ndarray] = None

    for cp in cube_paths:
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

    # Moving average along z
    z_mavg = float(z_mavg_window_ang)
    phi_mat_mavg = np.vstack([_moving_slab_average_z_periodic(row, z_ang_ref, z_mavg) for row in phi_mat])
    phi_mean_mavg = phi_mat_mavg.mean(axis=0)
    phi_std_mavg = phi_mat_mavg.std(axis=0, ddof=1) if phi_mat_mavg.shape[0] > 1 else np.zeros_like(phi_mean_mavg)
    phi_min_mavg = phi_mat_mavg.min(axis=0)
    phi_max_mavg = phi_mat_mavg.max(axis=0)
    _write_phi_z_csv(
        outdir / f"phi_z_planeavg_mavg_{z_mavg:g}A_stats.csv",
        z_ang_ref, phi_mean_mavg, phi_std_mavg, phi_min_mavg, phi_max_mavg,
    )

    # --- Plotting ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11, 8.5), dpi=160)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1.0], hspace=0.22)

    # Heatmap panel
    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(
        phi_mat,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        extent=[float(z_ang_ref[0]), float(z_ang_ref[-1]), float(steps_arr[0]), float(steps_arr[-1])],
        cmap="RdBu_r",
    )
    ax0.set_title("Plane-averaged Hartree potential φ(z) over snapshots")
    ax0.set_xlabel("z (Å)")
    ax0.set_ylabel("MD step (cube snapshots)")
    cbar = fig.colorbar(im, ax=ax0, pad=0.012)
    cbar.set_label("φ(z) (eV)  [plane-avg Hartree potential]")

    # Overlay panel
    ax1 = fig.add_subplot(gs[1, 0])
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

    if z_mavg > 0:
        ax1.fill_between(
            z_ang_ref,
            phi_mean_mavg - phi_std_mavg,
            phi_mean_mavg + phi_std_mavg,
            color="#d62728",
            alpha=0.12,
            label=f"moving avg ({z_mavg:g} Å) mean ± 1σ",
        )
        ax1.plot(z_ang_ref, phi_mean_mavg, color="#d62728", lw=2.2, label=f"moving avg ({z_mavg:g} Å) mean")

    ax1.set_xlabel("z (Å)")
    ax1.set_ylabel("φ(z) (eV)")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", frameon=True)

    fig.tight_layout()
    out_png = outdir / DEFAULT_PHI_Z_PNG_NAME
    fig.savefig(out_png)
    plt.close(fig)

    return out_png
