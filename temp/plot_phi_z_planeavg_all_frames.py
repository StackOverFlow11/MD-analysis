"""
Quick helper: plot plane-averaged Hartree potential profiles φ(z) for all CP2K cube snapshots.

It reads md-POTENTIAL-v_hartree-1_*.cube (or a user-provided glob), computes the xy plane
average at each z-slice, and then:
  - plots a heatmap (MD step vs z, colored by φ(z))
  - overlays all φ(z) curves + mean ± std on a second panel

Outputs are written into ./analysis_outputs by default.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np

import analyze_potential as ap


def _extract_step(path: Path) -> int:
    s = ap._extract_step_from_cube_filename(path)
    if s is None:
        # stable ordering fallback
        return 0
    return int(s)


def _plane_avg_phi_z_ev(header: ap.CubeHeader, values: np.ndarray) -> np.ndarray:
    # CP2K cube ordering: z fastest, then y, then x -> reshape (nx, ny, nz)
    field = values.reshape((header.nx, header.ny, header.nz))
    phi_z_ha = field.mean(axis=(0, 1))
    return phi_z_ha * ap.HA_TO_EV


def _z_coords_ang(header: ap.CubeHeader) -> np.ndarray:
    dz_bohr = float(np.linalg.norm(header.vz_bohr))
    dz_ang = dz_bohr * ap.BOHR_TO_ANG
    origin_z_ang = float(header.origin_bohr[2]) * ap.BOHR_TO_ANG
    # use slice centers (consistent with slab-average routine)
    return origin_z_ang + (np.arange(header.nz, dtype=float) + 0.5) * dz_ang


def moving_slab_average_z_periodic(phi_z: np.ndarray, z_ang: np.ndarray, window_ang: float) -> np.ndarray:
    """
    Moving average along z using a physical window size (Angstrom), with PBC on z.

    For each z position, average phi within +/- window_ang/2 using periodic distance.
    """
    phi_z = np.asarray(phi_z, dtype=float)
    z_ang = np.asarray(z_ang, dtype=float)
    if phi_z.ndim != 1 or z_ang.ndim != 1 or phi_z.shape != z_ang.shape:
        raise ValueError("phi_z and z_ang must be 1D arrays of the same shape")
    if window_ang <= 0:
        return phi_z.copy()

    # Uniform grid assumed (CP2K cube): estimate dz and Lz
    dz = float(np.mean(np.diff(z_ang)))
    if not np.isfinite(dz) or dz <= 0:
        raise ValueError(f"Invalid dz from z grid: {dz}")
    lz = dz * float(z_ang.size)

    # Use relative z in [0, Lz) for robust periodic distances
    z_rel = (z_ang - float(z_ang[0])) % lz
    half = 0.5 * float(window_ang)

    out = np.empty_like(phi_z, dtype=float)
    for i, zi in enumerate(z_rel):
        dist = np.abs(z_rel - zi)
        dist = np.minimum(dist, lz - dist)
        m = dist <= half
        out[i] = float(phi_z[m].mean())
    return out


def write_csv(path: Path, z_ang: np.ndarray, mean: np.ndarray, std: np.ndarray, vmin: np.ndarray, vmax: np.ndarray) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["z_ang", "phi_mean_ev", "phi_std_ev", "phi_min_ev", "phi_max_ev"])
        for zz, m, s, mn, mx in zip(z_ang, mean, std, vmin, vmax, strict=True):
            w.writerow([float(zz), float(m), float(s), float(mn), float(mx)])


def main() -> int:
    p = argparse.ArgumentParser(description="Plot plane-averaged Hartree potential φ(z) for all cube snapshots.")
    p.add_argument(
        "--cube-pattern",
        type=str,
        default="md-POTENTIAL-v_hartree-1_*.cube",
        help="Glob pattern for Hartree potential cube snapshots (default: md-POTENTIAL-v_hartree-1_*.cube).",
    )
    p.add_argument(
        "--outdir",
        type=str,
        default="analysis_outputs",
        help="Output directory for PNG/CSV files (default: analysis_outputs).",
    )
    p.add_argument(
        "--max-curves",
        type=int,
        default=0,
        help="If >0, randomly plot at most this many individual φ(z) curves on the overlay panel (default: 0 = plot all).",
    )
    p.add_argument(
        "--z-mavg-window",
        type=float,
        default=7.0,
        help="Moving-average window along z in Angstrom (periodic). Set <=0 to disable (default: 7).",
    )

    args = p.parse_args()
    workdir = Path(".").resolve()
    outdir = (workdir / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = [Path(pth) for pth in sorted(workdir.glob(args.cube_pattern))]
    if not cube_paths:
        raise SystemExit(f"No cube files matched pattern: {args.cube_pattern!r} in {workdir}")

    steps: list[int] = []
    phi_list: list[np.ndarray] = []
    z_ang_ref: Optional[np.ndarray] = None

    for cp in cube_paths:
        header, values = ap.read_cube_header_and_values(cp)
        z_ang = _z_coords_ang(header)
        phi_z_ev = _plane_avg_phi_z_ev(header, values)

        if z_ang_ref is None:
            z_ang_ref = z_ang
        else:
            # Should match in normal CP2K runs (fixed grid). If not, interpolate.
            if z_ang.shape != z_ang_ref.shape or not np.allclose(z_ang, z_ang_ref, rtol=0, atol=1e-6):
                phi_z_ev = np.interp(z_ang_ref, z_ang, phi_z_ev)

        steps.append(_extract_step(cp))
        phi_list.append(phi_z_ev)

    assert z_ang_ref is not None
    steps_arr = np.array(steps, dtype=int)
    phi_mat = np.vstack(phi_list)  # (n_frames, nz)

    # sort by step
    order = np.argsort(steps_arr)
    steps_arr = steps_arr[order]
    phi_mat = phi_mat[order, :]

    phi_mean = phi_mat.mean(axis=0)
    phi_std = phi_mat.std(axis=0, ddof=1) if phi_mat.shape[0] > 1 else np.zeros_like(phi_mean)
    phi_min = phi_mat.min(axis=0)
    phi_max = phi_mat.max(axis=0)

    write_csv(outdir / "phi_z_planeavg_stats.csv", z_ang_ref, phi_mean, phi_std, phi_min, phi_max)

    # ---- moving average along z (per-frame, then stats) ----
    z_mavg = float(args.z_mavg_window)
    phi_mat_mavg = np.vstack([moving_slab_average_z_periodic(row, z_ang_ref, z_mavg) for row in phi_mat])
    phi_mean_mavg = phi_mat_mavg.mean(axis=0)
    phi_std_mavg = phi_mat_mavg.std(axis=0, ddof=1) if phi_mat_mavg.shape[0] > 1 else np.zeros_like(phi_mean_mavg)
    phi_min_mavg = phi_mat_mavg.min(axis=0)
    phi_max_mavg = phi_mat_mavg.max(axis=0)
    write_csv(
        outdir / f"phi_z_planeavg_mavg_{z_mavg:g}A_stats.csv",
        z_ang_ref,
        phi_mean_mavg,
        phi_std_mavg,
        phi_min_mavg,
        phi_max_mavg,
    )

    # ---- plotting ----
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
    max_curves = int(args.max_curves)
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

    # Moving-average profile (7 Å by default)
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
    out_png = outdir / "phi_z_planeavg_all_frames.png"
    fig.savefig(out_png)
    plt.close(fig)

    print("Wrote:")
    print(f"  - {out_png}")
    print(f"  - {outdir / 'phi_z_planeavg_stats.csv'}")
    print(f"  - {outdir / f'phi_z_planeavg_mavg_{z_mavg:g}A_stats.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

