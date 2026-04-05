"""Matplotlib plotting helpers for potential analysis.

Kept separate from the computation modules (``CenterPotential.py``,
``PhiZProfile.py``) so that numerical analyses can be tested without a
matplotlib dependency.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def plot_series_with_cumavg(
    png_path: Path,
    x: np.ndarray,
    y: np.ndarray,
    y_cum: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    """Plot instantaneous values with cumulative average overlay."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(x, y, lw=1.0, alpha=0.65, label="instantaneous")
    ax.plot(x, y_cum, lw=2.0, label="cumulative average")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def plot_thickness_sensitivity(
    png_path: Path,
    thicknesses: np.ndarray,
    means: np.ndarray,
    spatial_stds: np.ndarray,
) -> None:
    """Dual-axis plot: mean U vs SHE (left) and spatial std φ(z) (right)
    as a function of slab averaging thickness."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(9, 4.8), dpi=160)
    color_mean = "tab:blue"
    ax1.plot(thicknesses, means, "o-", color=color_mean, lw=1.5, markersize=4)
    ax1.set_xlabel("Slab thickness (Å)")
    ax1.set_ylabel("Mean U vs SHE (V)", color=color_mean)
    ax1.tick_params(axis="y", labelcolor=color_mean)

    ax2 = ax1.twinx()
    color_std = "tab:red"
    ax2.plot(thicknesses, spatial_stds, "s--", color=color_std, lw=1.5, markersize=4)
    ax2.set_ylabel("Spatial std of φ(z) in slab (eV)", color=color_std)
    ax2.tick_params(axis="y", labelcolor=color_std)

    ax1.set_title("Electrode potential U vs SHE — thickness sensitivity")
    ax1.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def plot_phi_z_profile(
    png_path: Path,
    z_ang: np.ndarray,
    phi_mat: np.ndarray,
    phi_mean: np.ndarray,
    phi_std: np.ndarray,
    phi_min: np.ndarray,
    phi_max: np.ndarray,
    *,
    max_curves: int = 0,
) -> None:
    """Plot plane-averaged φ(z) profiles: raw per-frame overlay + mean ± std envelope.

    When ``max_curves > 0`` and more than *max_curves* frames are available,
    a deterministic random subset (seed 0) is shown.
    """
    png_path.parent.mkdir(parents=True, exist_ok=True)

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
        ax1.plot(z_ang, row, color="#1f77b4", alpha=0.05, lw=0.8)

    ax1.fill_between(z_ang, phi_mean - phi_std, phi_mean + phi_std,
                     color="k", alpha=0.15, label="mean ± 1σ")
    ax1.plot(z_ang, phi_mean, color="k", lw=2.0, label=f"mean {label_suffix}")
    ax1.plot(z_ang, phi_min, color="k", lw=1.0, ls="--", alpha=0.45,
             label="min/max envelope")
    ax1.plot(z_ang, phi_max, color="k", lw=1.0, ls="--", alpha=0.45)

    ax1.set_xlabel("z (Å)")
    ax1.set_ylabel("φ(z) (eV)")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", frameon=True)

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)
