"""Matplotlib plotting helpers for Bader charge analysis.

Kept separate from the computation modules (``SurfaceCharge.py``,
``AtomCharges.py``) so that calibration / charge calculations can be tested
without a matplotlib dependency.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def plot_surface_charge(
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


def plot_single_side_charge(
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


def plot_tracked_charges(
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


def plot_counterion_charges(
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
