"""Matplotlib plotting helpers for water analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _savgol_smooth_window5(values: np.ndarray) -> np.ndarray:
    """Savitzky-Golay smoothing with fixed 5-point window (polyorder=2).

    Coefficients: [-3, 12, 17, 12, -3] / 35
    """
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size < 5:
        return arr.copy()
    kernel = np.array([-3.0, 12.0, 17.0, 12.0, -3.0], dtype=float) / 35.0
    padded = np.pad(arr, (2, 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def plot_three_panel(
    out_png_path: Path,
    *,
    distance_A: np.ndarray,
    rho_ensemble: np.ndarray,
    orient_ensemble: np.ndarray,
    theta_centers: np.ndarray,
    theta_pdf: np.ndarray,
) -> None:
    """Render the three-panel water analysis figure.

    Panels
    ------
    1. Water mass density vs distance
    2. Orientation-weighted density vs distance
    3. Adsorbed-layer θ distribution PDF (0–180°)

    Input arrays are smoothed with a fixed 5-point Savitzky-Golay kernel
    before plotting; raw data is expected to be persisted separately via
    CSVs in the caller.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator, NullFormatter
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("matplotlib is required to generate plots.") from exc

    rho_smooth = _savgol_smooth_window5(rho_ensemble)
    orient_smooth = _savgol_smooth_window5(orient_ensemble)
    theta_pdf_smooth = _savgol_smooth_window5(theta_pdf)

    with plt.rc_context(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
        }
    ):
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.6, 10.8), sharex=False)

        # Panel 1: density profile.
        ax1.plot(distance_A, rho_smooth, color="tab:blue", lw=1.5)
        ax1.set_ylabel(r"$\rho_\mathrm{{H_2O}}\ \mathrm{(g\cdot cm^{-3})}$", fontsize=18)
        ax1.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax1.set_title("Water density distribution", fontsize=18)
        ax1.tick_params(axis="both", which="major", labelsize=16)

        # Panel 2: orientation-weighted profile.
        ax2.plot(distance_A, orient_smooth, color="tab:orange", lw=1.5)
        ax2.set_ylabel(r"$\rho_\mathrm{{H_2O}}\cdot\cos\varphi\ \mathrm{(g\cdot cm^{-3})}$", fontsize=18)
        ax2.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax2.set_title("Water dipole", fontsize=18)
        ax2.yaxis.set_major_locator(MultipleLocator(1.0))
        ax2.tick_params(axis="both", which="major", labelsize=16)

        x_max = float(np.max(distance_A))
        ax1.set_xlim(0.0, x_max)
        ax2.set_xlim(0.0, x_max)
        ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax2.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax1.xaxis.set_minor_formatter(NullFormatter())
        ax2.xaxis.set_minor_formatter(NullFormatter())

        # Panel 3: adsorbed-layer theta PDF.
        ax3.plot(theta_centers, theta_pdf_smooth, color="tab:green", lw=1.5)
        ax3.set_xlim(0.0, 180.0)
        ax3.set_xlabel(r"$\varphi\ \mathrm{deg}$", fontsize=18)
        ax3.set_ylabel(r"$\mathrm{P(\varphi)\ (a.u.)}$", fontsize=18)
        ax3.set_title("Probability distribution profiles", fontsize=18)
        ax3.xaxis.set_major_locator(MultipleLocator(45.0))
        ax3.xaxis.set_minor_locator(MultipleLocator(22.5))
        ax3.xaxis.set_minor_formatter(NullFormatter())
        ax3.set_yticks([0.0])
        ax3.tick_params(axis="both", which="major", labelsize=16)

        fig.tight_layout()
        fig.savefig(out_png_path, dpi=180)
        plt.close(fig)
