"""Calibration data and fit visualisation."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ._mapper import ChargePotentialMapper, FittingInfo
from .config import DEFAULT_CALIBRATION_PNG_NAME

logger = logging.getLogger(__name__)


def plot_calibration(
    sigma: np.ndarray,
    phi: np.ndarray,
    mapper: ChargePotentialMapper,
    fitting_info: FittingInfo,
    *,
    reference: str = "SHE",
    output_dir: Path | None = None,
) -> Path:
    """Plot calibration scatter + fitted curve with quality annotation.

    Style follows ``PhiZProfile.py``: ``figsize=(11, 4.8)``, ``dpi=160``,
    grid enabled, ``tight_layout``.

    Parameters
    ----------
    sigma, phi : (n,) arrays
        Calibration data points.
    mapper : ChargePotentialMapper
        Fitted mapper for drawing the curve.
    fitting_info : FittingInfo
        R² and equation string for annotation.
    reference : str
        Potential reference label for axis label.
    output_dir : Path, optional
        Directory for the PNG. Defaults to cwd.

    Returns
    -------
    Path
        Path to the saved PNG file.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 4.8), dpi=160)

    # Data points
    ax.scatter(
        sigma, phi,
        c="#1f77b4", s=60, zorder=5,
        edgecolors="k", linewidths=0.5,
        label="Calibration data",
    )

    # Fitted curve (smooth)
    margin = 0.05 * (sigma.max() - sigma.min()) if sigma.max() != sigma.min() else 1.0
    sigma_fine = np.linspace(sigma.min() - margin, sigma.max() + margin, 200)
    phi_fit = mapper.predict(sigma_fine)
    ax.plot(
        sigma_fine, phi_fit,
        color="#d62728", lw=2.0,
        label=f"Fit ({fitting_info.method})",
    )

    # Annotation box: equation + R² + RMSE
    lines = [fitting_info.equation_str]
    lines.append(f"R² = {fitting_info.r_squared:.6f}")
    lines.append(f"RMSE = {fitting_info.rmse:.4e}")
    annotation = "\n".join(lines)
    ax.annotate(
        annotation,
        xy=(0.05, 0.95), xycoords="axes fraction",
        va="top", fontsize=10, fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.3", fc="wheat", alpha=0.5),
    )

    # Extra annotations for differential_capacitance: draw interval boundaries
    # and label each segment with its capacitance value.
    if (
        fitting_info.method == "differential_capacitance"
        and getattr(mapper, "_sigma_nodes", None) is not None
    ):
        sigma_nodes = mapper._sigma_nodes
        phi_nodes = mapper._phi_nodes
        caps = mapper._capacitances
        for j in range(len(caps)):
            sigma_mid = 0.5 * (sigma_nodes[j] + sigma_nodes[j + 1])
            phi_mid = 0.5 * (phi_nodes[j] + phi_nodes[j + 1])
            ax.annotate(
                f"C={caps[j]:.1f} μF/cm²",
                xy=(sigma_mid, phi_mid),
                fontsize=8, ha="center", va="center", color="#2ca02c",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, ec="none"),
            )
        for s_node in sigma_nodes[1:-1]:
            ax.axvline(s_node, color="gray", ls="--", lw=0.8, alpha=0.5)

    ax.set_xlabel("σ (μC/cm²)")
    ax.set_ylabel(f"φ (V vs {reference})")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", frameon=True)

    fig.tight_layout()

    out_dir = Path(output_dir) if output_dir is not None else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / DEFAULT_CALIBRATION_PNG_NAME
    fig.savefig(out_png)
    plt.close(fig)

    logger.info("Calibration plot saved to %s", out_png)
    return out_png
