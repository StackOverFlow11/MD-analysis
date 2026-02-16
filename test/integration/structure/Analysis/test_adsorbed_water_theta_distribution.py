"""Integration test: theta distribution (0-180 deg) within adsorbed-water layer."""

from __future__ import annotations

from pathlib import Path
import re
import sys

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.structure.Analysis import ad_water_orientation_analysis
from scripts.structure.Analysis import compute_adsorbed_water_theta_distribution
from scripts.structure.Analysis import DEFAULT_START_INTERFACE


def _parse_range_txt(range_txt_path: Path) -> tuple[float, float]:
    text = range_txt_path.read_text(encoding="utf-8")
    m_start = re.search(r"adsorbed_layer_start_A=([0-9.eE+-]+)", text)
    m_end = re.search(r"adsorbed_layer_end_A=([0-9.eE+-]+)", text)
    if m_start is None or m_end is None:
        raise ValueError(f"Cannot parse adsorbed range from {range_txt_path}")
    return float(m_start.group(1)), float(m_end.group(1))


def test_adsorbed_water_theta_distribution_0_180_plot() -> None:
    """Compute and plot theta PDF (0-180 degree) inside detected adsorbed-water range."""
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: detect adsorbed-water distance range from density profile.
    _, range_txt_path = ad_water_orientation_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
    )
    d_start_A, d_end_A = _parse_range_txt(range_txt_path)
    assert d_end_A > d_start_A

    # Step 2: compute theta distribution in detected adsorbed range.
    theta_centers, theta_pdf, _ = compute_adsorbed_water_theta_distribution(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        adsorbed_range_A=(d_start_A, d_end_A),
        output_dir=out_dir,
        start_interface=DEFAULT_START_INTERFACE,
    )
    assert theta_centers.size > 0
    assert np.any(theta_pdf > 0.0)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4.8))
    ax.plot(theta_centers, theta_pdf, color="tab:green", lw=1.6)
    ax.set_xlim(0.0, 180.0)
    ax.set_xlabel("theta (degree)")
    ax.set_ylabel("PDF (degree^-1)")
    ax.set_title(
        f"Adsorbed-Water Theta Distribution (0-180 deg), range=[{d_start_A:.2f}, {d_end_A:.2f}] A"
    )
    ax.grid(True, ls="--", alpha=0.3)

    fig.tight_layout()
    png_path = out_dir / "adsorbed_water_theta_distribution_0_180.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    assert png_path.exists() and png_path.stat().st_size > 0


if __name__ == "__main__":
    test_adsorbed_water_theta_distribution_0_180_plot()
    print("saved: test/_tmp_preview/adsorbed_water_theta_distribution_0_180.png")
