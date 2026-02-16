"""Integration test: run WaterDensity analysis on potential data and plot output."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.structure.Analysis import water_mass_density_z_distribution_analysis


def test_water_density_analysis_and_plot_from_potential_data() -> None:
    """Generate CSV + PNG from data_example/potential using WaterDensity analysis."""
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = water_mass_density_z_distribution_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
        output_csv_name="water_density_interface_to_midpoint_from_potential.csv",
    )

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    path_fraction_center = data[:, 0]
    distance_A = data[:, 1]
    rho_g_cm3 = data[:, 2]

    assert path_fraction_center.size > 0
    assert np.all(np.isfinite(rho_g_cm3))

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=False)
    axes[0].plot(distance_A, rho_g_cm3, color="tab:blue", lw=1.5)
    axes[0].set_xlabel("distance from selected interface (Angstrom)")
    axes[0].set_ylabel("rho_ensemble_avg (g/cm^3)")
    axes[0].set_title("Water Density: interface -> midpoint (distance domain)")
    axes[0].grid(True, ls="--", alpha=0.3)

    axes[1].plot(path_fraction_center, rho_g_cm3, color="tab:green", lw=1.5)
    axes[1].set_xlabel("path fraction (0 to 1)")
    axes[1].set_ylabel("rho_ensemble_avg (g/cm^3)")
    axes[1].set_title("Water Density: interface -> midpoint (fraction domain)")
    axes[1].grid(True, ls="--", alpha=0.3)

    fig.tight_layout()
    png_path = out_dir / "water_density_interface_to_midpoint_from_potential.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    assert csv_path.exists()
    assert csv_path.stat().st_size > 0
    assert png_path.exists()
    assert png_path.stat().st_size > 0


if __name__ == "__main__":
    test_water_density_analysis_and_plot_from_potential_data()
    print("saved: test/_tmp_preview/water_density_interface_to_midpoint_from_potential.csv")
    print("saved: test/_tmp_preview/water_density_interface_to_midpoint_from_potential.png")
