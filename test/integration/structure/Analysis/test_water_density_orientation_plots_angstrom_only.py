"""Integration test: generate density/orientation plots in Angstrom domain only."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import MultipleLocator

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.structure.Analysis import water_mass_density_z_distribution_analysis
from scripts.structure.Analysis import water_orientation_weighted_density_z_distribution_analysis


def test_water_density_and_orientation_plots_angstrom_only() -> None:
    """Generate one combined figure with Angstrom x-axis only."""
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    density_csv_path = water_mass_density_z_distribution_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
        output_csv_name="water_density_interface_to_midpoint_from_potential.csv",
    )
    orientation_csv_path = water_orientation_weighted_density_z_distribution_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
        output_csv_name="water_orientation_interface_to_midpoint_from_potential.csv",
    )

    density_data = np.loadtxt(density_csv_path, delimiter=",", skiprows=1)
    if density_data.ndim == 1:
        density_data = density_data.reshape(1, -1)
    distance_A_density = density_data[:, 1]
    rho_g_cm3 = density_data[:, 2]

    orientation_data = np.loadtxt(orientation_csv_path, delimiter=",", skiprows=1)
    if orientation_data.ndim == 1:
        orientation_data = orientation_data.reshape(1, -1)
    distance_A_orientation = orientation_data[:, 1]
    orient_1_A3 = orientation_data[:, 2]

    assert np.all(np.isfinite(rho_g_cm3))
    assert np.all(np.isfinite(orient_1_A3))
    assert distance_A_density.size > 0
    assert distance_A_orientation.size > 0

    combined_png_path = out_dir / "water_density_orientation_interface_to_midpoint_angstrom_only.png"
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # Top panel: water mass density (major tick interval = 1).
    ax1.plot(distance_A_density, rho_g_cm3, color="tab:blue", lw=1.5)
    ax1.set_ylabel("rho_ensemble_avg (g/cm^3)")
    ax1.set_title("Water Density: interface -> midpoint")
    ax1.yaxis.set_major_locator(MultipleLocator(1.0))
    ax1.grid(True, ls="--", alpha=0.3)

    # Bottom panel: orientation profile (major ticks = 3, unit = a.u.).
    ax2.plot(distance_A_orientation, orient_1_A3, color="tab:orange", lw=1.5)
    ax2.set_xlabel("distance from selected interface (Angstrom)")
    ax2.set_ylabel("orientation_ensemble_avg (a.u.)")
    ax2.set_title("Water Orientation: interface -> midpoint")
    ax2.yaxis.set_major_locator(LinearLocator(3))
    ax2.grid(True, ls="--", alpha=0.3)

    # Align x-axis range for both panels.
    x_min = float(min(np.min(distance_A_density), np.min(distance_A_orientation)))
    x_max = float(max(np.max(distance_A_density), np.max(distance_A_orientation)))
    ax1.set_xlim(x_min, x_max)

    fig.tight_layout()
    fig.savefig(combined_png_path, dpi=180)
    plt.close(fig)

    assert combined_png_path.exists() and combined_png_path.stat().st_size > 0


if __name__ == "__main__":
    test_water_density_and_orientation_plots_angstrom_only()
    print("saved: test/_tmp_preview/water_density_orientation_interface_to_midpoint_angstrom_only.png")
