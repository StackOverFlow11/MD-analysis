"""Integrated three-panel water analysis plotting utilities."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import NullFormatter

from .config import DEFAULT_OUTPUT_DIR
from .config import DEFAULT_START_INTERFACE
from .config import DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME
from .WaterAnalysis import ad_water_orientation_analysis
from .WaterAnalysis import compute_adsorbed_water_theta_distribution
from .WaterAnalysis import water_mass_density_z_distribution_analysis
from .WaterAnalysis import water_orientation_weighted_density_z_distribution_analysis
from .WaterAnalysis.WaterDensity import StartInterface
from ..utils.config import DEFAULT_THETA_BIN_DEG
from ..utils.config import DEFAULT_Z_BIN_WIDTH_A


def _savgol_smooth_window5(values: np.ndarray) -> np.ndarray:
    """
    Savitzky-Golay smoothing with fixed 5-point window (polyorder=2).

    Coefficients: [-3, 12, 17, 12, -3] / 35
    """
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size < 5:
        return arr.copy()
    kernel = np.array([-3.0, 12.0, 17.0, 12.0, -3.0], dtype=float) / 35.0
    padded = np.pad(arr, (2, 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def _parse_adsorbed_range_txt(range_txt_path: Path) -> tuple[float, float]:
    text = range_txt_path.read_text(encoding="utf-8")
    m_start = re.search(r"adsorbed_layer_start_A=([0-9.eE+-]+)", text)
    m_end = re.search(r"adsorbed_layer_end_A=([0-9.eE+-]+)", text)
    if m_start is None or m_end is None:
        raise ValueError(f"Cannot parse adsorbed layer range from {range_txt_path}")
    d_start_A = float(m_start.group(1))
    d_end_A = float(m_end.group(1))
    if d_end_A <= d_start_A:
        raise ValueError("adsorbed layer range is invalid: end must be greater than start")
    return d_start_A, d_end_A


def plot_water_three_panel_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    output_png_name: str = DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    ndeg: float = DEFAULT_THETA_BIN_DEG,
) -> Path:
    """
    Create an integrated three-panel figure:
    1) water mass density vs distance
    2) orientation-weighted density vs distance (share x with panel 1)
    3) adsorbed-layer theta distribution PDF (0-180 degree)
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)

    if output_dir is None:
        output_dir_path = DEFAULT_OUTPUT_DIR
    else:
        output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    density_csv_path = water_mass_density_z_distribution_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=output_dir_path,
        output_csv_name=DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        start_interface=start_interface,
        dz_A=dz_A,
    )
    orientation_csv_path = water_orientation_weighted_density_z_distribution_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=output_dir_path,
        output_csv_name=DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        start_interface=start_interface,
        dz_A=dz_A,
    )
    _, range_txt_path = ad_water_orientation_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=output_dir_path,
        start_interface=start_interface,
        dz_A=dz_A,
    )

    d_start_A, d_end_A = _parse_adsorbed_range_txt(range_txt_path)
    theta_centers, theta_pdf, _ = compute_adsorbed_water_theta_distribution(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        adsorbed_range_A=(d_start_A, d_end_A),
        output_dir=output_dir_path,
        start_interface=start_interface,
        dz_A=dz_A,
        ndeg=ndeg,
    )

    density_data = np.loadtxt(density_csv_path, delimiter=",", skiprows=1)
    if density_data.ndim == 1:
        density_data = density_data.reshape(1, -1)
    orientation_data = np.loadtxt(orientation_csv_path, delimiter=",", skiprows=1)
    if orientation_data.ndim == 1:
        orientation_data = orientation_data.reshape(1, -1)

    distance_A_density = density_data[:, 1]
    rho_g_cm3 = density_data[:, 2]
    distance_A_orientation = orientation_data[:, 1]
    orient_1_A3 = orientation_data[:, 2]
    rho_smooth = _savgol_smooth_window5(rho_g_cm3)
    orient_smooth = _savgol_smooth_window5(orient_1_A3)
    theta_pdf_smooth = _savgol_smooth_window5(theta_pdf)

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("matplotlib is required to generate plots.") from exc

    with plt.rc_context(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
        }
    ):
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.6, 10.8), sharex=False)

        # Panel 1: density profile.
        ax1.plot(distance_A_density, rho_smooth, color="tab:blue", lw=1.5)
        ax1.set_ylabel(r"$\rho_\mathrm{{H_2O}}\ \mathrm{(g\cdot cm^3)}$", fontsize=18)
        ax1.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax1.set_title("Water density distribution", fontsize=18)
        ax1.tick_params(axis="both", which="major", labelsize=16)

        # Panel 2: orientation-weighted profile, sharing x-range with panel 1.
        ax2.plot(distance_A_orientation, orient_smooth, color="tab:orange", lw=1.5)
        ax2.set_ylabel(r"$\rho_\mathrm{{H_2O}}\cdot cos \varphi\ \mathrm{(a.u.)}$", fontsize=18)
        ax2.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax2.set_title("Water dipole", fontsize=18)
        ax2.set_yticks([0.0])
        ax2.tick_params(axis="both", which="major", labelsize=16)

        x_min = 0.0
        x_max = float(max(np.max(distance_A_density), np.max(distance_A_orientation)))
        ax1.set_xlim(x_min, x_max)
        ax2.set_xlim(x_min, x_max)
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
        out_png_path = output_dir_path / output_png_name
        fig.savefig(out_png_path, dpi=180)
        plt.close(fig)
    return out_png_path
