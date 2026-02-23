"""Integrated three-panel water analysis plotting utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import NullFormatter

from .config import DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
from .config import DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME
from .config import DEFAULT_OUTPUT_DIR
from .config import DEFAULT_START_INTERFACE
from .config import DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME
from .WaterAnalysis import compute_adsorbed_water_theta_distribution
from .WaterAnalysis import detect_adsorbed_layer_range_from_density_profile
from .WaterAnalysis._common import StartInterface, _compute_density_orientation_ensemble
from ..utils.config import DEFAULT_THETA_BIN_DEG, DEFAULT_Z_BIN_WIDTH_A


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
    2) orientation-weighted density vs distance (shared x range with panel 1)
    3) adsorbed-layer theta distribution PDF (0-180 degree)

    The trajectory is read exactly twice:
    - Read #1: compute density + orientation profiles (single combined pass).
    - Read #2: collect theta values for adsorbed-layer molecules.
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)
    output_dir_path = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # --- Trajectory read #1: density + orientation in one pass ---
    common_centers_u, mean_path_A, rho_ensemble, orient_ensemble = _compute_density_orientation_ensemble(
        xyz_path,
        md_inp_path,
        start_interface=start_interface,
        dz_A=dz_A,
    )
    distance_A = common_centers_u * mean_path_A

    # Save density CSV
    density_csv_path = output_dir_path / DEFAULT_WATER_MASS_DENSITY_CSV_NAME
    np.savetxt(
        density_csv_path,
        np.column_stack([common_centers_u, distance_A, rho_ensemble]),
        delimiter=",",
        header="path_fraction_center,distance_A,rho_ensemble_avg_g_cm3",
        comments="",
    )

    # Save orientation CSV
    orientation_csv_path = output_dir_path / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
    np.savetxt(
        orientation_csv_path,
        np.column_stack([common_centers_u, distance_A, orient_ensemble]),
        delimiter=",",
        header="path_fraction_center,distance_A,orientation_ensemble_avg_g_cm3",
        comments="",
    )

    # Detect adsorbed layer range in-memory (no trajectory read)
    d_start_A, d_end_A, d_peak_A = detect_adsorbed_layer_range_from_density_profile(
        distance_A, rho_ensemble
    )
    in_adsorbed = (distance_A >= d_start_A) & (distance_A <= d_end_A)

    # Save adsorbed profile CSV and range TXT (side effects expected by callers)
    adsorbed_profile_path = output_dir_path / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
    np.savetxt(
        adsorbed_profile_path,
        np.column_stack([distance_A, rho_ensemble, orient_ensemble, in_adsorbed.astype(int)]),
        delimiter=",",
        header="distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_g_cm3,is_adsorbed_layer_bin",
        comments="",
    )

    range_txt_path = output_dir_path / DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME
    range_txt_path.write_text(
        "\n".join([
            f"adsorbed_layer_start_A={d_start_A:.10f}",
            f"adsorbed_layer_end_A={d_end_A:.10f}",
            f"main_peak_distance_A={d_peak_A:.10f}",
            "near_zero_ratio=0.050000",
            "smoothing_window_bins=5",
        ]),
        encoding="utf-8",
    )

    # --- Trajectory read #2: adsorbed-layer theta distribution ---
    theta_centers, theta_pdf, _ = compute_adsorbed_water_theta_distribution(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        adsorbed_range_A=(d_start_A, d_end_A),
        output_dir=output_dir_path,
        start_interface=start_interface,
        dz_A=dz_A,
        ndeg=ndeg,
    )

    # --- Smooth and plot ---
    rho_smooth = _savgol_smooth_window5(rho_ensemble)
    orient_smooth = _savgol_smooth_window5(orient_ensemble)
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
        ax1.plot(distance_A, rho_smooth, color="tab:blue", lw=1.5)
        ax1.set_ylabel(r"$\rho_\mathrm{{H_2O}}\ \mathrm{(g\cdot cm^3)}$", fontsize=18)
        ax1.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax1.set_title("Water density distribution", fontsize=18)
        ax1.tick_params(axis="both", which="major", labelsize=16)

        # Panel 2: orientation-weighted profile.
        ax2.plot(distance_A, orient_smooth, color="tab:orange", lw=1.5)
        ax2.set_ylabel(r"$\rho_\mathrm{{H_2O}}\cdot\cos\varphi\ \mathrm{(g\cdot cm^{-3})}$", fontsize=18)
        ax2.set_xlabel(r"$\mathrm{distance(\AA)}$", fontsize=18)
        ax2.set_title("Water dipole", fontsize=18)
        ax2.set_yticks([0.0])
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
        out_png_path = output_dir_path / output_png_name
        fig.savefig(out_png_path, dpi=180)
        plt.close(fig)

    return out_png_path
