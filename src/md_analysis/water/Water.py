"""Integrated three-panel water analysis plotting utilities."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

from ..utils.io._io_helpers import _write_csv_from_arrays

from ._plot import plot_three_panel as _plot_three_panel
from .config import DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
from .config import DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME
from .config import DEFAULT_OUTPUT_DIR
from .config import DEFAULT_START_INTERFACE
from .config import DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME
from .WaterAnalysis import compute_adsorbed_water_theta_distribution
from .WaterAnalysis import detect_adsorbed_layer_range_from_density_profile
from .WaterAnalysis import StartInterface, _compute_density_orientation_ensemble
from ..utils.constants import DEFAULT_LAYER_TOL_A, DEFAULT_THETA_BIN_DEG, DEFAULT_Z_BIN_WIDTH_A


def plot_water_three_panel_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: str | Path | None = None,
    output_png_name: str = DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    ndeg: float = DEFAULT_THETA_BIN_DEG,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
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
    logger.info("Starting three-panel water analysis")

    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path) if md_inp_path is not None else None
    output_dir_path = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # --- Trajectory read #1: density + orientation in one pass ---
    common_centers_u, mean_path_A, rho_ensemble, orient_ensemble = _compute_density_orientation_ensemble(
        xyz_path,
        md_inp_path,
        cell_abc=cell_abc,
        start_interface=start_interface,
        dz_A=dz_A,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )
    distance_A = common_centers_u * mean_path_A

    # Save density CSV
    density_csv_path = output_dir_path / DEFAULT_WATER_MASS_DENSITY_CSV_NAME
    _write_csv_from_arrays(density_csv_path, {
        "path_fraction_center": common_centers_u,
        "distance_A": distance_A,
        "rho_ensemble_avg_g_cm3": rho_ensemble,
    })

    # Save orientation CSV
    orientation_csv_path = output_dir_path / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
    _write_csv_from_arrays(orientation_csv_path, {
        "path_fraction_center": common_centers_u,
        "distance_A": distance_A,
        "orientation_ensemble_avg_g_cm3": orient_ensemble,
    })

    # Detect adsorbed layer range in-memory (no trajectory read)
    d_start_A, d_end_A, d_peak_A = detect_adsorbed_layer_range_from_density_profile(
        distance_A, rho_ensemble
    )
    logger.info("Adsorbed layer range: %.3f - %.3f A", d_start_A, d_end_A)
    in_adsorbed = (distance_A >= d_start_A) & (distance_A <= d_end_A)

    # Save adsorbed profile CSV and range TXT (side effects expected by callers)
    adsorbed_profile_path = output_dir_path / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
    _write_csv_from_arrays(adsorbed_profile_path, {
        "distance_A": distance_A,
        "rho_ensemble_avg_g_cm3": rho_ensemble,
        "orientation_ensemble_avg_g_cm3": orient_ensemble,
        "is_adsorbed_layer_bin": in_adsorbed.astype(float),
    })

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
        cell_abc=cell_abc,
        adsorbed_range_A=(d_start_A, d_end_A),
        output_dir=output_dir_path,
        start_interface=start_interface,
        dz_A=dz_A,
        ndeg=ndeg,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    # --- Plotting ---
    out_png_path = output_dir_path / output_png_name
    _plot_three_panel(
        out_png_path,
        distance_A=distance_A,
        rho_ensemble=rho_ensemble,
        orient_ensemble=orient_ensemble,
        theta_centers=theta_centers,
        theta_pdf=theta_pdf,
    )

    logger.info("Three-panel PNG: %s", out_png_path)
    return out_png_path
