"""Adsorbed-water orientation analysis driven by water-density profile."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from ...utils.WaterParser import (
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
    _theta_bin_count_from_ndeg,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)
from ...utils.config import DEFAULT_THETA_BIN_DEG, DEFAULT_Z_BIN_WIDTH_A
from ..config import (
    DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
    DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
    DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_START_INTERFACE,
)
from ._common import (
    StartInterface,
    _compute_density_orientation_ensemble,
    _detect_low_high_interface_fractions,
    _iter_trajectory,
    _parse_abc_from_md_inp,
)

__all__ = [
    "detect_adsorbed_layer_range_from_density_profile",
    "ad_water_orientation_analysis",
    "compute_adsorbed_water_theta_distribution",
]


# ---------------------------------------------------------------------------
# Smoothing helper (local to this module)
# ---------------------------------------------------------------------------

def _smooth_1d(values: np.ndarray, window_bins: int) -> np.ndarray:
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size == 0:
        return values.copy()
    if window_bins <= 1:
        return values.copy()
    if window_bins % 2 == 0:
        window_bins += 1
    if window_bins > values.size:
        window_bins = values.size if values.size % 2 == 1 else max(1, values.size - 1)
    if window_bins <= 1:
        return values.copy()
    kernel = np.ones(window_bins, dtype=float) / float(window_bins)
    return np.convolve(values, kernel, mode="same")


# ---------------------------------------------------------------------------
# Adsorbed layer range detection
# ---------------------------------------------------------------------------

def detect_adsorbed_layer_range_from_density_profile(
    distance_A: np.ndarray,
    rho_g_cm3: np.ndarray,
    *,
    near_zero_ratio: float = 0.05,
    smoothing_window_bins: int = 5,
) -> tuple[float, float, float]:
    """
    Detect adsorbed-water layer range from density profile.

    Rule set:
    1) main peak position = distance of the maximum-density bin (user-defined)
    2) lower bound = last near-zero bin before main peak
    3) upper bound = first local minimum after main peak (on smoothed profile)
    """
    d = np.asarray(distance_A, dtype=float).reshape(-1)
    rho = np.asarray(rho_g_cm3, dtype=float).reshape(-1)
    if d.size != rho.size:
        raise ValueError("distance_A and rho_g_cm3 must have same length")
    if d.size < 3:
        raise ValueError("profile must contain at least 3 bins")
    if near_zero_ratio < 0.0:
        raise ValueError("near_zero_ratio must be >= 0")

    peak_idx = int(np.argmax(rho))
    peak_distance_A = float(d[peak_idx])
    peak_value = float(rho[peak_idx])

    near_zero_threshold = near_zero_ratio * peak_value
    left_candidates = np.where((np.arange(d.size) <= peak_idx) & (rho <= near_zero_threshold))[0]
    start_idx = int(left_candidates[-1]) if left_candidates.size > 0 else 0

    rho_smooth = _smooth_1d(rho, window_bins=smoothing_window_bins)
    end_idx: int | None = None
    for i in range(max(peak_idx + 1, 1), d.size - 1):
        if rho_smooth[i] <= rho_smooth[i - 1] and rho_smooth[i] <= rho_smooth[i + 1]:
            end_idx = i
            break
    if end_idx is None:
        if peak_idx < d.size - 1:
            end_idx = int(np.argmin(rho_smooth[peak_idx + 1 :])) + peak_idx + 1
        else:
            end_idx = peak_idx

    if end_idx <= peak_idx and peak_idx < d.size - 1:
        end_idx = peak_idx + 1
    if start_idx >= end_idx:
        start_idx = max(0, min(start_idx, end_idx - 1))

    start_distance_A = float(d[start_idx])
    end_distance_A = float(d[end_idx])
    return start_distance_A, end_distance_A, peak_distance_A


# ---------------------------------------------------------------------------
# Adsorbed water theta distribution
# ---------------------------------------------------------------------------

def compute_adsorbed_water_theta_distribution(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    adsorbed_range_A: tuple[float, float] | None = None,
    output_dir: str | Path | None = None,
    output_csv_name: str = DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    ndeg: float = DEFAULT_THETA_BIN_DEG,
    near_zero_ratio: float = 0.05,
    smoothing_window_bins: int = 5,
) -> tuple[np.ndarray, np.ndarray, Path]:
    """
    Compute theta distribution (0-180 degree) for waters in adsorbed layer.

    - If adsorbed_range_A is None, auto-detect range from density profile.
    - Returns (theta_centers_deg, theta_pdf_degree_inv, csv_path).
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)
    output_dir_path = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    if adsorbed_range_A is None:
        # One trajectory pass to get density, then derive range.
        common_centers_u, mean_path_A, rho_ensemble, _ = _compute_density_orientation_ensemble(
            xyz_path,
            md_inp_path,
            start_interface=start_interface,
            dz_A=dz_A,
        )
        distance_A_tmp = common_centers_u * mean_path_A
        d_start_A, d_end_A, _ = detect_adsorbed_layer_range_from_density_profile(
            distance_A_tmp,
            rho_ensemble,
            near_zero_ratio=near_zero_ratio,
            smoothing_window_bins=smoothing_window_bins,
        )
    else:
        d_start_A = float(adsorbed_range_A[0])
        d_end_A = float(adsorbed_range_A[1])
    if d_end_A <= d_start_A:
        raise ValueError("adsorbed_range_A must satisfy end > start")

    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)
    n_bins = _theta_bin_count_from_ndeg(float(ndeg))
    theta_values_deg: list[np.ndarray] = []

    for atoms in _iter_trajectory(xyz_path, a_A, b_A, c_A):
        low_c, high_c = _detect_low_high_interface_fractions(atoms)
        water_idx = detect_water_molecule_indices(atoms)
        oxygen_indices = get_water_oxygen_indices_array(water_idx).reshape(-1)
        o_to_h = _oxygen_to_hydrogen_map(water_idx)

        cell = np.asarray(atoms.cell.array, dtype=float)
        c_norm = float(np.linalg.norm(cell[2]))
        if c_norm <= 0.0:
            raise ValueError("cell c-axis norm must be positive")
        c_unit = cell[2] / c_norm

        scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
        oxygen_c = scaled[oxygen_indices, 2]
        if start_interface == "low_c":
            delta_c = np.mod(oxygen_c - low_c, 1.0)
        else:
            delta_c = np.mod(high_c - oxygen_c, 1.0)
        delta_A = delta_c * c_norm

        in_adsorbed = (delta_A >= d_start_A) & (delta_A <= d_end_A)
        selected_oxygen = oxygen_indices[in_adsorbed]
        if selected_oxygen.size == 0:
            continue

        # Vectorized bisector computation
        cos_theta = _compute_bisector_cos_theta_vec(atoms, selected_oxygen, o_to_h, c_unit)
        cos_theta_clipped = np.clip(cos_theta, -1.0, 1.0)
        theta_values_deg.append(np.degrees(np.arccos(cos_theta_clipped)))

    theta_edges = np.linspace(0.0, 180.0, n_bins + 1, dtype=float)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])

    if theta_values_deg:
        all_theta = np.concatenate(theta_values_deg)
        theta_pdf, _ = np.histogram(all_theta, bins=theta_edges, density=True)
        theta_pdf = theta_pdf.astype(float, copy=False)
    else:
        theta_pdf = np.zeros(n_bins, dtype=float)

    out_csv_path = output_dir_path / output_csv_name
    np.savetxt(
        out_csv_path,
        np.column_stack([theta_centers, theta_pdf]),
        delimiter=",",
        header="theta_degree,pdf_degree_inv",
        comments="",
    )
    return theta_centers, theta_pdf, out_csv_path


# ---------------------------------------------------------------------------
# Adsorbed water orientation profile
# ---------------------------------------------------------------------------

def ad_water_orientation_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    output_profile_csv_name: str = DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
    output_range_txt_name: str = DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    near_zero_ratio: float = 0.05,
    smoothing_window_bins: int = 5,
) -> tuple[Path, Path]:
    """
    Build adsorbed-water orientation profile based on detected adsorbed layer range.

    Single trajectory pass: density and orientation are computed simultaneously.

    Returns:
    - profile CSV path (distance, density, orientation, adsorbed-mask)
    - range TXT path (start/end/peak)
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)
    output_dir_path = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Single trajectory pass — replaces two separate analysis calls.
    common_centers_u, mean_path_A, rho_ensemble, orient_ensemble = _compute_density_orientation_ensemble(
        xyz_path,
        md_inp_path,
        start_interface=start_interface,
        dz_A=dz_A,
    )
    distance_A = common_centers_u * mean_path_A

    d_start, d_end, d_peak = detect_adsorbed_layer_range_from_density_profile(
        distance_A,
        rho_ensemble,
        near_zero_ratio=near_zero_ratio,
        smoothing_window_bins=smoothing_window_bins,
    )
    in_adsorbed = (distance_A >= d_start) & (distance_A <= d_end)

    profile_csv_path = output_dir_path / output_profile_csv_name
    np.savetxt(
        profile_csv_path,
        np.column_stack([distance_A, rho_ensemble, orient_ensemble, in_adsorbed.astype(int)]),
        delimiter=",",
        header="distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_g_cm3,is_adsorbed_layer_bin",
        comments="",
    )

    range_txt_path = output_dir_path / output_range_txt_name
    range_txt_path.write_text(
        "\n".join([
            f"adsorbed_layer_start_A={d_start:.10f}",
            f"adsorbed_layer_end_A={d_end:.10f}",
            f"main_peak_distance_A={d_peak:.10f}",
            f"near_zero_ratio={near_zero_ratio:.6f}",
            f"smoothing_window_bins={int(smoothing_window_bins)}",
        ]),
        encoding="utf-8",
    )
    return profile_csv_path, range_txt_path
