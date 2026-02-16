"""High-level water orientation-weighted density analysis along interface-to-midpoint direction."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from ...utils.WaterParser import WaterTopologyError
from ...utils.WaterParser import detect_water_molecule_indices
from ...utils.WaterParser import get_water_oxygen_indices_array
from ...utils.config import DEFAULT_Z_BIN_WIDTH_A
from ..config import DEFAULT_OUTPUT_DIR
from ..config import DEFAULT_START_INTERFACE
from ..config import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from .WaterDensity import StartInterface
from .WaterDensity import _detect_low_high_interface_fractions
from .WaterDensity import _parse_abc_from_md_inp


def _oxygen_to_hydrogen_map(water_molecule_indices: np.ndarray) -> dict[int, tuple[int, int]]:
    water_molecule_indices = np.asarray(water_molecule_indices, dtype=int)
    if water_molecule_indices.ndim != 2 or water_molecule_indices.shape[1] != 3:
        raise ValueError("water_molecule_indices must have shape (n_water, 3)")

    mapping: dict[int, tuple[int, int]] = {}
    for row in water_molecule_indices:
        o_idx = int(row[0])
        h1_idx = int(row[1])
        h2_idx = int(row[2])
        if o_idx in mapping:
            raise ValueError(f"duplicate oxygen index found in water_molecule_indices: {o_idx}")
        mapping[o_idx] = (h1_idx, h2_idx)
    return mapping


def _single_frame_orientation_profile_interface_to_midpoint(
    atoms,
    *,
    start_interface: StartInterface,
    dz_A: float,
    metal_symbols: Iterable[str] | None,
) -> tuple[np.ndarray, np.ndarray, float]:
    if dz_A <= 0.0:
        raise ValueError("dz_A must be > 0")
    if start_interface not in {"low_c", "high_c"}:
        raise ValueError(f"start_interface must be 'low_c' or 'high_c', got {start_interface!r}")

    low_c, high_c = _detect_low_high_interface_fractions(atoms, metal_symbols=metal_symbols)
    interface_separation_c = high_c - low_c
    if interface_separation_c <= 0.0:
        raise ValueError("interface_separation_c must be positive.")
    half_path_c = 0.5 * interface_separation_c

    water_indices = detect_water_molecule_indices(atoms)
    oxygen_indices = get_water_oxygen_indices_array(water_indices).reshape(-1)
    o_to_h = _oxygen_to_hydrogen_map(water_indices)

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    oxygen_c = scaled[oxygen_indices, 2]
    if start_interface == "low_c":
        delta_c = np.mod(oxygen_c - low_c, 1.0)
    else:
        delta_c = np.mod(high_c - oxygen_c, 1.0)

    in_path = (delta_c >= 0.0) & (delta_c <= half_path_c)
    selected_oxygen = oxygen_indices[in_path]
    selected_delta_c = delta_c[in_path]

    cell = np.asarray(atoms.cell.array, dtype=float)
    c_vec = cell[2]
    c_norm = float(np.linalg.norm(c_vec))
    if c_norm <= 0.0:
        raise ValueError("cell c-axis norm must be positive")
    c_unit = c_vec / c_norm

    path_length_A = float(half_path_c * c_norm)
    if path_length_A <= 0.0:
        raise ValueError("path_length_A must be positive")
    selected_delta_A = selected_delta_c * c_norm

    nbins = max(1, int(np.ceil(path_length_A / dz_A)))
    if nbins == 1:
        edges_A = np.array([0.0, path_length_A], dtype=float)
    else:
        interior = np.arange(1, nbins, dtype=float) * dz_A
        edges_A = np.concatenate(([0.0], interior, [path_length_A]))
    centers_A = 0.5 * (edges_A[:-1] + edges_A[1:])

    area_xy_A2 = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    if area_xy_A2 <= 0.0:
        raise ValueError("cell area_xy must be positive")
    bin_volumes_A3 = area_xy_A2 * np.diff(edges_A)

    if selected_oxygen.size == 0:
        return centers_A, np.zeros_like(centers_A, dtype=float), path_length_A

    cos_theta_values = np.empty(selected_oxygen.size, dtype=float)
    for i, o_idx in enumerate(selected_oxygen):
        o_int = int(o_idx)
        if o_int not in o_to_h:
            raise WaterTopologyError(f"Oxygen index {o_int} is not labeled as water molecule")
        h1_idx, h2_idx = o_to_h[o_int]
        vecs = np.asarray(atoms.get_distances(o_int, [int(h1_idx), int(h2_idx)], vector=True, mic=True), dtype=float)
        v1 = vecs[0]
        v2 = vecs[1]

        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 == 0.0 or n2 == 0.0:
            raise WaterTopologyError(f"Zero O-H vector norm for oxygen index {o_int}")

        bisector = v1 / n1 + v2 / n2
        nb = float(np.linalg.norm(bisector))
        if nb == 0.0:
            raise WaterTopologyError(f"Degenerate H-O-H bisector for oxygen index {o_int}")

        cos_theta_values[i] = float(np.dot(bisector, c_unit) / nb)

    cos_sum_per_bin = np.histogram(selected_delta_A, bins=edges_A, weights=cos_theta_values)[0].astype(float)
    orient_density_1_A3 = cos_sum_per_bin / bin_volumes_A3
    return centers_A, orient_density_1_A3, path_length_A


def water_orientation_weighted_density_z_distribution_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    output_csv_name: str = DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
) -> Path:
    """
    Ensemble-average orientation-weighted density profile from a selected metal interface
    toward the midpoint between two interfaces.

    Output CSV columns:
    - path_fraction_center: normalized distance in [0, 1] along interface->midpoint
    - distance_A: path_fraction_center * mean_path_length_A
    - orientation_ensemble_avg_1_A3: ensemble-averaged orientation-weighted density
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)

    if output_dir is None:
        output_dir_path = DEFAULT_OUTPUT_DIR
    else:
        output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    try:
        from ase.io import iread
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("ASE is required: please install `ase` to read trajectory frames.") from exc

    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)

    per_frame_profiles: list[tuple[np.ndarray, np.ndarray, float]] = []
    for atoms in iread(str(xyz_path), index=":"):
        atoms.set_cell([a_A, b_A, c_A])
        atoms.set_pbc([True, True, True])
        profile = _single_frame_orientation_profile_interface_to_midpoint(
            atoms,
            start_interface=start_interface,
            dz_A=dz_A,
            metal_symbols=metal_symbols,
        )
        per_frame_profiles.append(profile)

    if not per_frame_profiles:
        raise RuntimeError("No frames were read from trajectory.")

    mean_path_length_A = float(np.mean([p[2] for p in per_frame_profiles]))
    n_bins = max(1, int(np.ceil(mean_path_length_A / dz_A)))
    common_edges_u = np.linspace(0.0, 1.0, n_bins + 1, dtype=float)
    common_centers_u = 0.5 * (common_edges_u[:-1] + common_edges_u[1:])

    interpolated_orient = []
    for centers_A, orient_1_A3, path_length_A in per_frame_profiles:
        centers_u = centers_A / path_length_A
        frame_interp = np.interp(common_centers_u, centers_u, orient_1_A3, left=0.0, right=0.0)
        interpolated_orient.append(frame_interp)

    # A口径：每帧等权（equal-weight over frames）系综平均。
    orient_ensemble = np.mean(np.asarray(interpolated_orient, dtype=float), axis=0)
    distance_A = common_centers_u * mean_path_length_A

    out_csv_path = output_dir_path / output_csv_name
    out_data = np.column_stack([common_centers_u, distance_A, orient_ensemble])
    np.savetxt(
        out_csv_path,
        out_data,
        delimiter=",",
        header="path_fraction_center,distance_A,orientation_ensemble_avg_1_A3",
        comments="",
    )
    return out_csv_path
