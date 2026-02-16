"""High-level water density analysis along interface-to-midpoint direction."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Literal

import numpy as np

from ...utils.LayerParser import SurfaceGeometryError
from ...utils.LayerParser import detect_interface_layers
from ...utils.WaterParser import detect_water_molecule_indices
from ...utils.WaterParser import get_water_oxygen_indices_array
from ...utils.config import DEFAULT_Z_BIN_WIDTH_A
from ...utils.config import WATER_MOLAR_MASS_G_PER_MOL

from ..config import DEFAULT_OUTPUT_DIR
from ..config import DEFAULT_START_INTERFACE
from ..config import DEFAULT_WATER_MASS_DENSITY_CSV_NAME

StartInterface = Literal["low_c", "high_c"]
AVOGADRO_NUMBER = 6.022_140_76e23
ANGSTROM3_TO_CM3 = 1.0e-24


def _parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def _circular_mean_fractional(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size == 0:
        raise ValueError("cannot compute circular mean for empty values")
    angles = 2.0 * np.pi * values
    z = np.mean(np.cos(angles)) + 1j * np.mean(np.sin(angles))
    if z == 0:
        return float(np.mod(np.mean(values), 1.0))
    mean_angle = float(np.angle(z))
    if mean_angle < 0.0:
        mean_angle += 2.0 * np.pi
    return float(mean_angle / (2.0 * np.pi))


def _detect_low_high_interface_fractions(
    atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
) -> tuple[float, float]:
    detection = detect_interface_layers(
        atoms,
        normal="c",
        metal_symbols=metal_symbols,
        n_interface_layers=1,
    )
    interface_layers = detection.interface_layers()
    if len(interface_layers) < 2:
        raise SurfaceGeometryError(
            f"Expected at least 2 interface layers, got {len(interface_layers)}."
        )

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    interface_cfractions = []
    for layer in interface_layers:
        frac_values = scaled[list(layer.atom_indices), 2]
        interface_cfractions.append(_circular_mean_fractional(frac_values))

    interface_cfractions.sort()

    # `detect_interface_layers` now returns the two directly water-facing layers.
    # Keep compatibility fallback for legacy/atypical cases with more layers.
    if len(interface_cfractions) == 2:
        low_c = float(interface_cfractions[0])
        high_c = float(interface_cfractions[1])
    elif len(interface_cfractions) >= 4:
        low_c = float(interface_cfractions[1])
        high_c = float(interface_cfractions[-2])
    else:
        raise SurfaceGeometryError(
            "Cannot determine water-facing interface pair. "
            f"Need at least 2 layers and prefer 4, got {len(interface_cfractions)}."
        )

    if np.isclose(low_c, high_c):
        raise SurfaceGeometryError("Detected interface c-fractions are identical; cannot build analysis window.")
    return low_c, high_c


def _single_frame_density_profile_interface_to_midpoint(
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
        raise SurfaceGeometryError("interface_separation_c must be positive.")
    half_path_c = 0.5 * interface_separation_c

    water_indices = detect_water_molecule_indices(atoms)
    oxygen_indices = get_water_oxygen_indices_array(water_indices).reshape(-1)
    if oxygen_indices.size == 0:
        raise ValueError("No water oxygen indices found in current frame.")

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    oxygen_c = scaled[oxygen_indices, 2]

    if start_interface == "low_c":
        start_c = low_c
        delta_c = np.mod(oxygen_c - start_c, 1.0)
    else:
        start_c = high_c
        delta_c = np.mod(start_c - oxygen_c, 1.0)

    # Keep oxygens between selected interface and the midpoint of the two interfaces.
    in_path = (delta_c >= 0.0) & (delta_c <= half_path_c)
    selected_delta_c = delta_c[in_path]

    c_length_A = float(np.linalg.norm(np.asarray(atoms.cell.array, dtype=float)[2]))
    if c_length_A <= 0.0:
        raise ValueError("cell c-axis norm must be positive")
    selected_delta_A = selected_delta_c * c_length_A
    path_length_A = float(half_path_c * c_length_A)
    if path_length_A <= 0.0:
        raise ValueError("path_length_A must be positive")

    nbins = max(1, int(np.ceil(path_length_A / dz_A)))
    if nbins == 1:
        edges_A = np.array([0.0, path_length_A], dtype=float)
    else:
        interior = np.arange(1, nbins, dtype=float) * dz_A
        edges_A = np.concatenate(([0.0], interior, [path_length_A]))
    centers_A = 0.5 * (edges_A[:-1] + edges_A[1:])

    counts = np.histogram(selected_delta_A, bins=edges_A)[0].astype(float)

    cell = np.asarray(atoms.cell.array, dtype=float)
    area_xy_A2 = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    if area_xy_A2 <= 0.0:
        raise ValueError("cell area_xy must be positive")
    bin_volumes_A3 = area_xy_A2 * np.diff(edges_A)

    mass_per_water_g = WATER_MOLAR_MASS_G_PER_MOL / AVOGADRO_NUMBER
    mass_per_bin_g = counts * mass_per_water_g
    rho_g_cm3 = mass_per_bin_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)
    return centers_A, rho_g_cm3, path_length_A


def water_mass_density_z_distribution_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    output_csv_name: str = DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
) -> Path:
    """
    Ensemble-average water mass density profile from a selected metal interface
    toward the midpoint between two interfaces.

    Output CSV columns:
    - path_fraction_center: normalized distance in [0, 1] along interface->midpoint
    - distance_A: path_fraction_center * mean_path_length_A
    - rho_ensemble_avg_g_cm3: ensemble-averaged water mass density
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
        profile = _single_frame_density_profile_interface_to_midpoint(
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

    interpolated_rho = []
    for centers_A, rho_g_cm3, path_length_A in per_frame_profiles:
        centers_u = centers_A / path_length_A
        frame_interp = np.interp(common_centers_u, centers_u, rho_g_cm3, left=0.0, right=0.0)
        interpolated_rho.append(frame_interp)

    # A口径：每帧等权（equal-weight over frames）系综平均。
    rho_ensemble = np.mean(np.asarray(interpolated_rho, dtype=float), axis=0)
    distance_A = common_centers_u * mean_path_length_A

    out_csv_path = output_dir_path / output_csv_name
    out_data = np.column_stack([common_centers_u, distance_A, rho_ensemble])
    np.savetxt(
        out_csv_path,
        out_data,
        delimiter=",",
        header="path_fraction_center,distance_A,rho_ensemble_avg_g_cm3",
        comments="",
    )
    return out_csv_path
