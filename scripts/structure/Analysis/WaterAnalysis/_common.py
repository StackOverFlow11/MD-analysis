"""
Shared private utilities for the WaterAnalysis layer.

Provides:
- Cell parameter parsing from CP2K md.inp
- Interface fraction detection
- Trajectory frame iterator
- Combined single-frame density + orientation profile
- Combined ensemble-average computation (single trajectory pass)

All names here are private (underscore-prefixed). External callers should use
the public API in WaterAnalysis/__init__.py.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Generator, Iterable, Literal

import numpy as np

try:
    from ase import Atoms
    from ase.io import iread
except Exception:  # pragma: no cover
    Atoms = object  # type: ignore
    iread = None  # type: ignore

from ...utils.LayerParser import SurfaceGeometryError, _circular_mean_fractional, detect_interface_layers
from ...utils.WaterParser import (
    AVOGADRO_NUMBER,
    ANGSTROM3_TO_CM3,
    WaterTopologyError,
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)
from ...utils.config import DEFAULT_Z_BIN_WIDTH_A, WATER_MOLAR_MASS_G_PER_MOL

StartInterface = Literal["low_c", "high_c"]


# ---------------------------------------------------------------------------
# Cell parameter parsing
# ---------------------------------------------------------------------------

def _parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    """Parse orthogonal cell lengths (Å) from a CP2K md.inp file."""
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(
        r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$",
        text,
        re.MULTILINE,
    )
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


# ---------------------------------------------------------------------------
# Interface detection
# ---------------------------------------------------------------------------

def _detect_low_high_interface_fractions(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
) -> tuple[float, float]:
    """
    Return (low_c_fraction, high_c_fraction) of the two water-facing metal layers.

    Uses LayerParser to find interface layers, then computes their circular-mean
    fractional c-coordinate.
    """
    detection = detect_interface_layers(
        atoms,
        normal="c",
        metal_symbols=metal_symbols,
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
        raise SurfaceGeometryError(
            "Detected interface c-fractions are identical; cannot build analysis window."
        )
    return low_c, high_c


# ---------------------------------------------------------------------------
# Trajectory iterator
# ---------------------------------------------------------------------------

def _iter_trajectory(
    xyz_path: Path,
    a_A: float,
    b_A: float,
    c_A: float,
) -> Generator[Atoms, None, None]:
    """Yield ASE Atoms frames with orthogonal cell set and PBC enabled."""
    if iread is None:  # pragma: no cover
        raise RuntimeError("ASE is required: please install `ase` to read trajectory frames.")
    for atoms in iread(str(xyz_path), index=":"):
        atoms.set_cell([a_A, b_A, c_A])
        atoms.set_pbc([True, True, True])
        yield atoms


# ---------------------------------------------------------------------------
# Combined per-frame computation
# ---------------------------------------------------------------------------

def _single_frame_density_and_orientation(
    atoms: Atoms,
    *,
    start_interface: StartInterface,
    dz_A: float,
    metal_symbols: Iterable[str] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Compute both density and orientation profiles for a single frame.

    Returns
    -------
    centers_A : np.ndarray, shape (nbins,)
        Bin center positions from interface (Å).
    rho_g_cm3 : np.ndarray, shape (nbins,)
        Water mass density per bin (g/cm³).
    orient_g_cm3 : np.ndarray, shape (nbins,)
        Orientation-weighted mass density per bin (g/cm³).
    path_length_A : float
        Physical half-path length from interface to midpoint (Å).
    """
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
    c_norm = float(np.linalg.norm(cell[2]))
    if c_norm <= 0.0:
        raise ValueError("cell c-axis norm must be positive")
    c_unit = cell[2] / c_norm

    path_length_A = float(half_path_c * c_norm)
    if path_length_A <= 0.0:
        raise ValueError("path_length_A must be positive")
    selected_delta_A = selected_delta_c * c_norm

    # Build bin geometry
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

    # --- Density ---
    counts = np.histogram(selected_delta_A, bins=edges_A)[0].astype(float)
    mass_per_water_g = WATER_MOLAR_MASS_G_PER_MOL / AVOGADRO_NUMBER
    rho_g_cm3 = counts * mass_per_water_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)

    # --- Orientation (vectorized) ---
    if selected_oxygen.size == 0:
        orient_g_cm3 = np.zeros_like(centers_A, dtype=float)
    else:
        cos_theta = _compute_bisector_cos_theta_vec(atoms, selected_oxygen, o_to_h, c_unit)
        cos_sum = np.histogram(selected_delta_A, bins=edges_A, weights=cos_theta)[0].astype(float)
        orient_g_cm3 = cos_sum * mass_per_water_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)

    return centers_A, rho_g_cm3, orient_g_cm3, path_length_A


# ---------------------------------------------------------------------------
# Combined ensemble average (single trajectory pass)
# ---------------------------------------------------------------------------

def _compute_density_orientation_ensemble(
    xyz_path: Path,
    md_inp_path: Path,
    *,
    start_interface: StartInterface = "low_c",
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """
    Read the trajectory **once** and compute ensemble-averaged density and orientation.

    Returns
    -------
    common_centers_u : np.ndarray, shape (nbins,)
        Normalized bin centers on [0, 1].
    mean_path_length_A : float
        Mean physical half-path length (Å), averaged over frames.
    rho_ensemble : np.ndarray, shape (nbins,)
        Equal-weight ensemble-averaged water mass density (g/cm³).
    orient_ensemble : np.ndarray, shape (nbins,)
        Equal-weight ensemble-averaged orientation-weighted mass density (g/cm³).
    """
    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)

    per_frame: list[tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []
    for atoms in _iter_trajectory(xyz_path, a_A, b_A, c_A):
        per_frame.append(
            _single_frame_density_and_orientation(
                atoms,
                start_interface=start_interface,
                dz_A=dz_A,
                metal_symbols=metal_symbols,
            )
        )

    if not per_frame:
        raise RuntimeError("No frames were read from trajectory.")

    mean_path_length_A = float(np.mean([d[3] for d in per_frame]))
    n_bins = max(1, int(np.ceil(mean_path_length_A / dz_A)))
    common_edges_u = np.linspace(0.0, 1.0, n_bins + 1, dtype=float)
    common_centers_u = 0.5 * (common_edges_u[:-1] + common_edges_u[1:])

    interp_rho: list[np.ndarray] = []
    interp_orient: list[np.ndarray] = []
    for centers_A, rho_g_cm3, orient_g_cm3, path_length_A in per_frame:
        centers_u = centers_A / path_length_A
        interp_rho.append(np.interp(common_centers_u, centers_u, rho_g_cm3, left=0.0, right=0.0))
        interp_orient.append(np.interp(common_centers_u, centers_u, orient_g_cm3, left=0.0, right=0.0))

    rho_ensemble = np.mean(np.asarray(interp_rho, dtype=float), axis=0)
    orient_ensemble = np.mean(np.asarray(interp_orient, dtype=float), axis=0)

    return common_centers_u, mean_path_length_A, rho_ensemble, orient_ensemble
