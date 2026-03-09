"""
Shared private utilities for the WaterAnalysis layer.

Provides:
- Interface fraction detection
- Trajectory frame iterator
- Combined single-frame density + orientation profile
- Combined ensemble-average computation (single trajectory pass)

All names here are private (underscore-prefixed). External callers should use
the public API in WaterAnalysis/__init__.py.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Generator, Iterable, Literal

import numpy as np

logger = logging.getLogger(__name__)

try:
    from ase import Atoms
    from ase.io import iread
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]
    iread = None  # type: ignore[assignment]

from ...utils import (
    SurfaceGeometryError,
    detect_interface_layers,
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
)
from ...utils.StructureParser.WaterParser import (
    AVOGADRO_NUMBER,
    ANGSTROM3_TO_CM3,
    WaterTopologyError,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)
from ...utils.RestartParser.CellParser import parse_abc_from_md_inp as _parse_abc_from_md_inp
from ...utils.config import (
    DEFAULT_Z_BIN_WIDTH_A,
    INTERFACE_NORMAL_ALIGNED,
    INTERFACE_NORMAL_OPPOSED,
    WATER_MOLAR_MASS_G_PER_MOL,
)

StartInterface = Literal["normal_aligned", "normal_opposed"]


# ---------------------------------------------------------------------------
# Interface detection
# ---------------------------------------------------------------------------

def _detect_interface_fractions(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
) -> tuple[float, float]:
    """
    Return (aligned_frac, opposed_frac) of the two water-facing metal layers.

    Uses LayerParser to find interface layers and reads their ``center_frac``
    directly (already the circular-mean fractional coordinate).

    - ``aligned_frac`` = normal_aligned interface (outward normal = +axis)
    - ``opposed_frac`` = normal_opposed interface (outward normal = -axis)
    """
    detection = detect_interface_layers(
        atoms,
        normal="c",
        metal_symbols=metal_symbols,
    )

    aligned_frac = detection.interface_normal_aligned().center_frac
    opposed_frac = detection.interface_normal_opposed().center_frac

    if np.isclose(aligned_frac, opposed_frac):
        raise SurfaceGeometryError(
            "Detected interface c-fractions are identical; cannot build analysis window."
        )
    return aligned_frac, opposed_frac


# ---------------------------------------------------------------------------
# Trajectory iterator
# ---------------------------------------------------------------------------

def _iter_trajectory(
    xyz_path: Path,
    a_A: float,
    b_A: float,
    c_A: float,
    *,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
) -> Generator[Atoms, None, None]:
    """Yield ASE Atoms frames with orthogonal cell set and PBC enabled.

    Parameters
    ----------
    frame_start : int | None
        0-based index of first frame to yield (default: 0).
    frame_end : int | None
        0-based exclusive upper bound (default: no limit).
    frame_step : int | None
        Step between yielded frames (default: 1).
    """
    if iread is None:  # pragma: no cover
        raise RuntimeError("ASE is required: please install `ase` to read trajectory frames.")
    start = frame_start or 0
    stop = frame_end  # None = no upper bound
    step = frame_step or 1
    next_yield = start
    for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
        if stop is not None and idx >= stop:
            break
        if idx == next_yield:
            atoms.set_cell([a_A, b_A, c_A])
            atoms.set_pbc([True, True, True])
            yield atoms
            next_yield += step


# ---------------------------------------------------------------------------
# Combined per-frame computation
# ---------------------------------------------------------------------------

def _single_frame_density_and_orientation(
    atoms: Atoms,
    *,
    start_interface: StartInterface,
    dz_A: float,
    metal_symbols: Iterable[str] | None,
    compute_orientation: bool = True,
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
    if start_interface not in {INTERFACE_NORMAL_ALIGNED, INTERFACE_NORMAL_OPPOSED}:
        raise ValueError(f"start_interface must be 'normal_aligned' or 'normal_opposed', got {start_interface!r}")

    aligned_frac, opposed_frac = _detect_interface_fractions(atoms, metal_symbols=metal_symbols)
    gap_frac = (opposed_frac - aligned_frac) % 1.0
    if gap_frac <= 0.0:
        raise SurfaceGeometryError("Water gap fractional width must be positive.")
    half_path_c = 0.5 * gap_frac

    water_indices = detect_water_molecule_indices(atoms)
    oxygen_indices = get_water_oxygen_indices_array(water_indices).reshape(-1)
    if oxygen_indices.size == 0:
        raise ValueError("No water oxygen indices found in current frame.")
    o_to_h = _oxygen_to_hydrogen_map(water_indices)

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    oxygen_c = scaled[oxygen_indices, 2]

    if start_interface == INTERFACE_NORMAL_ALIGNED:
        delta_c = np.mod(oxygen_c - aligned_frac, 1.0)
    else:
        delta_c = np.mod(opposed_frac - oxygen_c, 1.0)

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

    # --- Orientation (vectorized, skipped when compute_orientation=False) ---
    if compute_orientation and selected_oxygen.size > 0:
        cos_theta = _compute_bisector_cos_theta_vec(atoms, selected_oxygen, o_to_h, c_unit)
        cos_sum = np.histogram(selected_delta_A, bins=edges_A, weights=cos_theta)[0].astype(float)
        orient_g_cm3 = cos_sum * mass_per_water_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)
    else:
        orient_g_cm3 = np.zeros_like(centers_A, dtype=float)

    return centers_A, rho_g_cm3, orient_g_cm3, path_length_A


# ---------------------------------------------------------------------------
# Combined ensemble average (single trajectory pass)
# ---------------------------------------------------------------------------

def _compute_density_orientation_ensemble(
    xyz_path: Path,
    md_inp_path: Path,
    *,
    start_interface: StartInterface = "normal_aligned",
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
    compute_orientation: bool = True,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
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
    logger.info("Computing density+orientation ensemble")

    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)

    iterator = _iter_trajectory(
        xyz_path, a_A, b_A, c_A,
        frame_start=frame_start, frame_end=frame_end, frame_step=frame_step,
    )
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(iterator, desc="Water density+orientation", unit="frame", ascii=" =")

    per_frame: list[tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []
    for atoms in iterator:
        per_frame.append(
            _single_frame_density_and_orientation(
                atoms,
                start_interface=start_interface,
                dz_A=dz_A,
                metal_symbols=metal_symbols,
                compute_orientation=compute_orientation,
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
