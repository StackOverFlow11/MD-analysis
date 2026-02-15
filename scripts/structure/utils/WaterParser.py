"""
Water parsing and z-axis distribution utilities for a single ASE frame.

This module provides:
1) water-molecule labeling as an (n_water, 3) integer array [O, H1, H2]
2) water mass-density z distribution (g/cm^3) from oxygen indices
3) plane-density-weighted orientation z distribution from oxygen indices
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

try:
    from ase import Atoms
except Exception:  # pragma: no cover
    Atoms = object  # type: ignore

try:
    from .config import DEFAULT_THETA_BIN_DEG
    from .config import DEFAULT_WATER_OH_CUTOFF_A
    from .config import DEFAULT_Z_BIN_WIDTH_A
    from .config import WATER_MOLAR_MASS_G_PER_MOL
except Exception:  # pragma: no cover
    from config import DEFAULT_THETA_BIN_DEG  # type: ignore
    from config import DEFAULT_WATER_OH_CUTOFF_A  # type: ignore
    from config import DEFAULT_Z_BIN_WIDTH_A  # type: ignore
    from config import WATER_MOLAR_MASS_G_PER_MOL  # type: ignore


AVOGADRO_NUMBER = 6.022_140_76e23
ANGSTROM3_TO_CM3 = 1.0e-24


class WaterTopologyError(RuntimeError):
    """Raised when water topology cannot be inferred robustly."""


def _as_flat_index_array(indices: Sequence[int] | np.ndarray, *, name: str) -> np.ndarray:
    arr = np.asarray(indices)
    if arr.size == 0:
        raise ValueError(f"{name} must be non-empty")
    arr = arr.reshape(-1)
    if np.issubdtype(arr.dtype, np.integer):
        return arr.astype(int, copy=False)

    rounded = np.rint(arr).astype(int)
    if not np.allclose(arr, rounded):
        raise ValueError(f"{name} must contain integer-like values")
    return rounded


def _validate_atom_indices(indices: np.ndarray, natoms: int, *, name: str) -> None:
    if np.any(indices < 0) or np.any(indices >= natoms):
        raise ValueError(f"{name} contains out-of-range indices for natoms={natoms}")


def _build_z_bin_geometry(atoms: Atoms, dz_A: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    if dz_A <= 0.0:
        raise ValueError("dz_A must be > 0")

    cell = np.asarray(atoms.cell.array, dtype=float)
    area_xy_A2 = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    lz_A = float(np.linalg.norm(cell[2]))
    if area_xy_A2 <= 0.0 or lz_A <= 0.0:
        raise ValueError("cell area/height must be positive to build z bins")

    nbins = max(1, int(np.ceil(lz_A / dz_A)))
    if nbins == 1:
        z_edges_A = np.array([0.0, lz_A], dtype=float)
    else:
        interior = np.arange(1, nbins, dtype=float) * dz_A
        z_edges_A = np.concatenate(([0.0], interior, [lz_A]))

    bin_widths_A = np.diff(z_edges_A)
    bin_volumes_A3 = area_xy_A2 * bin_widths_A
    return z_edges_A, bin_widths_A, bin_volumes_A3, lz_A


def _assign_z_bins(z_values_A: np.ndarray, z_edges_A: np.ndarray) -> np.ndarray:
    nbins = z_edges_A.size - 1
    bin_ids = np.searchsorted(z_edges_A, z_values_A, side="right") - 1
    return np.clip(bin_ids, 0, nbins - 1).astype(int, copy=False)


def _fraction_window_mask(values: np.ndarray, *, start: float, end: float) -> np.ndarray:
    """
    Build mask for a periodic fractional interval on [0, 1).

    The interval is left-closed and right-open. If start == end (mod 1), treat as full range.
    """
    values = np.asarray(values, dtype=float)
    start_mod = float(np.mod(start, 1.0))
    end_mod = float(np.mod(end, 1.0))

    if np.isclose(start_mod, end_mod):
        return np.ones(values.shape, dtype=bool)
    if start_mod < end_mod:
        return (values >= start_mod) & (values < end_mod)
    return (values >= start_mod) | (values < end_mod)


def _theta_bin_count_from_ndeg(ndeg: float) -> int:
    ndeg = float(ndeg)
    if ndeg <= 0.0:
        raise ValueError("ndeg must be > 0")
    n_bins_float = 180.0 / ndeg
    n_bins = int(round(n_bins_float))
    if not np.isclose(n_bins_float, float(n_bins), rtol=0.0, atol=1.0e-12):
        raise ValueError("ndeg must divide 180 exactly")
    return n_bins


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


def detect_water_molecule_indices(
    atoms: Atoms,
    *,
    oxygen_symbol: str = "O",
    hydrogen_symbol: str = "H",
    oh_cutoff_A: float = DEFAULT_WATER_OH_CUTOFF_A,
) -> np.ndarray:
    """
    Label water molecules as an (n_water, 3) array [O_idx, H1_idx, H2_idx].

    Connectivity is inferred from O-H distances under MIC with a distance cutoff.
    Each H can belong to at most one O in the final assignment.
    """
    if oh_cutoff_A <= 0.0:
        raise ValueError("oh_cutoff_A must be > 0")

    symbols = np.asarray(atoms.get_chemical_symbols())
    oxygen_indices = np.where(symbols == oxygen_symbol)[0]
    hydrogen_indices = np.where(symbols == hydrogen_symbol)[0]

    if oxygen_indices.size == 0:
        raise WaterTopologyError(f"No oxygen atoms found with symbol={oxygen_symbol!r}")
    if hydrogen_indices.size < 2:
        raise WaterTopologyError(f"Not enough hydrogen atoms found with symbol={hydrogen_symbol!r}")

    # Build candidate O-H assignments under cutoff, then greedily assign shortest pairs.
    candidate_pairs: list[tuple[float, int, int]] = []
    hydrogen_list = hydrogen_indices.tolist()
    for local_o, o_idx in enumerate(oxygen_indices):
        distances = np.asarray(atoms.get_distances(int(o_idx), hydrogen_list, mic=True), dtype=float)
        for local_h, distance_A in enumerate(distances):
            if distance_A <= oh_cutoff_A:
                candidate_pairs.append((float(distance_A), int(local_o), int(hydrogen_indices[local_h])))

    if not candidate_pairs:
        raise WaterTopologyError("No O-H pairs found under the given oh_cutoff_A")

    candidate_pairs.sort(key=lambda x: x[0])
    assigned_h: set[int] = set()
    h_by_local_o: list[list[int]] = [[] for _ in range(oxygen_indices.size)]
    for _, local_o, h_idx in candidate_pairs:
        if len(h_by_local_o[local_o]) >= 2:
            continue
        if h_idx in assigned_h:
            continue
        h_by_local_o[local_o].append(h_idx)
        assigned_h.add(h_idx)

    water_rows: list[tuple[int, int, int]] = []
    for local_o, h_list in enumerate(h_by_local_o):
        if len(h_list) != 2:
            continue
        o_idx = int(oxygen_indices[local_o])
        h1_idx, h2_idx = sorted(h_list)
        water_rows.append((o_idx, h1_idx, h2_idx))

    if not water_rows:
        raise WaterTopologyError("No complete H2O molecules could be identified")

    return np.asarray(water_rows, dtype=int)


def get_water_oxygen_indices_array(water_molecule_indices: np.ndarray) -> np.ndarray:
    """
    Return water oxygen indices as an (n_water, 1) integer array.
    """
    water_molecule_indices = np.asarray(water_molecule_indices, dtype=int)
    if water_molecule_indices.ndim != 2 or water_molecule_indices.shape[1] != 3:
        raise ValueError("water_molecule_indices must have shape (n_water, 3)")
    return water_molecule_indices[:, 0:1]


def compute_water_mass_density_z_distribution(
    atoms: Atoms,
    oxygen_indices: Sequence[int] | np.ndarray,
    *,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
) -> np.ndarray:
    """
    Compute water mass-density distribution along z, shape (nbins, 1), unit g/cm^3.

    Each oxygen is treated as one H2O molecule in the corresponding z bin.
    """
    oxygen_indices = _as_flat_index_array(oxygen_indices, name="oxygen_indices")
    _validate_atom_indices(oxygen_indices, len(atoms), name="oxygen_indices")

    z_edges_A, _, bin_volumes_A3, lz_A = _build_z_bin_geometry(atoms, dz_A=dz_A)

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    z_all_A = scaled[:, 2] * lz_A
    z_oxygen_A = z_all_A[oxygen_indices]
    bin_ids = _assign_z_bins(z_oxygen_A, z_edges_A)

    counts_per_bin = np.bincount(bin_ids, minlength=bin_volumes_A3.size).astype(float)
    mass_per_water_g = WATER_MOLAR_MASS_G_PER_MOL / AVOGADRO_NUMBER
    mass_per_bin_g = counts_per_bin * mass_per_water_g

    density_g_cm3 = mass_per_bin_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)
    return density_g_cm3.reshape(-1, 1)


def compute_water_orientation_weighted_density_z_distribution(
    atoms: Atoms,
    oxygen_indices: Sequence[int] | np.ndarray,
    *,
    water_molecule_indices: np.ndarray | None = None,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
) -> np.ndarray:
    """
    Compute plane-density-weighted orientation distribution along z, shape (nbins, 1).

    For each oxygen, theta is the angle between +c (cell c-axis direction) and the H-O-H bisector direction
    (vector starts at O and points along the molecular angle bisector). The profile
    value per bin is sum(cos(theta_i)) / V_bin, i.e. density-weighted orientation.
    """
    oxygen_indices = _as_flat_index_array(oxygen_indices, name="oxygen_indices")
    _validate_atom_indices(oxygen_indices, len(atoms), name="oxygen_indices")

    if water_molecule_indices is None:
        water_molecule_indices = detect_water_molecule_indices(atoms)
    o_to_h = _oxygen_to_hydrogen_map(water_molecule_indices)

    missing_o = [int(o_idx) for o_idx in oxygen_indices if int(o_idx) not in o_to_h]
    if missing_o:
        raise WaterTopologyError(
            "Some oxygen indices are not labeled as water molecules: "
            f"{missing_o[:8]}{'...' if len(missing_o) > 8 else ''}"
        )

    z_edges_A, _, bin_volumes_A3, lz_A = _build_z_bin_geometry(atoms, dz_A=dz_A)
    cell = np.asarray(atoms.cell.array, dtype=float)
    c_vec = cell[2]
    c_norm = float(np.linalg.norm(c_vec))
    if c_norm <= 0.0:
        raise ValueError("cell c-axis norm must be positive for orientation projection")
    c_unit = c_vec / c_norm

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    z_all_A = scaled[:, 2] * lz_A
    z_oxygen_A = z_all_A[oxygen_indices]
    bin_ids = _assign_z_bins(z_oxygen_A, z_edges_A)

    cos_theta_values = np.empty(oxygen_indices.size, dtype=float)
    for i, o_idx in enumerate(oxygen_indices):
        h1_idx, h2_idx = o_to_h[int(o_idx)]
        vecs = np.asarray(atoms.get_distances(int(o_idx), [int(h1_idx), int(h2_idx)], vector=True, mic=True), dtype=float)
        v1 = vecs[0]
        v2 = vecs[1]

        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 == 0.0 or n2 == 0.0:
            raise WaterTopologyError(f"Zero O-H vector norm for oxygen index {int(o_idx)}")

        bisector = v1 / n1 + v2 / n2
        nb = float(np.linalg.norm(bisector))
        if nb == 0.0:
            raise WaterTopologyError(f"Degenerate H-O-H bisector for oxygen index {int(o_idx)}")

        cos_theta_values[i] = float(np.dot(bisector, c_unit) / nb)

    cos_sum_per_bin = np.bincount(bin_ids, weights=cos_theta_values, minlength=bin_volumes_A3.size).astype(float)
    weighted_orientation_density = cos_sum_per_bin / bin_volumes_A3
    return weighted_orientation_density.reshape(-1, 1)


def compute_water_orientation_theta_pdf_in_c_fraction_window(
    atoms: Atoms,
    oxygen_indices: Sequence[int] | np.ndarray,
    c_fraction_range: Sequence[float],
    *,
    water_molecule_indices: np.ndarray | None = None,
    ndeg: float = DEFAULT_THETA_BIN_DEG,
) -> np.ndarray:
    """
    Compute orientation-theta PDF in a given c-fraction window, shape (180 / ndeg,).

    The angle theta is between the H-O-H bisector and +c direction of the cell.
    Returns probability density in degree^-1. The histogram domain is [0, 180] degree.
    """
    oxygen_indices = _as_flat_index_array(oxygen_indices, name="oxygen_indices")
    _validate_atom_indices(oxygen_indices, len(atoms), name="oxygen_indices")

    if len(c_fraction_range) != 2:
        raise ValueError("c_fraction_range must contain exactly two values: [start, end]")
    start_c_fraction = float(c_fraction_range[0])
    end_c_fraction = float(c_fraction_range[1])

    n_theta_bins = _theta_bin_count_from_ndeg(float(ndeg))

    if water_molecule_indices is None:
        water_molecule_indices = detect_water_molecule_indices(atoms)
    o_to_h = _oxygen_to_hydrogen_map(water_molecule_indices)

    missing_o = [int(o_idx) for o_idx in oxygen_indices if int(o_idx) not in o_to_h]
    if missing_o:
        raise WaterTopologyError(
            "Some oxygen indices are not labeled as water molecules: "
            f"{missing_o[:8]}{'...' if len(missing_o) > 8 else ''}"
        )

    cell = np.asarray(atoms.cell.array, dtype=float)
    c_vec = cell[2]
    c_norm = float(np.linalg.norm(c_vec))
    if c_norm <= 0.0:
        raise ValueError("cell c-axis norm must be positive for orientation projection")
    c_unit = c_vec / c_norm

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    oxygen_frac_c = scaled[oxygen_indices, 2]
    in_window = _fraction_window_mask(oxygen_frac_c, start=start_c_fraction, end=end_c_fraction)
    selected_oxygen = oxygen_indices[in_window]

    if selected_oxygen.size == 0:
        return np.zeros(n_theta_bins, dtype=float)

    theta_deg_values = np.empty(selected_oxygen.size, dtype=float)
    for i, o_idx in enumerate(selected_oxygen):
        h1_idx, h2_idx = o_to_h[int(o_idx)]
        vecs = np.asarray(atoms.get_distances(int(o_idx), [int(h1_idx), int(h2_idx)], vector=True, mic=True), dtype=float)
        v1 = vecs[0]
        v2 = vecs[1]

        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 == 0.0 or n2 == 0.0:
            raise WaterTopologyError(f"Zero O-H vector norm for oxygen index {int(o_idx)}")

        bisector = v1 / n1 + v2 / n2
        nb = float(np.linalg.norm(bisector))
        if nb == 0.0:
            raise WaterTopologyError(f"Degenerate H-O-H bisector for oxygen index {int(o_idx)}")

        cos_theta = float(np.dot(bisector, c_unit) / nb)
        cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
        theta_deg_values[i] = float(np.degrees(np.arccos(cos_theta)))

    theta_edges_deg = np.linspace(0.0, 180.0, n_theta_bins + 1, dtype=float)
    theta_pdf, _ = np.histogram(theta_deg_values, bins=theta_edges_deg, density=True)
    return theta_pdf.astype(float, copy=False)
