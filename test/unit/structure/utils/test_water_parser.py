"""Unit tests for water parsing and z-distribution utilities."""

from __future__ import annotations

import math

import numpy as np
import pytest

ase = pytest.importorskip("ase")
from ase import Atoms

from src.structure.utils.WaterParser import _compute_water_mass_density_z_distribution as compute_water_mass_density_z_distribution
from src.structure.utils.WaterParser import _compute_water_orientation_theta_pdf_in_c_fraction_window as compute_water_orientation_theta_pdf_in_c_fraction_window
from src.structure.utils.WaterParser import _compute_water_orientation_weighted_density_z_distribution as compute_water_orientation_weighted_density_z_distribution
from src.structure.utils.WaterParser import detect_water_molecule_indices
from src.structure.utils.WaterParser import get_water_oxygen_indices_array
from src.structure.utils.config import DEFAULT_THETA_BIN_DEG
from src.structure.utils.config import DEFAULT_Z_BIN_WIDTH_A


@pytest.fixture
def two_water_atoms() -> Atoms:
    """
    Build a tiny periodic box with two isolated water molecules.
    """
    symbols = ["O", "H", "H", "O", "H", "H"]
    positions = np.array(
        [
            [2.00, 2.00, 2.00],  # O1
            [2.96, 2.00, 2.00],  # H1
            [1.76, 2.93, 2.00],  # H2
            [7.00, 7.00, 7.00],  # O2
            [7.96, 7.00, 7.00],  # H3
            [6.76, 7.93, 7.00],  # H4
        ],
        dtype=float,
    )
    return Atoms(symbols=symbols, positions=positions, cell=[10.0, 10.0, 10.0], pbc=[True, True, True])


def test_detect_water_molecule_indices_shape(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    assert water_idx.ndim == 2
    assert water_idx.shape == (2, 3)


def test_get_water_oxygen_indices_array_shape(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    assert oxygen_n1.shape == (2, 1)
    assert set(oxygen_n1.reshape(-1).tolist()) == {0, 3}


def test_compute_water_mass_density_z_distribution_shape(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    rho = compute_water_mass_density_z_distribution(two_water_atoms, oxygen_n1)
    expected_nbins = math.ceil(10.0 / DEFAULT_Z_BIN_WIDTH_A)
    assert rho.ndim == 2
    assert rho.shape == (expected_nbins, 1)
    assert np.all(np.isfinite(rho))
    assert float(np.sum(rho)) > 0.0


def test_compute_water_orientation_weighted_density_z_distribution_shape(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    orient = compute_water_orientation_weighted_density_z_distribution(
        two_water_atoms,
        oxygen_n1,
        water_molecule_indices=water_idx,
    )
    expected_nbins = math.ceil(10.0 / DEFAULT_Z_BIN_WIDTH_A)
    assert orient.ndim == 2
    assert orient.shape == (expected_nbins, 1)
    assert np.all(np.isfinite(orient))


def test_compute_water_orientation_theta_pdf_in_c_fraction_window_shape_and_norm(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    pdf = compute_water_orientation_theta_pdf_in_c_fraction_window(
        two_water_atoms,
        oxygen_n1,
        c_fraction_range=[0.0, 1.0],
        water_molecule_indices=water_idx,
    )

    expected_bins = int(round(180.0 / DEFAULT_THETA_BIN_DEG))
    assert pdf.shape == (expected_bins,)
    assert np.all(np.isfinite(pdf))
    assert np.all(pdf >= 0.0)
    assert np.isclose(float(np.sum(pdf) * DEFAULT_THETA_BIN_DEG), 1.0, rtol=1.0e-12, atol=1.0e-12)


def test_compute_water_orientation_theta_pdf_in_c_fraction_window_empty_returns_zero(two_water_atoms: Atoms) -> None:
    water_idx = detect_water_molecule_indices(two_water_atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    pdf = compute_water_orientation_theta_pdf_in_c_fraction_window(
        two_water_atoms,
        oxygen_n1,
        c_fraction_range=[0.30, 0.40],
        water_molecule_indices=water_idx,
    )

    expected_bins = int(round(180.0 / DEFAULT_THETA_BIN_DEG))
    assert pdf.shape == (expected_bins,)
    assert np.allclose(pdf, 0.0)
