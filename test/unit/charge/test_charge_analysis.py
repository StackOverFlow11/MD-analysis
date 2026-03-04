"""Unit tests for md_analysis.charge.ChargeAnalysis."""

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from md_analysis.charge.ChargeAnalysis import (
    AtomSelector,
    ElementSelector,
    IndexSelector,
    TrajectoryChargeResult,
    compute_frame_surface_charge,
)
from md_analysis.charge.config import E_PER_A2_TO_UC_PER_CM2
from md_analysis.utils.BaderParser import load_bader_atoms

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "bader_work_dir"


# ---------------------------------------------------------------------------
# Atom selectors
# ---------------------------------------------------------------------------

class TestAtomSelectors:
    """Tests for ElementSelector and IndexSelector."""

    def test_element_selector_basic(self):
        atoms = Atoms("CuAgOH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        sel = ElementSelector({"Ag"})
        idx = sel.select(atoms)
        np.testing.assert_array_equal(idx, [1])

    def test_element_selector_multiple(self):
        atoms = Atoms("CuAgOH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        sel = ElementSelector({"Cu", "Ag"})
        idx = sel.select(atoms)
        np.testing.assert_array_equal(idx, [0, 1])

    def test_element_selector_no_match(self):
        atoms = Atoms("CuCu", positions=[[0, 0, 0], [1, 0, 0]])
        sel = ElementSelector({"Ag"})
        idx = sel.select(atoms)
        assert len(idx) == 0

    def test_element_selector_empty_symbols_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            ElementSelector(set())

    def test_index_selector_basic(self):
        atoms = Atoms("CuAgO", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        sel = IndexSelector([0, 2])
        idx = sel.select(atoms)
        np.testing.assert_array_equal(idx, [0, 2])

    def test_index_selector_out_of_range(self):
        atoms = Atoms("CuAg", positions=[[0, 0, 0], [1, 0, 0]])
        sel = IndexSelector([0, 5, 10])
        idx = sel.select(atoms)
        np.testing.assert_array_equal(idx, [0])

    def test_index_selector_deduplicates(self):
        atoms = Atoms("CuAgO", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        sel = IndexSelector([1, 1, 2])
        idx = sel.select(atoms)
        np.testing.assert_array_equal(idx, [1, 2])

    def test_abc_not_instantiable(self):
        with pytest.raises(TypeError):
            AtomSelector()  # type: ignore[abstract]


# ---------------------------------------------------------------------------
# compute_frame_surface_charge
# ---------------------------------------------------------------------------

class TestComputeFrameSurfaceCharge:
    """Tests for single-frame surface charge computation."""

    @pytest.fixture()
    def bader_atoms(self):
        return load_bader_atoms(
            DATA_DIR / "POSCAR", DATA_DIR / "ACF.dat", DATA_DIR / "POTCAR"
        )

    def test_basic_sigma(self, bader_atoms):
        result = compute_frame_surface_charge(bader_atoms)
        assert "surface_charge_density_e_A2" in result.info
        assert "surface_charge_density_uC_cm2" in result.info
        sigma_e = result.info["surface_charge_density_e_A2"]
        sigma_uc = result.info["surface_charge_density_uC_cm2"]
        assert len(sigma_e) == 2
        assert len(sigma_uc) == 2
        # Verify unit conversion
        for i in range(2):
            assert sigma_uc[i] == pytest.approx(
                sigma_e[i] * E_PER_A2_TO_UC_PER_CM2, rel=1e-10
            )

    def test_sigma_is_finite(self, bader_atoms):
        result = compute_frame_surface_charge(bader_atoms)
        sigma = result.info["surface_charge_density_uC_cm2"]
        assert all(np.isfinite(v) for v in sigma)

    def test_with_element_selector(self, bader_atoms):
        sel = ElementSelector({"Ag"})
        result = compute_frame_surface_charge(bader_atoms, atom_selector=sel)
        assert "selected_atom_indices" in result.info
        assert "selected_atom_net_charges" in result.info
        # There are 2 Ag atoms in the structure
        assert len(result.info["selected_atom_indices"]) == 2
        assert len(result.info["selected_atom_net_charges"]) == 2

    def test_without_selector_no_selected_keys(self, bader_atoms):
        result = compute_frame_surface_charge(bader_atoms)
        assert "selected_atom_indices" not in result.info
        assert "selected_atom_net_charges" not in result.info

    def test_missing_bader_net_charge_raises(self):
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        with pytest.raises(ValueError, match="bader_net_charge"):
            compute_frame_surface_charge(atoms)


# ---------------------------------------------------------------------------
# TrajectoryChargeResult dataclass
# ---------------------------------------------------------------------------

class TestTrajectoryChargeResult:
    """Tests for the frozen dataclass."""

    def test_frozen(self):
        result = TrajectoryChargeResult(
            frame_labels=("f1", "f2"),
            surface_charge_density_uC_cm2=np.zeros((2, 2)),
            mean_surface_charge_density_uC_cm2=np.zeros(2),
            std_surface_charge_density_uC_cm2=np.zeros(2),
            selected_atom_net_charges=np.empty((2, 0)),
            mean_selected_atom_net_charges=np.empty(0),
            selected_atom_indices=(),
        )
        with pytest.raises(AttributeError):
            result.frame_labels = ("x",)  # type: ignore[misc]

    def test_fields_accessible(self):
        sigma = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = TrajectoryChargeResult(
            frame_labels=("a", "b"),
            surface_charge_density_uC_cm2=sigma,
            mean_surface_charge_density_uC_cm2=sigma.mean(axis=0),
            std_surface_charge_density_uC_cm2=sigma.std(axis=0),
            selected_atom_net_charges=np.empty((2, 0)),
            mean_selected_atom_net_charges=np.empty(0),
            selected_atom_indices=(),
        )
        assert result.frame_labels == ("a", "b")
        assert result.surface_charge_density_uC_cm2.shape == (2, 2)
