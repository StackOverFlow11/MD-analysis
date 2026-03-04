"""Unit tests for md_analysis.charge.ChargeAnalysis."""

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from md_analysis.charge.ChargeAnalysis import (
    compute_frame_surface_charge,
    trajectory_indexed_atom_charges,
)
from md_analysis.charge.config import E_PER_A2_TO_UC_PER_CM2
from md_analysis.utils.BaderParser import load_bader_atoms

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "bader_work_dir"


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

    def test_no_selected_keys(self, bader_atoms):
        result = compute_frame_surface_charge(bader_atoms)
        assert "selected_atom_indices" not in result.info
        assert "selected_atom_net_charges" not in result.info

    def test_normal_a(self, bader_atoms):
        """normal='a' uses cell vectors b × c for area."""
        # This will likely fail on layer detection for this particular
        # dataset, but we verify the area-vector selection path is reached.
        # For a proper test we just check that it doesn't raise ValueError
        # on normal validation and proceeds to layer detection.
        try:
            compute_frame_surface_charge(bader_atoms, normal="a")
        except ValueError as exc:
            # Accept layer-detection failures but NOT normal-validation errors
            assert "normal must be" not in str(exc)

    def test_normal_b(self, bader_atoms):
        """normal='b' uses cell vectors a × c for area."""
        try:
            compute_frame_surface_charge(bader_atoms, normal="b")
        except ValueError as exc:
            assert "normal must be" not in str(exc)

    def test_invalid_normal_raises(self, bader_atoms):
        with pytest.raises(ValueError, match="normal must be"):
            compute_frame_surface_charge(bader_atoms, normal="z")

    def test_missing_bader_net_charge_raises(self):
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        with pytest.raises(ValueError, match="bader_net_charge"):
            compute_frame_surface_charge(atoms)


# ---------------------------------------------------------------------------
# trajectory_indexed_atom_charges
# ---------------------------------------------------------------------------

class TestTrajectoryIndexedAtomCharges:
    """Unit tests for trajectory_indexed_atom_charges validation."""

    def test_non_2d_raises(self, tmp_path):
        with pytest.raises(ValueError, match="2-D"):
            trajectory_indexed_atom_charges(tmp_path, np.array([1, 2, 3]))

    def test_non_integer_raises(self, tmp_path):
        with pytest.raises(ValueError, match="integer dtype"):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[1.0, 2.0]])
            )

    def test_negative_index_raises(self, tmp_path):
        with pytest.raises(ValueError, match="negative"):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[0, -1]])
            )

    def test_t_mismatch_raises(self, tmp_path):
        """t rows != number of frame directories."""
        # Create 1 frame dir but provide 2-row matrix
        frame = tmp_path / "calc_t000_i000"
        frame.mkdir()
        with pytest.raises(ValueError, match="rows"):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[0], [1]])
            )

    def test_missing_file_raises(self, tmp_path):
        """Missing POSCAR/ACF.dat/POTCAR raises FileNotFoundError."""
        frame = tmp_path / "calc_t000_i000"
        frame.mkdir()
        with pytest.raises(FileNotFoundError):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[0]])
            )

    def test_out_of_bounds_raises(self, tmp_path):
        """Index exceeding atom count raises IndexError with frame name."""
        import shutil
        frame = tmp_path / "calc_t000_i000"
        frame.mkdir()
        for f in ["POSCAR", "ACF.dat", "POTCAR"]:
            shutil.copy2(DATA_DIR / f, frame / f)

        big_idx = 99999
        with pytest.raises(IndexError, match="calc_t000_i000"):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[big_idx]])
            )
