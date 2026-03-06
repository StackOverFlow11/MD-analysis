"""Unit tests for md_analysis.charge.BaderAnalysis."""

import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from md_analysis.charge.BaderAnalysis import (
    _extract_t_value,
    _sorted_frame_dirs,
    compute_frame_surface_charge,
    frame_indexed_atom_charges,
    surface_charge_analysis,
    trajectory_indexed_atom_charges,
    trajectory_surface_charge,
)
from md_analysis.charge.config import E_PER_A2_TO_UC_PER_CM2
from md_analysis.utils.BaderParser import load_bader_atoms
from md_analysis.utils.LayerParser import detect_interface_layers

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
        """Test that σ results are present and well-formed."""
        result = compute_frame_surface_charge(bader_atoms)
        assert "surface_charge_density_e_A2" in result.info
        assert "surface_charge_density_uC_cm2" in result.info
        sigma_e = result.info["surface_charge_density_e_A2"]
        sigma_uc = result.info["surface_charge_density_uC_cm2"]
        assert len(sigma_e) == 2
        assert len(sigma_uc) == 2
        assert all(np.isfinite(v) for v in sigma_e)
        assert all(np.isfinite(v) for v in sigma_uc)

    def test_sigma_is_finite(self, bader_atoms):
        result = compute_frame_surface_charge(bader_atoms)
        sigma = result.info["surface_charge_density_uC_cm2"]
        assert all(np.isfinite(v) for v in sigma)

    def test_info_keys_present(self, bader_atoms):
        """All expected info keys are set after computation."""
        result = compute_frame_surface_charge(bader_atoms)
        assert "surface_charge_density_uC_cm2" in result.info
        assert "n_charged_atoms_per_surface" in result.info
        assert "charge_per_surface_e" in result.info
        n_ch = result.info["n_charged_atoms_per_surface"]
        assert len(n_ch) == 2
        assert all(isinstance(v, int) for v in n_ch)

    def test_with_synthetic_charged_atoms(self):
        """Build a mock slab + charged atoms → σ reflects counterion charges only."""
        # 4-layer Cu slab along c with vacuum
        cell_c = 30.0
        cell_a = 5.0
        # Place 4 Cu layers at frac c = 0.1, 0.2, 0.3, 0.4 (3–12 Å)
        n_per_layer = 4
        positions = []
        symbols = []
        for layer_frac in [0.1, 0.2, 0.3, 0.4]:
            z = layer_frac * cell_c
            for ix in range(2):
                for iy in range(2):
                    positions.append([ix * cell_a / 2, iy * cell_a / 2, z])
                    symbols.append("Cu")
        # Add one O + two H (water) at frac c = 0.6 (18 Å, in water region)
        positions.append([2.5, 2.5, 0.6 * cell_c])
        symbols.append("O")
        positions.append([2.5, 3.2, 0.6 * cell_c + 0.5])
        symbols.append("H")
        positions.append([2.5, 1.8, 0.6 * cell_c + 0.5])
        symbols.append("H")
        # Add one K near bottom surface (frac c = 0.05)
        positions.append([2.5, 2.5, 0.05 * cell_c])
        symbols.append("K")

        atoms = Atoms(
            symbols=symbols,
            positions=positions,
            cell=[cell_a, cell_a, cell_c],
            pbc=True,
        )
        # Assign fake Bader net charges
        net_charge = np.zeros(len(atoms))
        # Cu atoms: small positive charge
        for i in range(n_per_layer * 4):
            net_charge[i] = 0.01
        # O: slightly negative
        net_charge[-4] = -0.1
        # H: slightly positive
        net_charge[-3] = 0.05
        net_charge[-2] = 0.05
        # K: +0.8e
        net_charge[-1] = 0.8
        atoms.arrays["bader_net_charge"] = net_charge

        result = compute_frame_surface_charge(
            atoms, metal_symbols={"Cu"}, normal="c",
        )

        sigma = result.info["surface_charge_density_uC_cm2"]
        n_ch = result.info["n_charged_atoms_per_surface"]
        q_ch = result.info["charge_per_surface_e"]
        # K is near bottom surface → contributes to σ_bottom
        assert n_ch[0] == 1
        assert q_ch[0] == pytest.approx(0.8)
        assert sigma[0] != 0.0
        # Water and metal are excluded → no atoms assigned to top
        assert n_ch[1] == 0
        assert q_ch[1] == 0.0
        # Charges are finite
        assert all(np.isfinite(v) for v in sigma)
        assert all(np.isfinite(v) for v in q_ch)

    def test_layer_method_is_invariant_under_periodic_shift(self, bader_atoms):
        """Periodic translation must not swap the two surface-charge columns."""
        det = detect_interface_layers(bader_atoms, normal="c")
        net_charge = bader_atoms.arrays["bader_net_charge"].copy()

        # Charge the two detected interface layers asymmetrically so a column
        # swap becomes visible after a periodic translation.
        net_charge[list(det.interface_normal_aligned().atom_indices)] += 1.0
        net_charge[list(det.interface_normal_opposed().atom_indices)] -= 1.0
        bader_atoms.arrays["bader_net_charge"] = net_charge

        result = compute_frame_surface_charge(
            bader_atoms, normal="c", method="layer",
        )
        sigma_ref = np.array(result.info["surface_charge_density_e_A2"], dtype=float)

        shifted = bader_atoms.copy()
        shifted.translate(np.array([0.0, 0.0, 0.45]) @ shifted.cell.array)
        shifted.wrap()
        shifted_result = compute_frame_surface_charge(
            shifted, normal="c", method="layer",
        )
        sigma_shifted = np.array(
            shifted_result.info["surface_charge_density_e_A2"], dtype=float
        )

        np.testing.assert_allclose(sigma_shifted, sigma_ref, atol=1e-12)

    def test_normal_a(self, bader_atoms):
        """normal='a' uses cell vectors b × c for area."""
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


# ---------------------------------------------------------------------------
# frame_indexed_atom_charges
# ---------------------------------------------------------------------------

class TestFrameIndexedAtomCharges:
    """Tests for single-frame indexed atom charge extraction."""

    @pytest.fixture()
    def bader_atoms(self):
        return load_bader_atoms(
            DATA_DIR / "POSCAR", DATA_DIR / "ACF.dat", DATA_DIR / "POTCAR"
        )

    def test_basic_shape_and_values(self, bader_atoms):
        idx = np.array([0, 1])
        result = frame_indexed_atom_charges(bader_atoms, idx)
        assert result.shape == (2, 2)
        np.testing.assert_array_equal(result[:, 0], idx)

    def test_missing_bader_raises(self):
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        with pytest.raises(ValueError, match="bader_net_charge"):
            frame_indexed_atom_charges(atoms, np.array([0]))

    def test_non_1d_raises(self, bader_atoms):
        with pytest.raises(ValueError, match="1-D"):
            frame_indexed_atom_charges(bader_atoms, np.array([[0, 1]]))

    def test_non_integer_raises(self, bader_atoms):
        with pytest.raises(ValueError, match="integer dtype"):
            frame_indexed_atom_charges(bader_atoms, np.array([0.0, 1.0]))

    def test_negative_raises(self, bader_atoms):
        with pytest.raises(ValueError, match="negative"):
            frame_indexed_atom_charges(bader_atoms, np.array([0, -1]))

    def test_out_of_bounds_raises(self, bader_atoms):
        with pytest.raises(IndexError, match="out of bounds"):
            frame_indexed_atom_charges(bader_atoms, np.array([99999]))


# ---------------------------------------------------------------------------
# trajectory_surface_charge (unit-level error paths)
# ---------------------------------------------------------------------------

class TestTrajectorySurfaceCharge:
    """Unit tests for trajectory_surface_charge validation."""

    def test_missing_dir_raises(self, tmp_path):
        fake = tmp_path / "nonexistent"
        with pytest.raises(FileNotFoundError, match="root_dir does not exist"):
            trajectory_surface_charge(fake)

    def test_invalid_normal_raises(self, tmp_path):
        with pytest.raises(ValueError, match="normal must be"):
            trajectory_surface_charge(tmp_path, normal="z")


# ---------------------------------------------------------------------------
# _extract_t_value / _sorted_frame_dirs
# ---------------------------------------------------------------------------

class TestSortedFrameDirs:
    """Tests for numeric sorting of frame directories."""

    def test_extract_t_value_basic(self):
        assert _extract_t_value("calc_t50_i0") == 50
        assert _extract_t_value("calc_t1000_i0") == 1000
        assert _extract_t_value("calc_t0_i0") == 0

    def test_extract_t_value_no_match(self):
        assert _extract_t_value("no_match") == 0

    def test_sorted_frame_dirs_numeric_order(self, tmp_path):
        """Non-zero-padded directory names sort numerically, not lexically."""
        for t in [50, 1000, 200, 5]:
            (tmp_path / f"calc_t{t}_i0").mkdir()
        result = _sorted_frame_dirs(tmp_path, "calc_t*_i*")
        names = [p.name for p in result]
        assert names == ["calc_t5_i0", "calc_t50_i0", "calc_t200_i0", "calc_t1000_i0"]

    def test_sorted_frame_dirs_empty_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No subdirectories"):
            _sorted_frame_dirs(tmp_path, "calc_t*_i*")


# ---------------------------------------------------------------------------
# surface_charge_analysis
# ---------------------------------------------------------------------------

_FRAME_FILES = ["POSCAR", "ACF.dat", "POTCAR"]


def _build_fake_trajectory(tmp_path: Path, n_frames: int = 2) -> Path:
    for i in range(n_frames):
        frame_dir = tmp_path / f"calc_t{i:03d}_i000"
        frame_dir.mkdir()
        for fname in _FRAME_FILES:
            shutil.copy2(DATA_DIR / fname, frame_dir / fname)
    return tmp_path


class TestSurfaceChargeAnalysis:
    """Unit tests for surface_charge_analysis."""

    def test_csv_and_png_created(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        out = tmp_path / "output"
        csv_path = surface_charge_analysis(root, output_dir=out)
        assert csv_path.exists()
        png_path = out / "surface_charge.png"
        assert png_path.exists()

    def test_csv_columns(self, tmp_path):
        import csv
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        out = tmp_path / "output"
        csv_path = surface_charge_analysis(root, output_dir=out)
        with csv_path.open(encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 2
        expected_cols = {
            "step", "sigma_bottom_uC_cm2", "sigma_top_uC_cm2",
            "sigma_bottom_cumavg_uC_cm2", "sigma_top_cumavg_uC_cm2",
        }
        assert set(rows[0].keys()) == expected_cols

    def test_frame_slicing(self, tmp_path):
        import csv
        root = _build_fake_trajectory(tmp_path, n_frames=4)
        out = tmp_path / "output"
        csv_path = surface_charge_analysis(
            root, output_dir=out, frame_start=1, frame_end=3,
        )
        with csv_path.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 2

    def test_invalid_normal_raises(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=1)
        with pytest.raises(ValueError, match="normal must be"):
            surface_charge_analysis(root, normal="z")
