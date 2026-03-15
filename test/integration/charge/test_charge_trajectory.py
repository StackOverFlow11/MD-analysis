"""Integration tests for trajectory charge functions."""

import csv
import shutil
from pathlib import Path

import numpy as np
import pytest

from md_analysis.electrochemical.charge.Bader.AtomCharges import (
    trajectory_indexed_atom_charges,
)
from md_analysis.electrochemical.charge.Bader.SurfaceCharge import (
    compute_frame_surface_charge,
    surface_charge_analysis,
    trajectory_surface_charge,
)
from md_analysis.utils.BaderParser import load_bader_atoms

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "bader" / "bader_work_dir"

# Files needed per frame
_FRAME_FILES = ["POSCAR", "ACF.dat", "POTCAR"]


def _build_fake_trajectory(tmp_path: Path, n_frames: int = 2) -> Path:
    """Copy bader_work_dir data into bader_t*_i* subdirectories."""
    for i in range(n_frames):
        frame_dir = tmp_path / f"bader_t{i:03d}_i000"
        frame_dir.mkdir()
        for fname in _FRAME_FILES:
            shutil.copy2(DATA_DIR / fname, frame_dir / fname)
    return tmp_path


class TestTrajectoryIndexedAtomCharges:
    """Integration tests using real Bader data arranged as a fake trajectory."""

    def test_basic_shape_and_values(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        # Query atoms 0 and 1 in both frames
        idx_matrix = np.array([[0, 1], [0, 1]])
        result = trajectory_indexed_atom_charges(root, idx_matrix)

        assert result.shape == (2, 2, 2)
        # Channel 0 echoes back input indices
        np.testing.assert_array_equal(result[:, :, 0], idx_matrix)
        # Channel 1 contains net charges (finite floats)
        assert np.all(np.isfinite(result[:, :, 1]))

    def test_indices_echoed(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        idx_matrix = np.array([[3, 5, 7], [3, 5, 7]])
        result = trajectory_indexed_atom_charges(root, idx_matrix)

        np.testing.assert_array_equal(result[:, :, 0], idx_matrix)

    def test_identical_frames_same_charges(self, tmp_path):
        """All frames are identical copies, so charges should match."""
        root = _build_fake_trajectory(tmp_path, n_frames=3)
        idx_matrix = np.array([[0, 1]] * 3)
        result = trajectory_indexed_atom_charges(root, idx_matrix)

        # Charges across frames should be identical
        np.testing.assert_allclose(
            result[0, :, 1], result[1, :, 1], atol=1e-10
        )
        np.testing.assert_allclose(
            result[0, :, 1], result[2, :, 1], atol=1e-10
        )

    def test_missing_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No subdirectories"):
            trajectory_indexed_atom_charges(
                tmp_path, np.array([[0]])
            )


class TestTrajectorySurfaceCharge:
    """Integration tests for trajectory_surface_charge."""

    def test_basic_shape_and_values(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        result = trajectory_surface_charge(root)

        assert result.shape == (2, 2)
        assert np.all(np.isfinite(result))

    def test_identical_frames_same_sigma(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        result = trajectory_surface_charge(root)

        np.testing.assert_allclose(result[0], result[1], atol=1e-10)

    def test_consistent_with_single_frame(self, tmp_path):
        """trajectory_surface_charge matches per-frame compute_frame_surface_charge."""
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        traj_result = trajectory_surface_charge(root)

        # Manually compute for each frame
        for i, frame_dir in enumerate(sorted(root.glob("bader_t*_i*"))):
            atoms = load_bader_atoms(
                frame_dir / "POSCAR",
                frame_dir / "ACF.dat",
                frame_dir / "POTCAR",
            )
            compute_frame_surface_charge(atoms)
            expected = atoms.info["surface_charge_density_uC_cm2"]
            np.testing.assert_allclose(traj_result[i], expected, atol=1e-10)

    def test_numeric_sort_non_zero_padded(self, tmp_path):
        """Non-zero-padded directory names are sorted numerically."""
        for t in [50, 1000, 200, 5]:
            frame_dir = tmp_path / f"bader_t{t}_i0"
            frame_dir.mkdir()
            for fname in _FRAME_FILES:
                shutil.copy2(DATA_DIR / fname, frame_dir / fname)
        result = trajectory_surface_charge(tmp_path, dir_pattern="bader_t*_i*")
        assert result.shape == (4, 2)
        # All frames are identical data, so values should match
        np.testing.assert_allclose(result[0], result[3], atol=1e-10)


class TestSurfaceChargeAnalysis:
    """Integration tests for surface_charge_analysis end-to-end."""

    def test_end_to_end(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=3)
        out = tmp_path / "output"
        csv_path = surface_charge_analysis(root, output_dir=out)

        assert csv_path.exists()
        assert (out / "surface_charge.png").exists()

        with csv_path.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 3

        # All frames identical → cumulative average converges to inst. value
        bot_0 = float(rows[0]["sigma_aligned_uC_cm2"])
        bot_cum_2 = float(rows[2]["sigma_aligned_cumavg_uC_cm2"])
        assert bot_cum_2 == pytest.approx(bot_0, rel=1e-10)

    def test_numeric_sort_in_analysis(self, tmp_path):
        """surface_charge_analysis uses numeric sort for non-zero-padded dirs."""
        for t in [1000, 50]:
            frame_dir = tmp_path / f"bader_t{t}_i0"
            frame_dir.mkdir()
            for fname in _FRAME_FILES:
                shutil.copy2(DATA_DIR / fname, frame_dir / fname)
        out = tmp_path / "output"
        csv_path = surface_charge_analysis(
            tmp_path, output_dir=out, dir_pattern="bader_t*_i*",
        )
        with csv_path.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        # First row should be t=50, second t=1000
        assert int(rows[0]["step"]) == 50
        assert int(rows[1]["step"]) == 1000
