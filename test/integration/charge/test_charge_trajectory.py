"""Integration tests for trajectory_indexed_atom_charges."""

import shutil
from pathlib import Path

import numpy as np
import pytest

from md_analysis.charge.ChargeAnalysis import trajectory_indexed_atom_charges

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "bader_work_dir"

# Files needed per frame
_FRAME_FILES = ["POSCAR", "ACF.dat", "POTCAR"]


def _build_fake_trajectory(tmp_path: Path, n_frames: int = 2) -> Path:
    """Copy bader_work_dir data into calc_t*_i* subdirectories."""
    for i in range(n_frames):
        frame_dir = tmp_path / f"calc_t{i:03d}_i000"
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
