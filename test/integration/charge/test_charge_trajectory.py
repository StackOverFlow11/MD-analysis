"""Integration tests for trajectory_charge_analysis."""

import shutil
from pathlib import Path

import numpy as np
import pytest

from md_analysis.charge.ChargeAnalysis import (
    ElementSelector,
    trajectory_charge_analysis,
)
from md_analysis.charge.config import (
    DEFAULT_SELECTED_ATOM_CHARGES_CSV_NAME,
    DEFAULT_SURFACE_CHARGE_CSV_NAME,
)

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


class TestTrajectoryChargeAnalysis:
    """Integration tests using real Bader data arranged as a fake trajectory."""

    def test_basic_trajectory(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        result = trajectory_charge_analysis(root)

        assert len(result.frame_labels) == 2
        assert result.surface_charge_density_uC_cm2.shape == (2, 2)
        assert result.mean_surface_charge_density_uC_cm2.shape == (2,)
        assert result.std_surface_charge_density_uC_cm2.shape == (2,)
        assert result.selected_atom_indices == ()
        assert result.selected_atom_net_charges.shape == (2, 0)

    def test_labels_sorted(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=3)
        result = trajectory_charge_analysis(root)
        assert result.frame_labels == (
            "calc_t000_i000",
            "calc_t001_i000",
            "calc_t002_i000",
        )

    def test_csv_output(self, tmp_path):
        root = _build_fake_trajectory(tmp_path)
        trajectory_charge_analysis(root)

        csv_path = root / DEFAULT_SURFACE_CHARGE_CSV_NAME
        assert csv_path.exists()
        lines = csv_path.read_text().splitlines()
        # header + 2 frames + mean + std = 5 lines
        assert len(lines) == 5
        assert lines[0].startswith("frame")

    def test_with_element_selector(self, tmp_path):
        root = _build_fake_trajectory(tmp_path)
        sel = ElementSelector({"Ag"})
        result = trajectory_charge_analysis(root, atom_selector=sel)

        assert len(result.selected_atom_indices) == 2  # 2 Ag atoms
        assert result.selected_atom_net_charges.shape == (2, 2)
        assert result.mean_selected_atom_net_charges.shape == (2,)

        csv_path = root / DEFAULT_SELECTED_ATOM_CHARGES_CSV_NAME
        assert csv_path.exists()

    def test_identical_frames_zero_std(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=3)
        result = trajectory_charge_analysis(root)

        # All frames are identical copies, so std should be ~0
        np.testing.assert_allclose(
            result.std_surface_charge_density_uC_cm2, 0.0, atol=1e-10
        )

    def test_custom_output_dir(self, tmp_path):
        root = _build_fake_trajectory(tmp_path)
        out = tmp_path / "output"
        trajectory_charge_analysis(root, output_dir=out)

        assert (out / DEFAULT_SURFACE_CHARGE_CSV_NAME).exists()

    def test_missing_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No subdirectories"):
            trajectory_charge_analysis(tmp_path)

    def test_skip_frame_without_acf(self, tmp_path):
        root = _build_fake_trajectory(tmp_path, n_frames=2)
        # Remove ACF.dat from one frame
        (root / "calc_t001_i000" / "ACF.dat").unlink()

        with pytest.warns(UserWarning, match="Skipping"):
            result = trajectory_charge_analysis(root)

        assert len(result.frame_labels) == 1
