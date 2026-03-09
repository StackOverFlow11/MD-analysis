"""Integration tests for md_analysis.main programmatic API."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import pytest

from md_analysis.main import (
    run_all,
    run_charge_analysis,
    run_potential_analysis,
    run_water_analysis,
)

# ---------------------------------------------------------------------------
# Data directories
# ---------------------------------------------------------------------------

_DATA_DIR = Path(__file__).resolve().parents[2] / "data_example" / "potential"
_BADER_DIR = Path(__file__).resolve().parents[2] / "data_example" / "bader_work_dir"

pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"data_example/potential/ not found at {_DATA_DIR}",
)

# Files needed per fake charge frame
_FRAME_FILES = ["POSCAR", "ACF.dat", "POTCAR"]


def _build_fake_trajectory(tmp_path: Path, n_frames: int = 2) -> Path:
    """Copy bader_work_dir data into bader_t*_i* subdirectories."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    for i in range(n_frames):
        frame_dir = tmp_path / f"bader_t{i:03d}_i000"
        frame_dir.mkdir()
        for fname in _FRAME_FILES:
            shutil.copy2(_BADER_DIR / fname, frame_dir / fname)
    return tmp_path


# ===========================================================================
# TestRunWaterAnalysis
# ===========================================================================


class TestRunWaterAnalysis:

    def test_happy_path(self, tmp_path: Path):
        results = run_water_analysis(
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
            md_inp_path=_DATA_DIR / "md.inp",
            output_dir=tmp_path,
        )

        assert len(results) == 6
        expected_keys = {
            "density_csv",
            "orientation_csv",
            "adsorbed_profile_csv",
            "adsorbed_range_txt",
            "adsorbed_theta_csv",
            "plot_png",
        }
        assert set(results.keys()) == expected_keys

        # All files exist and water/ subdirectory was created
        assert (tmp_path / "water").is_dir()
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_frame_slicing(self, tmp_path: Path):
        results = run_water_analysis(
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
            md_inp_path=_DATA_DIR / "md.inp",
            output_dir=tmp_path,
            frame_start=0,
            frame_end=5,
            frame_step=2,
        )

        # Should succeed without error; outputs exist
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"


# ===========================================================================
# TestRunPotentialAnalysis
# ===========================================================================


class TestRunPotentialAnalysis:

    def test_full_electrode(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            results = run_potential_analysis(
                output_dir=tmp_path,
                md_out_path=_DATA_DIR / "md.out",
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
                compute_u=True,
                compute_phi_z=True,
                thickness_ang=7.0,
                center_mode="interface",
            )
        finally:
            os.chdir(old_cwd)

        assert "electrode_csv" in results
        assert "phi_z_png" in results
        assert "thickness_sensitivity_csv" in results
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_separate_center_fermi(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            results = run_potential_analysis(
                output_dir=tmp_path,
                md_out_path=_DATA_DIR / "md.out",
                compute_u=False,
                compute_phi_z=False,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "center_csv" in results
        assert "fermi_csv" in results
        assert "electrode_csv" not in results
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_no_md_out(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            results = run_potential_analysis(
                output_dir=tmp_path,
                md_out_path=None,
                compute_phi_z=True,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "center_csv" in results
        assert "phi_z_png" in results
        assert "fermi_csv" not in results
        assert "electrode_csv" not in results
        assert "thickness_sensitivity_csv" not in results
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_no_phi_z(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            results = run_potential_analysis(
                output_dir=tmp_path,
                md_out_path=None,
                compute_phi_z=False,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "phi_z_png" not in results


# ===========================================================================
# TestRunChargeAnalysis
# ===========================================================================


class TestRunChargeAnalysis:

    @pytest.mark.skipif(
        not _BADER_DIR.exists(),
        reason=f"data_example/bader_work_dir/ not found at {_BADER_DIR}",
    )
    def test_happy_path(self, tmp_path: Path):
        root = _build_fake_trajectory(tmp_path / "traj", n_frames=2)
        out = tmp_path / "output"

        results = run_charge_analysis(
            output_dir=out,
            root_dir=root,
            method="counterion",
        )

        assert "charge_csv" in results
        assert "charge_png" in results
        assert results["charge_csv"].exists()
        # charge/counterion/ subdirectory was created
        assert (out / "charge" / "counterion").is_dir()


# ===========================================================================
# TestRunAll
# ===========================================================================


class TestRunAll:

    def test_happy_path(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            results = run_all(
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
                md_inp_path=_DATA_DIR / "md.inp",
                output_dir=tmp_path,
                md_out_path=_DATA_DIR / "md.out",
                compute_u=True,
                compute_phi_z=True,
                thickness_ang=7.0,
                center_mode="interface",
            )
        finally:
            os.chdir(old_cwd)

        # Should contain both water and potential keys
        assert "density_csv" in results
        assert "electrode_csv" in results
        for key, path in results.items():
            assert path.exists(), f"{key} not found: {path}"
