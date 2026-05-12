"""Integration tests for the ``md_analysis.main`` re-export facade.

Phase 7a of the entrance refactor removed the legacy
``run_*_analysis`` and ``run_all`` shims from ``md_analysis.main``.
This module now re-exports the canonical workflow functions from
``md_analysis.workflows``; the tests below cover the same
scientific paths as before but against the new ``WorkflowResult``
return contract.

Fixtures live under ``data_example/potential/dense/`` (continuous
mode) — the legacy ``data_example/potential/`` root-level files were
moved into ``dense/`` during the Phase 0 utils/engines refactor.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import pytest

from md_analysis.main import (
    WorkflowResult,
    run_interface_analysis,
    run_potential_full,
    run_surface_charge,
    run_water_three_panel,
)

# ---------------------------------------------------------------------------
# Data directories
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATA_DIR = _REPO_ROOT / "data_example" / "potential" / "dense"
_BADER_DIR = _REPO_ROOT / "data_example" / "bader" / "single_frame"

pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"data_example/potential/dense/ not found at {_DATA_DIR}",
)

# Files needed per fake charge frame
_FRAME_FILES = ["POSCAR", "ACF.dat", "POTCAR"]


def _build_fake_trajectory(tmp_path: Path, n_frames: int = 2) -> Path:
    """Copy single_frame data into bader_t*_i* subdirectories."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    for i in range(n_frames):
        frame_dir = tmp_path / f"bader_t{i:03d}_i000"
        frame_dir.mkdir()
        for fname in _FRAME_FILES:
            shutil.copy2(_BADER_DIR / fname, frame_dir / fname)
    return tmp_path


# ===========================================================================
# run_water_three_panel
# ===========================================================================


class TestRunWaterThreePanel:

    def test_happy_path(self, tmp_path: Path):
        water_dir = tmp_path / "water"
        result = run_water_three_panel(
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
            md_inp_path=_DATA_DIR / "md.inp",
            output_dir=water_dir,
        )

        assert isinstance(result, WorkflowResult)
        assert result.name == "water_three_panel"

        expected_keys = {
            "density_csv",
            "orientation_csv",
            "adsorbed_profile_csv",
            "adsorbed_range_txt",
            "adsorbed_theta_csv",
            "plot_png",
        }
        assert set(result.artifacts.keys()) == expected_keys

        # All files exist and water_dir was created
        assert water_dir.is_dir()
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_frame_slicing(self, tmp_path: Path):
        result = run_water_three_panel(
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
            md_inp_path=_DATA_DIR / "md.inp",
            output_dir=tmp_path,
            frame_start=0,
            frame_end=5,
            frame_step=2,
        )

        assert isinstance(result, WorkflowResult)
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"


# ===========================================================================
# run_potential_full
# ===========================================================================


class TestRunPotentialFull:

    def test_full_electrode(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            result = run_potential_full(
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

        assert isinstance(result, WorkflowResult)
        assert "electrode_csv" in result.artifacts
        assert "phi_z_png" in result.artifacts
        assert "thickness_sensitivity_csv" in result.artifacts
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_separate_center_fermi(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            result = run_potential_full(
                output_dir=tmp_path,
                md_out_path=_DATA_DIR / "md.out",
                compute_u=False,
                compute_phi_z=False,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "center_csv" in result.artifacts
        assert "fermi_csv" in result.artifacts
        assert "electrode_csv" not in result.artifacts
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_no_md_out(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            result = run_potential_full(
                output_dir=tmp_path,
                md_out_path=None,
                compute_phi_z=True,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "center_csv" in result.artifacts
        assert "phi_z_png" in result.artifacts
        assert "fermi_csv" not in result.artifacts
        assert "electrode_csv" not in result.artifacts
        assert "thickness_sensitivity_csv" not in result.artifacts
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"

    def test_no_phi_z(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            result = run_potential_full(
                output_dir=tmp_path,
                md_out_path=None,
                compute_phi_z=False,
                thickness_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert "phi_z_png" not in result.artifacts


# ===========================================================================
# run_surface_charge
# ===========================================================================


class TestRunSurfaceCharge:

    @pytest.mark.skipif(
        not _BADER_DIR.exists(),
        reason=f"data_example/bader/single_frame/ not found at {_BADER_DIR}",
    )
    def test_happy_path(self, tmp_path: Path):
        root = _build_fake_trajectory(tmp_path / "traj", n_frames=2)
        # run_surface_charge writes into output_dir/<method>/.
        out = tmp_path / "output" / "electrochemical" / "charge"

        result = run_surface_charge(
            output_dir=out,
            root_dir=root,
            method="counterion",
        )

        assert isinstance(result, WorkflowResult)
        assert "charge_csv" in result.artifacts
        assert "charge_png" in result.artifacts
        assert result.artifacts["charge_csv"].exists()
        # method subdirectory was created inside output_dir
        assert (out / "counterion").is_dir()


# ===========================================================================
# run_interface_analysis (replaces legacy run_all)
# ===========================================================================


class TestRunInterfaceAnalysis:

    def test_happy_path(self, tmp_path: Path):
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            result = run_interface_analysis(
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

        assert isinstance(result, WorkflowResult)
        assert result.name == "interface_analysis"
        # Should contain both water and potential keys.
        assert "density_csv" in result.artifacts
        assert "electrode_csv" in result.artifacts
        for key, path in result.artifacts.items():
            assert path.exists(), f"{key} not found: {path}"
        # Composite metadata
        assert result.metadata["sub_workflows"] == [
            "water_three_panel", "potential_full",
        ]
        assert result.metadata["water_output_dir"].endswith("water")
        assert result.metadata["potential_output_dir"].endswith(
            "electrochemical/potential"
        )
