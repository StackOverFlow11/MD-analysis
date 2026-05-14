"""Unit tests for ``md_analysis.workflows.potential`` single-step facades.

Phase 6.2 adds 5 single-step run_* functions that wrap the existing
electrochemical.potential business functions:

  - run_center_potential        -> center_csv
  - run_fermi_energy            -> fermi_csv
  - run_electrode_potential     -> electrode_csv
  - run_phi_z_profile           -> phi_z_png
  - run_thickness_sensitivity   -> thickness_sensitivity_csv

These tests monkeypatch the underlying business functions (mock,
no I/O) and pin the facade contract: WorkflowResult shape, artifact
keys, output_dir mkdir, metadata fields, and key kwargs forwarding
in both continuous and distributed input modes.

The composite ``run_potential_full`` is exercised end-to-end by the
existing integration tests under ``test/integration/potential`` and is
not re-covered here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from md_analysis.workflows import (
    WorkflowResult,
    run_center_potential,
    run_electrode_potential,
    run_fermi_energy,
    run_phi_z_profile,
    run_thickness_sensitivity,
)


# ---------------------------------------------------------------------------
# run_center_potential
# ---------------------------------------------------------------------------


class TestRunCenterPotential:
    def test_continuous_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(cube_pattern: str, **kwargs: Any) -> Path:
            captured["cube_pattern"] = cube_pattern
            captured.update(kwargs)
            return tmp_path / "center" / "center.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.center_slab_potential_analysis",
            fake_business,
        )

        result = run_center_potential(
            output_dir=tmp_path / "center",
            cube_pattern="md-POTENTIAL-v_hartree-1_*.cube",
            xyz_path=tmp_path / "md-pos-1.xyz",
            thickness_ang=7.0,
            center_mode="interface",
            metal_elements={"Cu"},
            layer_tol_ang=0.5,
            frame_start=0,
            frame_end=10,
            frame_step=1,
            verbose=True,
        )

        # WorkflowResult shape
        assert isinstance(result, WorkflowResult)
        assert result.name == "center_potential"
        assert result.output_dir == tmp_path / "center"
        assert (tmp_path / "center").exists()
        assert set(result.artifacts) == {"center_csv"}
        assert result.artifacts["center_csv"] == tmp_path / "center" / "center.csv"

        # metadata contract
        assert result.metadata["input_mode"] == "continuous"
        assert result.metadata["center_mode"] == "interface"
        assert result.metadata["thickness_ang"] == 7.0

        # business call kwargs forwarded
        assert captured["cube_pattern"] == "md-POTENTIAL-v_hartree-1_*.cube"
        assert captured["thickness_ang"] == 7.0
        assert captured["xyz_path"] == tmp_path / "md-pos-1.xyz"
        assert captured["input_mode"] == "continuous"

    def test_distributed_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(cube_pattern: str, **kwargs: Any) -> Path:
            captured["cube_pattern"] = cube_pattern
            captured.update(kwargs)
            return tmp_path / "center" / "center.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.center_slab_potential_analysis",
            fake_business,
        )

        result = run_center_potential(
            output_dir=tmp_path / "center",
            input_mode="distributed",
            sp_root_dir=tmp_path / "sp_dirs",
            sp_dir_pattern="potential_t*_i*",
        )

        assert result.name == "center_potential"
        assert result.metadata["input_mode"] == "distributed"
        # business receives input_mode="distributed" + sp_root_dir Path
        assert captured["input_mode"] == "distributed"
        assert captured["sp_root_dir"] == tmp_path / "sp_dirs"


# ---------------------------------------------------------------------------
# run_fermi_energy
# ---------------------------------------------------------------------------


class TestRunFermiEnergy:
    def test_continuous_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(md_out: Path | None, **kwargs: Any) -> Path:
            captured["md_out"] = md_out
            captured.update(kwargs)
            return tmp_path / "fermi" / "fermi.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.fermi_energy_analysis",
            fake_business,
        )

        result = run_fermi_energy(
            output_dir=tmp_path / "fermi",
            md_out_path=tmp_path / "md.out",
            fermi_unit="au",
        )

        assert result.name == "fermi_energy"
        assert set(result.artifacts) == {"fermi_csv"}
        assert result.metadata["fermi_unit"] == "au"
        assert captured["md_out"] == tmp_path / "md.out"
        assert captured["input_mode"] == "continuous"

    def test_distributed_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(md_out: Path | None, **kwargs: Any) -> Path:
            captured["md_out"] = md_out
            captured.update(kwargs)
            return tmp_path / "fermi" / "fermi.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.fermi_energy_analysis",
            fake_business,
        )

        result = run_fermi_energy(
            output_dir=tmp_path / "fermi",
            input_mode="distributed",
            sp_root_dir=tmp_path / "sp_dirs",
        )

        assert result.metadata["input_mode"] == "distributed"
        assert captured["md_out"] is None
        assert captured["input_mode"] == "distributed"


# ---------------------------------------------------------------------------
# run_electrode_potential
# ---------------------------------------------------------------------------


class TestRunElectrodePotential:
    def test_continuous_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(
            cube_pattern: str, md_out: Path | None, **kwargs: Any
        ) -> Path:
            captured["cube_pattern"] = cube_pattern
            captured["md_out"] = md_out
            captured.update(kwargs)
            return tmp_path / "electrode" / "electrode.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.electrode_potential_analysis",
            fake_business,
        )

        result = run_electrode_potential(
            output_dir=tmp_path / "electrode",
            cube_pattern="md-POTENTIAL-*.cube",
            md_out_path=tmp_path / "md.out",
            xyz_path=tmp_path / "md-pos-1.xyz",
        )

        assert result.name == "electrode_potential"
        assert set(result.artifacts) == {"electrode_csv"}
        assert captured["cube_pattern"] == "md-POTENTIAL-*.cube"
        assert captured["md_out"] == tmp_path / "md.out"


# ---------------------------------------------------------------------------
# run_phi_z_profile
# ---------------------------------------------------------------------------


class TestRunPhiZProfile:
    def test_continuous_mode_png_artifact(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(cube_pattern: str, **kwargs: Any) -> Path:
            captured["cube_pattern"] = cube_pattern
            captured.update(kwargs)
            return tmp_path / "phi_z" / "phi_z.png"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.phi_z_planeavg_analysis",
            fake_business,
        )

        result = run_phi_z_profile(
            output_dir=tmp_path / "phi_z",
            cube_pattern="md-POTENTIAL-*.cube",
            max_curves=5,
        )

        # PhiZ is the only single-step potential facade returning a PNG
        assert result.name == "phi_z_profile"
        assert set(result.artifacts) == {"phi_z_png"}
        assert result.artifacts["phi_z_png"].suffix == ".png"
        assert result.metadata["max_curves"] == 5


# ---------------------------------------------------------------------------
# run_thickness_sensitivity
# ---------------------------------------------------------------------------


class TestRunThicknessSensitivity:
    def test_distributed_mode_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(
            cube_pattern: str, md_out: Path | None, **kwargs: Any
        ) -> Path:
            captured["cube_pattern"] = cube_pattern
            captured["md_out"] = md_out
            captured.update(kwargs)
            return tmp_path / "ts" / "thickness_sensitivity.csv"

        monkeypatch.setattr(
            "md_analysis.electrochemical.potential.thickness_sensitivity_analysis",
            fake_business,
        )

        result = run_thickness_sensitivity(
            output_dir=tmp_path / "ts",
            thickness_end=12.0,
            input_mode="distributed",
            sp_root_dir=tmp_path / "sp_dirs",
        )

        assert result.name == "thickness_sensitivity"
        assert set(result.artifacts) == {"thickness_sensitivity_csv"}
        assert result.metadata["thickness_end"] == 12.0
        assert result.metadata["input_mode"] == "distributed"
        assert captured["input_mode"] == "distributed"
        assert captured["sp_root_dir"] == tmp_path / "sp_dirs"
