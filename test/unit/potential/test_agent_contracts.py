"""Batch 4: contract-backed ``potential_full``.

Covers the ``run_potential_analysis_with_report`` wrapper plus the
agent dispatch path.  Continuous-mode integration tests have to
``chdir()`` into the cube directory because ``run_potential_analysis``
resolves ``cube_pattern`` relative to the process current working
directory — the wrapper preserves that historical semantic.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.main import (
    PotentialFullReport,
    run_potential_analysis_with_report,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "data_example" / "potential" / "dense"
MD_OUT = DATA_DIR / "md.out"
XYZ_PATH = DATA_DIR / "md-pos-1.xyz"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(),
    reason=f"data_example/potential/dense/ not found at {DATA_DIR}",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _Chdir:
    """Context manager around ``os.chdir`` — the wrapper preserves the
    legacy continuous-mode behaviour where ``cube_pattern`` is resolved
    relative to the process working directory."""

    def __init__(self, target: Path):
        self._target = target
        self._saved: str | None = None

    def __enter__(self):
        self._saved = os.getcwd()
        os.chdir(self._target)
        return self

    def __exit__(self, *_exc):
        assert self._saved is not None
        os.chdir(self._saved)


# ---------------------------------------------------------------------------
# Wrapper-level
# ---------------------------------------------------------------------------


class TestRunPotentialAnalysisWithReport:
    def test_continuous_success_has_fermi(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_potential_analysis_with_report(
                output_dir=tmp_path,
                md_out_path=MD_OUT,
                xyz_path=XYZ_PATH,
            )
        assert isinstance(r, PotentialFullReport)
        assert r.input_mode == "continuous"
        # Every advertised artifact must exist on disk
        for key, path in r.artifacts.items():
            assert path.is_file(), f"missing artifact {key}: {path}"
        # Fermi-dependent artifacts present → has_fermi is True
        assert r.has_fermi is True
        assert r.ran_electrode is True
        assert r.ran_thickness_sensitivity is True
        assert r.n_artifacts == len(r.artifacts)

    def test_compute_phi_z_false_omits_phi_z_png(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_potential_analysis_with_report(
                output_dir=tmp_path,
                md_out_path=MD_OUT,
                xyz_path=XYZ_PATH,
                compute_phi_z=False,
            )
        assert "phi_z_png" not in r.artifacts
        assert r.ran_phi_z is False

    def test_no_md_out_omits_fermi_dependent_artifacts(self, tmp_path):
        """When md_out_path is omitted, the legacy function skips all
        Fermi-dependent sub-analyses.  The wrapper's ``has_fermi`` must
        reflect that."""
        with _Chdir(DATA_DIR):
            r = run_potential_analysis_with_report(
                output_dir=tmp_path,
                md_out_path=None,
                xyz_path=XYZ_PATH,
            )
        assert r.has_fermi is False
        assert r.ran_electrode is False
        assert r.ran_fermi is False
        assert r.ran_thickness_sensitivity is False
        # center_csv should still have been produced
        assert r.ran_center is True

    def test_report_to_dict_is_json_serializable(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_potential_analysis_with_report(
                output_dir=tmp_path,
                md_out_path=MD_OUT,
                xyz_path=XYZ_PATH,
            )
        json.dumps(r.to_dict())

    def test_invalid_input_mode_is_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="input_mode"):
            run_potential_analysis_with_report(
                output_dir=tmp_path,
                input_mode="bogus",
            )

    def test_invalid_center_mode_is_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="center_mode"):
            run_potential_analysis_with_report(
                output_dir=tmp_path,
                center_mode="bogus",
            )

    def test_invalid_fermi_unit_is_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="fermi_unit"):
            run_potential_analysis_with_report(
                output_dir=tmp_path,
                fermi_unit="bogus",
            )


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


class TestPotentialFullSchema:
    def test_schema_keys(self):
        s = get_task_schema("potential_full")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("potential_full")
        assert s["parameters"]["required"] == ["output_dir"]

    def test_input_mode_enum(self):
        s = get_task_schema("potential_full")
        assert set(
            s["parameters"]["properties"]["input_mode"].get("enum", [])
        ) == {"continuous", "distributed"}

    def test_center_mode_enum(self):
        s = get_task_schema("potential_full")
        assert set(
            s["parameters"]["properties"]["center_mode"].get("enum", [])
        ) == {"interface", "cell"}

    def test_fermi_unit_enum(self):
        s = get_task_schema("potential_full")
        assert set(
            s["parameters"]["properties"]["fermi_unit"].get("enum", [])
        ) == {"au", "ev"}

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["potential_full"].contract is not None


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


class TestPotentialFullDispatch:
    def test_continuous_success(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = dispatch("potential_full", {
                "output_dir": str(tmp_path),
                "md_out_path": str(MD_OUT),
                "xyz_path": str(XYZ_PATH),
            })
        assert r.success, r.errors
        # All advertised outputs exist on disk
        for k, path in r.outputs.items():
            assert Path(path).is_file(), (
                f"outputs[{k!r}] points to a non-existent file: {path}"
            )
        assert r.summary["input_mode"] == "continuous"
        assert r.summary["has_fermi"] is True
        # Summary booleans match outputs keys presence
        assert r.summary["ran_electrode"] == ("electrode_csv" in r.outputs)
        assert r.summary["ran_center"] == ("center_csv" in r.outputs)
        assert r.summary["ran_fermi"] == ("fermi_csv" in r.outputs)
        assert r.summary["ran_phi_z"] == ("phi_z_png" in r.outputs)
        assert r.summary["ran_thickness_sensitivity"] == (
            "thickness_sensitivity_csv" in r.outputs
        )
        json.dumps(r.summary)

    def test_invalid_input_mode_is_validation(self, tmp_path):
        r = dispatch("potential_full", {
            "output_dir": str(tmp_path),
            "input_mode": "bogus",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_invalid_center_mode_is_validation(self, tmp_path):
        r = dispatch("potential_full", {
            "output_dir": str(tmp_path),
            "center_mode": "bogus",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_invalid_fermi_unit_is_validation(self, tmp_path):
        r = dispatch("potential_full", {
            "output_dir": str(tmp_path),
            "fermi_unit": "bogus",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_continuous_no_matching_cube_is_file_not_found(self, tmp_path):
        """Run continuous mode from an empty dir — no cube files match
        the default pattern so the lower-level analysis raises
        FileNotFoundError."""
        empty = tmp_path / "empty"
        empty.mkdir()
        with _Chdir(empty):
            r = dispatch("potential_full", {
                "output_dir": str(tmp_path / "out"),
                "md_out_path": str(MD_OUT),
                "xyz_path": str(XYZ_PATH),
            })
        assert not r.success
        assert r.error_type == "file_not_found"
