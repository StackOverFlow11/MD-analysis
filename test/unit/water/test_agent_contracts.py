"""Batch 4: contract-backed ``water_three_panel``.

Covers the ``run_water_analysis_with_report`` wrapper and its
end-to-end agent dispatch.  Uses ``data_example/potential/`` (which
already serves the integration tests) as the real MD fixture so we
exercise the full downstream pipeline instead of stubbing it out.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.main import (
    WaterThreePanelReport,
    run_water_analysis_with_report,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "data_example" / "potential"
XYZ_PATH = DATA_DIR / "md-pos-1.xyz"
MD_INP_PATH = DATA_DIR / "md.inp"

pytestmark = pytest.mark.skipif(
    not (XYZ_PATH.is_file() and MD_INP_PATH.is_file()),
    reason=(
        f"fixture inputs missing: {XYZ_PATH} or {MD_INP_PATH} "
        "(expected under data_example/potential/)"
    ),
)


# ---------------------------------------------------------------------------
# Wrapper-level
# ---------------------------------------------------------------------------


class TestRunWaterAnalysisWithReport:
    def test_success_all_six_artifacts(self, tmp_path):
        r = run_water_analysis_with_report(
            xyz_path=XYZ_PATH,
            md_inp_path=MD_INP_PATH,
            output_dir=tmp_path,
        )
        assert isinstance(r, WaterThreePanelReport)
        expected_keys = {
            "density_csv",
            "orientation_csv",
            "adsorbed_profile_csv",
            "adsorbed_range_txt",
            "adsorbed_theta_csv",
            "plot_png",
        }
        assert set(r.artifacts.keys()) == expected_keys
        # Every advertised artifact must exist on disk
        for key, path in r.artifacts.items():
            assert path.is_file(), f"missing artifact {key}: {path}"
        assert r.n_artifacts == 6

    def test_report_to_dict_is_json_serializable(self, tmp_path):
        r = run_water_analysis_with_report(
            xyz_path=XYZ_PATH,
            md_inp_path=MD_INP_PATH,
            output_dir=tmp_path,
        )
        json.dumps(r.to_dict())

    def test_slice_echoes_into_summary(self, tmp_path):
        r = run_water_analysis_with_report(
            xyz_path=XYZ_PATH,
            md_inp_path=MD_INP_PATH,
            output_dir=tmp_path,
            frame_start=0,
            frame_end=5,
            frame_step=1,
        )
        assert r.frame_start == 0
        assert r.frame_end == 5
        assert r.frame_step == 1


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


class TestWaterThreePanelSchema:
    def test_schema_keys(self):
        s = get_task_schema("water_three_panel")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("water_three_panel")
        assert set(s["parameters"]["required"]) == {"xyz_path", "output_dir"}

    def test_cell_abc_is_nullable_length_3(self):
        s = get_task_schema("water_three_panel")
        ca = s["parameters"]["properties"]["cell_abc"]
        assert ca["type"] == ["array", "null"]
        assert ca["minItems"] == 3 and ca["maxItems"] == 3
        assert ca["items"]["type"] == "number"

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["water_three_panel"].contract is not None


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


class TestWaterThreePanelDispatch:
    def test_success(self, tmp_path):
        r = dispatch("water_three_panel", {
            "xyz_path": str(XYZ_PATH),
            "md_inp_path": str(MD_INP_PATH),
            "output_dir": str(tmp_path),
        })
        assert r.success, r.errors
        for key in (
            "density_csv",
            "orientation_csv",
            "adsorbed_profile_csv",
            "adsorbed_range_txt",
            "adsorbed_theta_csv",
            "plot_png",
        ):
            assert key in r.outputs
            assert Path(r.outputs[key]).is_file()
        assert r.summary["n_artifacts"] == 6
        json.dumps(r.summary)

    def test_missing_xyz_is_file_not_found(self, tmp_path):
        r = dispatch("water_three_panel", {
            "xyz_path": str(tmp_path / "nope.xyz"),
            "md_inp_path": str(MD_INP_PATH),
            "output_dir": str(tmp_path),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_cell_abc_length_is_validation(self, tmp_path):
        r = dispatch("water_three_panel", {
            "xyz_path": str(XYZ_PATH),
            "cell_abc": [10.0, 10.0],  # too short
            "output_dir": str(tmp_path),
        })
        assert not r.success
        assert r.error_type == "validation"
