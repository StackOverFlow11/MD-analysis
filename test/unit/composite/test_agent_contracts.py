"""Batch 5: contract-backed ``run_all``.

``run_all`` is a convenience composite of water + potential.  These
tests reuse the Batch 4 leaf test fixture (``data_example/potential/``)
and ``chdir`` into it so continuous-mode potential cube discovery
succeeds — same pattern as the potential agent tests.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.main import RunAllReport, run_all_with_report

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "data_example" / "potential"
XYZ_PATH = DATA_DIR / "md-pos-1.xyz"
MD_INP_PATH = DATA_DIR / "md.inp"
MD_OUT = DATA_DIR / "md.out"

pytestmark = pytest.mark.skipif(
    not (
        XYZ_PATH.is_file()
        and MD_INP_PATH.is_file()
        and MD_OUT.is_file()
    ),
    reason=f"fixture inputs missing under {DATA_DIR}",
)


class _Chdir:
    """Context manager around ``os.chdir`` — continuous-mode potential
    resolves ``cube_pattern`` relative to the process cwd."""

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


_WATER_KEYS = {
    "density_csv",
    "orientation_csv",
    "adsorbed_profile_csv",
    "adsorbed_range_txt",
    "adsorbed_theta_csv",
    "plot_png",
}


# ---------------------------------------------------------------------------
# Wrapper-level
# ---------------------------------------------------------------------------


class TestRunAllWithReport:
    def test_success_returns_report(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_all_with_report(
                xyz_path=XYZ_PATH,
                md_inp_path=MD_INP_PATH,
                output_dir=tmp_path,
                md_out_path=MD_OUT,
            )
        assert isinstance(r, RunAllReport)
        # All artifact paths exist on disk
        for key, path in r.artifacts.items():
            assert path.is_file(), f"missing artifact {key}: {path}"
        # Water leaf always produces 6 artifacts on success
        assert r.water_n_artifacts == 6
        assert _WATER_KEYS <= set(r.artifacts.keys())
        # Totals align
        assert r.n_artifacts == r.water_n_artifacts + r.potential_n_artifacts
        # Flags mirror produced artifacts, not input flags
        assert r.ran_water is True
        assert r.ran_potential is True
        assert r.potential_has_fermi is True  # md.out supplied
        assert r.potential_ran_electrode is True
        assert r.potential_ran_thickness_sensitivity is True

    def test_preserves_output_layout(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_all_with_report(
                xyz_path=XYZ_PATH,
                md_inp_path=MD_INP_PATH,
                output_dir=tmp_path,
                md_out_path=MD_OUT,
            )
        assert r.water_output_dir == tmp_path / "water"
        assert (
            r.potential_output_dir
            == tmp_path / "electrochemical" / "potential"
        )
        # Artifact paths reside under the correct sub-dir
        for key in _WATER_KEYS:
            assert tmp_path / "water" in r.artifacts[key].parents

    def test_report_to_dict_is_json_serializable(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = run_all_with_report(
                xyz_path=XYZ_PATH,
                md_inp_path=MD_INP_PATH,
                output_dir=tmp_path,
                md_out_path=MD_OUT,
            )
        json.dumps(r.to_dict())

    def test_no_md_out_skips_fermi_dependent_potential(self, tmp_path):
        """Without ``md_out_path``, potential produces only center_csv +
        phi_z_png (Fermi-dependent artifacts are omitted)."""
        with _Chdir(DATA_DIR):
            r = run_all_with_report(
                xyz_path=XYZ_PATH,
                md_inp_path=MD_INP_PATH,
                output_dir=tmp_path,
                md_out_path=None,
            )
        assert r.ran_water is True
        assert r.ran_potential is True
        assert r.potential_has_fermi is False
        assert r.potential_ran_electrode is False
        assert r.potential_ran_fermi is False
        assert r.potential_ran_thickness_sensitivity is False
        # center_csv survives; Fermi-dependent keys absent
        assert "center_csv" in r.artifacts
        assert "electrode_csv" not in r.artifacts
        assert "fermi_csv" not in r.artifacts
        assert "thickness_sensitivity_csv" not in r.artifacts

    def test_unknown_potential_kwarg_is_ignored(self, tmp_path):
        """Matches legacy ``run_all`` behaviour: non-whitelisted kwargs
        are silently dropped when forwarded to potential."""
        with _Chdir(DATA_DIR):
            r = run_all_with_report(
                xyz_path=XYZ_PATH,
                md_inp_path=MD_INP_PATH,
                output_dir=tmp_path,
                md_out_path=MD_OUT,
                # silently ignored
                bogus_key="nope",
            )
        assert isinstance(r, RunAllReport)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


class TestRunAllSchema:
    def test_schema_keys(self):
        s = get_task_schema("run_all")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("run_all")
        assert set(s["parameters"]["required"]) == {"xyz_path", "output_dir"}

    def test_exact_public_input_fields(self):
        """Conservative Batch 5 surface: only legacy signature fields."""
        s = get_task_schema("run_all")
        assert set(s["parameters"]["properties"].keys()) == {
            "xyz_path", "md_inp_path", "cell_abc", "output_dir",
            "cube_pattern", "md_out_path",
            "frame_start", "frame_end", "frame_step",
            "verbose",
        }

    def test_cell_abc_is_nullable_length_3(self):
        s = get_task_schema("run_all")
        ca = s["parameters"]["properties"]["cell_abc"]
        assert ca["type"] == ["array", "null"]
        assert ca["minItems"] == 3 and ca["maxItems"] == 3

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["run_all"].contract is not None


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


class TestRunAllDispatch:
    def test_success(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = dispatch("run_all", {
                "xyz_path": str(XYZ_PATH),
                "md_inp_path": str(MD_INP_PATH),
                "output_dir": str(tmp_path),
                "md_out_path": str(MD_OUT),
            })
        assert r.success, r.errors
        # Every advertised output actually exists on disk
        for k, path in r.outputs.items():
            assert Path(path).is_file(), (
                f"outputs[{k!r}] points to a non-existent file: {path}"
            )
        # Both leaves contributed artifacts
        assert _WATER_KEYS <= set(r.outputs.keys())
        assert r.summary["ran_water"] is True
        assert r.summary["ran_potential"] is True
        assert r.summary["water_n_artifacts"] == 6
        assert r.summary["n_artifacts"] == (
            r.summary["water_n_artifacts"]
            + r.summary["potential_n_artifacts"]
        )
        # Summary booleans mirror actual artifact presence
        assert r.summary["potential_ran_electrode"] == (
            "electrode_csv" in r.outputs
        )
        assert r.summary["potential_ran_phi_z"] == (
            "phi_z_png" in r.outputs
        )
        json.dumps(r.summary)

    def test_missing_xyz_is_file_not_found(self, tmp_path):
        r = dispatch("run_all", {
            "xyz_path": str(tmp_path / "nope.xyz"),
            "md_inp_path": str(MD_INP_PATH),
            "output_dir": str(tmp_path / "out"),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_missing_cube_files_is_file_not_found(self, tmp_path):
        """Run from an empty dir — continuous mode finds no cube
        matches for the default pattern."""
        empty = tmp_path / "empty"
        empty.mkdir()
        with _Chdir(empty):
            r = dispatch("run_all", {
                "xyz_path": str(XYZ_PATH),
                "md_inp_path": str(MD_INP_PATH),
                "output_dir": str(tmp_path / "out"),
                "md_out_path": str(MD_OUT),
            })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_cell_abc_length_is_validation(self, tmp_path):
        with _Chdir(DATA_DIR):
            r = dispatch("run_all", {
                "xyz_path": str(XYZ_PATH),
                "cell_abc": [10.0, 10.0],  # too short
                "output_dir": str(tmp_path),
                "md_out_path": str(MD_OUT),
            })
        assert not r.success
        assert r.error_type == "validation"
