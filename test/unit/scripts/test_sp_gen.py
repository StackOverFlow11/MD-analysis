"""Tests for SpGen agent-facing wrapper + contract-backed dispatch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.scripts.SpGen import (
    SpGenBatchReport,
    SpGenError,
    generate_sp_batch_with_report,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _write_trajectory(xyz_path: Path, n_frames: int = 3) -> None:
    """Write a tiny multi-frame XYZ trajectory with CP2K-style info fields."""
    frames = []
    for i in range(n_frames):
        atoms = Atoms("H2", positions=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.7)])
        atoms.info["i"] = i * 10
        atoms.info["time"] = float(i * 5)
        frames.append(atoms)
    write(str(xyz_path), frames, format="xyz")


def _minimal_sp_template() -> str:
    """Minimum CP2K sp.inp the modify_inp_for_sp helper will accept."""
    return (
        "&FORCE_EVAL\n"
        "  &SUBSYS\n"
        "    &CELL\n"
        "      ABC 1.0 1.0 1.0\n"
        "    &END CELL\n"
        "    &TOPOLOGY\n"
        "      COORD_FILE_NAME init.xyz\n"
        "      COORD_FILE_FORMAT XYZ\n"
        "    &END TOPOLOGY\n"
        "  &END SUBSYS\n"
        "&END FORCE_EVAL\n"
    )


@pytest.fixture
def sp_fixture(tmp_path):
    xyz = tmp_path / "md.xyz"
    tmpl = tmp_path / "sp.inp"
    outdir = tmp_path / "out"
    _write_trajectory(xyz, n_frames=3)
    tmpl.write_text(_minimal_sp_template(), encoding="utf-8")
    return {"xyz": xyz, "tmpl": tmpl, "outdir": outdir, "tmp": tmp_path}


# ---------------------------------------------------------------------------
# Wrapper-level behaviour
# ---------------------------------------------------------------------------


class TestGenerateSpBatchWithReport:
    def test_success_path_writes_expected_files(self, sp_fixture):
        r = generate_sp_batch_with_report(
            xyz_path=sp_fixture["xyz"],
            cell_abc=(5.0, 5.0, 10.0),
            output_dir=sp_fixture["outdir"],
            inp_template_path=sp_fixture["tmpl"],
        )
        assert isinstance(r, SpGenBatchReport)
        assert r.n_frames == 3
        assert len(r.workdirs) == 3
        for wd in r.workdirs:
            assert (wd / "init.xyz").is_file()
            assert (wd / "sp.inp").is_file()

    def test_report_to_dict_is_json_serializable(self, sp_fixture):
        r = generate_sp_batch_with_report(
            xyz_path=sp_fixture["xyz"],
            cell_abc=(5.0, 5.0, 10.0),
            output_dir=sp_fixture["outdir"],
            inp_template_path=sp_fixture["tmpl"],
        )
        # Must not raise — round trip via json
        json.dumps(r.to_dict())

    def test_missing_xyz_raises_file_not_found(self, sp_fixture):
        with pytest.raises(FileNotFoundError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["tmp"] / "missing.xyz",
                cell_abc=(5.0, 5.0, 10.0),
                output_dir=sp_fixture["outdir"],
                inp_template_path=sp_fixture["tmpl"],
            )
        assert not sp_fixture["outdir"].exists()

    def test_bad_cell_abc_length_raises_value_error(self, sp_fixture):
        with pytest.raises(ValueError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["xyz"],
                cell_abc=(5.0, 5.0),  # only 2 components
                output_dir=sp_fixture["outdir"],
                inp_template_path=sp_fixture["tmpl"],
            )
        assert not sp_fixture["outdir"].exists()

    def test_non_positive_cell_component_raises_value_error(self, sp_fixture):
        with pytest.raises(ValueError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["xyz"],
                cell_abc=(5.0, 0.0, 10.0),
                output_dir=sp_fixture["outdir"],
                inp_template_path=sp_fixture["tmpl"],
            )
        assert not sp_fixture["outdir"].exists()

    def test_no_template_no_config_fallback_raises_sp_gen_error(
        self, sp_fixture, monkeypatch,
    ):
        """SpGenError when neither inp_template_path nor config supplies one."""
        monkeypatch.setattr(
            "md_analysis.scripts.SpGen.get_config", lambda *a, **kw: None,
        )
        with pytest.raises(SpGenError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["xyz"],
                cell_abc=(5.0, 5.0, 10.0),
                output_dir=sp_fixture["outdir"],
                inp_template_path=None,
            )

    def test_missing_template_file_raises_file_not_found(self, sp_fixture):
        with pytest.raises(FileNotFoundError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["xyz"],
                cell_abc=(5.0, 5.0, 10.0),
                output_dir=sp_fixture["outdir"],
                inp_template_path=sp_fixture["tmp"] / "nope.inp",
            )

    def test_missing_script_path_raises_file_not_found(self, sp_fixture):
        with pytest.raises(FileNotFoundError):
            generate_sp_batch_with_report(
                xyz_path=sp_fixture["xyz"],
                cell_abc=(5.0, 5.0, 10.0),
                output_dir=sp_fixture["outdir"],
                inp_template_path=sp_fixture["tmpl"],
                script_path=sp_fixture["tmp"] / "nope.sh",
            )

    def test_existing_same_name_workdir_is_overwritten(self, sp_fixture):
        """Documented behavior: init.xyz / sp.inp are overwritten on re-run."""
        out = sp_fixture["outdir"]
        r1 = generate_sp_batch_with_report(
            xyz_path=sp_fixture["xyz"],
            cell_abc=(5.0, 5.0, 10.0),
            output_dir=out,
            inp_template_path=sp_fixture["tmpl"],
        )
        mtime_before = (r1.workdirs[0] / "init.xyz").stat().st_mtime_ns
        # Rewrite must not error even though workdirs already exist.
        r2 = generate_sp_batch_with_report(
            xyz_path=sp_fixture["xyz"],
            cell_abc=(5.0, 5.0, 10.0),
            output_dir=out,
            inp_template_path=sp_fixture["tmpl"],
        )
        assert r1.workdirs == r2.workdirs
        mtime_after = (r2.workdirs[0] / "init.xyz").stat().st_mtime_ns
        assert mtime_after >= mtime_before


# ---------------------------------------------------------------------------
# Contract / dispatch-level behaviour
# ---------------------------------------------------------------------------


class TestSpGenBatchSchema:
    def test_schema_keys(self):
        s = get_task_schema("sp_gen_batch")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required_fields(self):
        s = get_task_schema("sp_gen_batch")
        assert set(s["parameters"]["required"]) == {
            "xyz_path", "cell_abc", "output_dir",
        }

    def test_includes_all_wrapper_params(self):
        s = get_task_schema("sp_gen_batch")
        props = s["parameters"]["properties"]
        expected = {
            "xyz_path", "cell_abc", "output_dir",
            "inp_template_path",
            "mode", "frame_start", "frame_end", "frame_step",
            "time_start_fs", "time_end_fs", "time_step_fs",
            "script_path", "verbose",
        }
        assert expected == set(props.keys())

    def test_mode_enum(self):
        s = get_task_schema("sp_gen_batch")
        assert set(
            s["parameters"]["properties"]["mode"].get("enum", [])
        ) == {"index", "time"}

    def test_cell_abc_is_length_3_numeric(self):
        s = get_task_schema("sp_gen_batch")
        ca = s["parameters"]["properties"]["cell_abc"]
        assert ca["type"] == "array"
        assert ca["minItems"] == 3 and ca["maxItems"] == 3
        assert ca["items"]["type"] == "number"

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["sp_gen_batch"].contract is not None


class TestSpGenBatchDispatch:
    def test_missing_xyz_returns_file_not_found(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["tmp"] / "missing.xyz"),
            "cell_abc": [5.0, 5.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmpl"]),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_no_template_no_config_returns_validation(
        self, sp_fixture, monkeypatch,
    ):
        monkeypatch.setattr(
            "md_analysis.scripts.SpGen.get_config", lambda *a, **kw: None,
        )
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 5.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_missing_template_returns_file_not_found(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 5.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmp"] / "nope.inp"),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_mode_returns_validation(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 5.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmpl"]),
            "mode": "bogus",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_bad_cell_abc_length_returns_validation(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 5.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmpl"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_non_positive_cell_abc_returns_validation(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 0.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmpl"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_success_path_returns_summary_metrics(self, sp_fixture):
        r = dispatch("sp_gen_batch", {
            "xyz_path": str(sp_fixture["xyz"]),
            "cell_abc": [5.0, 5.0, 10.0],
            "output_dir": str(sp_fixture["outdir"]),
            "inp_template_path": str(sp_fixture["tmpl"]),
        })
        assert r.success, r.errors
        for k in ("n_frames", "frame_indices", "steps", "times_fs"):
            assert k in r.summary
        assert r.summary["n_frames"] == 3
        # At least one workdir surfaced in outputs.
        assert any(k.startswith("workdir_") for k in r.outputs)
        # summary is JSON-serializable
        json.dumps(r.summary)
