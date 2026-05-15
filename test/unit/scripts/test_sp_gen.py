"""Tests for the SpGen batch wrapper (``generate_sp_batch_with_report``)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write

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
