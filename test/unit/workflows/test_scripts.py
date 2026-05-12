"""Unit tests for ``md_analysis.workflows.scripts``.

Covers all eight work-directory generation workflows in both
single-frame and batch modes:

- Bader (VASP) — single + batch
- TI (CP2K constrained-MD) — single + batch
- Potential (CP2K Hartree SP) — single + batch
- SP (CP2K DeePMD training-data SP) — single + batch

All tests run with ``generate_potcar=False`` so ``vaspkit`` is not
required. None of the workflows submit jobs; that is enforced by
the underlying ``scripts.*`` API and is documented at the workflow
boundary.

For TI we synthesise a minimal SG-style restart + ColvarMDInfo
indirectly by reusing the bundled SG fixture under
``data_example/sg/distance``. When the bundled fixture is missing
the TI tests skip.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write

from md_analysis.scripts.BaderGen import BaderGenBatchReport
from md_analysis.scripts.SpGen import SpGenBatchReport
from md_analysis.scripts.TIGen import TIGenBatchReport
from md_analysis.workflows import (
    WorkflowResult,
    require_artifacts_exist,
    run_bader_batch,
    run_bader_single,
    run_potential_batch,
    run_potential_single,
    run_sp_batch,
    run_sp_single,
    run_ti_batch,
    run_ti_single,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


def _write_minimal_xyz(xyz_path: Path, n_frames: int = 3) -> None:
    """Minimal H2 trajectory with CP2K-style ``i = ... , time = ...`` info."""
    frames = []
    for i in range(n_frames):
        atoms = Atoms("H2", positions=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.7)])
        atoms.info["i"] = i * 10
        atoms.info["time"] = float(i * 5)
        frames.append(atoms)
    write(str(xyz_path), frames, format="xyz")


def _write_slab_xyz(xyz_path: Path, n_frames: int = 3) -> None:
    """Cu-O-H slab fixture compatible with Bader IndexMapper."""
    symbols = ["Cu", "Cu", "O", "H", "H", "O", "H", "H"]
    positions_template = [
        "0.0 0.0 0.0",
        "1.8 1.8 0.0",
        "0.9 0.9 2.5",
        "0.9 0.2 3.1",
        "0.9 1.6 3.1",
        "2.7 0.9 2.5",
        "2.7 0.2 3.1",
        "2.7 1.6 3.1",
    ]
    lines: list[str] = []
    for fi in range(n_frames):
        step = fi * 5
        time_fs = float(step)
        lines.append(f"{len(symbols)}")
        lines.append(
            f" i = {step:>8d}, time = {time_fs:>12.3f}, E = -1000.0"
        )
        for sym, pos in zip(symbols, positions_template):
            lines.append(f"{sym}  {pos}")
    xyz_path.write_text("\n".join(lines) + "\n")


def _minimal_sp_template() -> str:
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
def slab_xyz(tmp_path):
    xyz = tmp_path / "traj.xyz"
    _write_slab_xyz(xyz, n_frames=3)
    return xyz


@pytest.fixture
def h2_xyz(tmp_path):
    xyz = tmp_path / "h2.xyz"
    _write_minimal_xyz(xyz, n_frames=3)
    return xyz


@pytest.fixture
def sp_inp(tmp_path):
    p = tmp_path / "sp.inp"
    p.write_text(_minimal_sp_template(), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Bader single
# ---------------------------------------------------------------------------


class TestRunBaderSingle:
    def test_creates_workdir_with_metadata(
        self, slab_xyz: Path, tmp_path: Path
    ) -> None:
        result = run_bader_single(
            xyz_path=slab_xyz,
            cell_abc=(3.6, 3.6, 10.0),
            output_dir=tmp_path / "out",
            frame=1,
            generate_potcar=False,
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "bader_single"
        assert "workdir" in result.artifacts
        require_artifacts_exist(result)
        assert result.artifacts["workdir"].is_dir()
        # POSCAR/INCAR/KPOINTS are written by the underlying generator;
        # we verify the workflow recorded the directory itself.
        assert (result.artifacts["workdir"] / "POSCAR").is_file()

        meta = result.metadata
        assert meta["frame_index"] == 1
        assert meta["mode"] == "index"
        assert meta["generate_potcar"] is False
        assert meta["source_xyz"] == str(slab_xyz)


# ---------------------------------------------------------------------------
# Bader batch
# ---------------------------------------------------------------------------


class TestRunBaderBatch:
    def test_creates_all_workdirs_with_batch_metadata(
        self, slab_xyz: Path, tmp_path: Path
    ) -> None:
        result = run_bader_batch(
            xyz_path=slab_xyz,
            cell_abc=(3.6, 3.6, 10.0),
            output_dir=tmp_path / "out",
            generate_potcar=False,
        )
        assert result.name == "bader_batch"
        require_artifacts_exist(result)
        # Three frames → three workdir_<i> artifacts
        workdir_keys = sorted(k for k in result.artifacts if k.startswith("workdir_"))
        assert workdir_keys == ["workdir_0", "workdir_1", "workdir_2"]
        for k in workdir_keys:
            assert result.artifacts[k].is_dir()

        meta = result.metadata
        # Mandatory batch contract fields
        assert meta["n_successful"] == 3
        assert meta["n_skipped"] == 0
        assert meta["n_failed"] == 0
        assert len(meta["workdir_paths"]) == 3
        assert meta["mode"] == "index"
        assert meta["frame_indices"] == [0, 1, 2]
        assert meta["n_frames"] == 3
        assert meta["generate_potcar"] is False

        assert isinstance(result.extra, BaderGenBatchReport)

    def test_to_dict_serialises_paths(
        self, slab_xyz: Path, tmp_path: Path
    ) -> None:
        result = run_bader_batch(
            xyz_path=slab_xyz,
            cell_abc=(3.6, 3.6, 10.0),
            output_dir=tmp_path / "out",
            generate_potcar=False,
        )
        dumped = result.to_dict()
        # The strongly-typed BaderGenBatchReport on extra is JSON-friendly
        json.dumps(dumped)


# ---------------------------------------------------------------------------
# TI single / batch
#   Uses the bundled SG distance fixture; skip when missing.
# ---------------------------------------------------------------------------

SG_DIR = REPO_ROOT / "data_example" / "sg" / "distance_combinedCV"
SG_INP = SG_DIR / "sg.inp"
SG_RESTART = SG_DIR / "slowgrowth-1.restart"

requires_sg_distance = pytest.mark.skipif(
    not (SG_INP.is_file() and SG_RESTART.is_file()),
    reason="SG distance_combinedCV fixture missing",
)


def _make_ti_xyz_from_restart(xyz_path: Path) -> None:
    """Write an XYZ that lines up with the SG distance_combinedCV restart.

    The restart starts at step 0 with timestep 0.5 fs by default in
    the bundled fixture; we generate frames at the steps the restart
    expects so ``_load_trajectory_cv`` does not error. The CV value
    at each frame is irrelevant for TI workdir generation as long as
    the snap-to-nearest-frame logic picks something.
    """
    # Use the simplest possible 2-atom system so iread succeeds; only
    # the i= / time= fields matter for frame sniffing. The actual
    # constraint geometry is not validated by generate_ti_workdir.
    frames = []
    for k in range(5):
        atoms = Atoms("H2", positions=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.7)])
        atoms.info["i"] = k
        atoms.info["time"] = float(k) * 0.5
        frames.append(atoms)
    write(str(xyz_path), frames, format="xyz")


@requires_sg_distance
class TestRunTiSingle:
    def test_creates_workdir(self, tmp_path: Path) -> None:
        # Build a synthetic xyz aligned with the SG restart.
        xyz = tmp_path / "sg.xyz"
        _make_ti_xyz_from_restart(xyz)

        # _load_trajectory_cv reconstructs CV from restart, not from
        # xyz positions, so any target_au within the CV range will be
        # snapped to one of the synthetic frames.
        result = run_ti_single(
            inp_path=SG_INP,
            xyz_path=xyz,
            restart_path=SG_RESTART,
            target_au=0.0,
            output_dir=tmp_path / "out",
        )
        assert result.name == "ti_single"
        require_artifacts_exist(result)
        workdir = result.artifacts["workdir"]
        assert workdir.is_dir()
        assert (workdir / "cMD.inp").is_file()
        assert (workdir / "init.xyz").is_file()
        assert result.metadata["requested_target_au"] == 0.0
        assert result.metadata["steps"] == 10000


@requires_sg_distance
class TestRunTiBatch:
    def test_targets_au_path(self, tmp_path: Path) -> None:
        xyz = tmp_path / "sg.xyz"
        _make_ti_xyz_from_restart(xyz)

        # Build targets in the small synthetic CV window. The
        # _plan_ti_targets snap maps each to the nearest synthetic frame.
        result = run_ti_batch(
            inp_path=SG_INP,
            xyz_path=xyz,
            restart_path=SG_RESTART,
            output_dir=tmp_path / "out",
            targets_au=[0.0, 0.001],
        )
        assert result.name == "ti_batch"
        require_artifacts_exist(result)
        meta = result.metadata
        assert meta["n_successful"] == 2
        assert meta["n_skipped"] == 0
        assert meta["n_failed"] == 0
        assert len(meta["workdir_paths"]) == 2
        assert meta["input_source"] == "targets_au"
        assert len(meta["requested_targets_au"]) == 2
        assert isinstance(result.extra, TIGenBatchReport)


# ---------------------------------------------------------------------------
# Potential single / batch (no upstream *_with_report)
# ---------------------------------------------------------------------------


class TestRunPotentialSingle:
    def test_creates_workdir(
        self, h2_xyz: Path, sp_inp: Path, tmp_path: Path
    ) -> None:
        result = run_potential_single(
            xyz_path=h2_xyz,
            output_dir=tmp_path / "out",
            inp_template_path=sp_inp,
            cell_abc=(5.0, 5.0, 5.0),
        )
        assert result.name == "potential_single"
        require_artifacts_exist(result)
        workdir = result.artifacts["workdir"]
        assert workdir.is_dir()
        assert (workdir / "sp.inp").is_file()
        assert (workdir / "init.xyz").is_file()
        assert result.metadata["inp_template_path_resolved_from"] == "argument"


class TestRunPotentialBatch:
    def test_creates_all_workdirs(
        self, h2_xyz: Path, sp_inp: Path, tmp_path: Path
    ) -> None:
        result = run_potential_batch(
            xyz_path=h2_xyz,
            cell_abc=(5.0, 5.0, 5.0),
            output_dir=tmp_path / "out",
            inp_template_path=sp_inp,
        )
        assert result.name == "potential_batch"
        require_artifacts_exist(result)
        meta = result.metadata
        # H2 fixture has 3 frames at index mode, default step=1.
        assert meta["n_successful"] == 3
        assert meta["n_skipped"] == 0
        assert meta["n_failed"] == 0
        assert len(meta["workdir_paths"]) == 3
        # Potential workflow has no upstream *BatchReport
        assert result.extra is None


# ---------------------------------------------------------------------------
# SP single / batch
# ---------------------------------------------------------------------------


class TestRunSpSingle:
    def test_creates_workdir(
        self, h2_xyz: Path, sp_inp: Path, tmp_path: Path
    ) -> None:
        result = run_sp_single(
            xyz_path=h2_xyz,
            cell_abc=(5.0, 5.0, 5.0),
            output_dir=tmp_path / "out",
            inp_template_path=sp_inp,
            workdir_name="sp",
        )
        assert result.name == "sp_single"
        require_artifacts_exist(result)
        workdir = result.artifacts["workdir"]
        assert workdir.is_dir()
        assert (workdir / "init.xyz").is_file()
        assert (workdir / "sp.inp").is_file()


class TestRunSpBatch:
    def test_creates_all_workdirs(
        self, h2_xyz: Path, sp_inp: Path, tmp_path: Path
    ) -> None:
        result = run_sp_batch(
            xyz_path=h2_xyz,
            cell_abc=(5.0, 5.0, 5.0),
            output_dir=tmp_path / "out",
            inp_template_path=sp_inp,
        )
        assert result.name == "sp_batch"
        require_artifacts_exist(result)
        meta = result.metadata
        assert meta["n_successful"] == 3
        assert meta["n_skipped"] == 0
        assert meta["n_failed"] == 0
        assert len(meta["workdir_paths"]) == 3
        assert meta["frame_indices"] == [0, 1, 2]
        assert isinstance(result.extra, SpGenBatchReport)


# ---------------------------------------------------------------------------
# Cross-cutting: artifact map shape contract
# ---------------------------------------------------------------------------


class TestBatchArtifactContract:
    def test_keys_match_workdir_paths_order(
        self, slab_xyz: Path, tmp_path: Path
    ) -> None:
        """``workdir_<i>`` artifact keys must align with
        ``metadata['workdir_paths']`` indices so failing-artifact
        diagnostics map cleanly between the two views.
        """
        result = run_bader_batch(
            xyz_path=slab_xyz,
            cell_abc=(3.6, 3.6, 10.0),
            output_dir=tmp_path / "out",
            generate_potcar=False,
        )
        n = result.metadata["n_successful"]
        for i in range(n):
            assert (
                str(result.artifacts[f"workdir_{i}"])
                == result.metadata["workdir_paths"][i]
            )
