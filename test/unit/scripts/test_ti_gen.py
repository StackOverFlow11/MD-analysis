"""Tests for md_analysis.scripts.TIGen — TI work directory generation."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from md_analysis.scripts.TIGen import (
    TIGenError,
    _format_target_dirname,
    _load_trajectory_cv,
    _modify_collective_block,
    _modify_inp_for_ti,
    _snap_to_nearest_frame,
    batch_generate_ti_workdirs,
    generate_ti_workdir,
)
from md_analysis.utils.RestartParser.ColvarParser import (
    ColvarRestart,
    ColvarInfo,
    ConstraintInfo,
)


# ---------------------------------------------------------------------------
# Fixtures — synthetic SG data
# ---------------------------------------------------------------------------

# A minimal SG inp file with [unit] annotations (distance case)
_SAMPLE_INP = """\
&GLOBAL
  PROJECT slowgrowth
  RUN_TYPE MD
  PRINT_LEVEL low
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS    8000
    TIMESTEP 1.0
    TEMPERATURE 300.15
  &END MD

  &CONSTRAINT
     &COLLECTIVE
       COLVAR  1
       INTERMOLECULAR .TRUE.
       TARGET [angstrom] 5.64
       TARGET_GROWTH [fs^-1*angstrom] -0.5E-3
     &END COLLECTIVE

     &LAGRANGE_MULTIPLIERS  SILENT
       COMMON_ITERATION_LEVELS 1
       FILENAME constraint_force.dat
     &END LAGRANGE_MULTIPLIERS
   &END CONSTRAINT

  &PRINT
     &TRAJECTORY  SILENT
       &EACH
         MD               5
       &END EACH
     &END TRAJECTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD Quickstep
  &DFT
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom]   10.0   10.0   20.0
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME init.xyz
      COORD_FILE_FORMAT xyz
    &END TOPOLOGY
    &COLVAR
      &DISTANCE
        ATOMS  51 119
      &END DISTANCE
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL
"""

# Inp without COORD_FILE_NAME/FORMAT (continuation run scenario)
_INP_NO_COORD = """\
&GLOBAL
  PROJECT slowgrowth
  RUN_TYPE MD
  PRINT_LEVEL low
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS    8000
    TIMESTEP 1.0
  &END MD

  &CONSTRAINT
     &COLLECTIVE
       COLVAR  1
       TARGET 1.0E+1
       TARGET_GROWTH 5.0E-5
     &END COLLECTIVE
     &LAGRANGE_MULTIPLIERS  SILENT
       FILENAME constraint_force.dat
     &END LAGRANGE_MULTIPLIERS
   &END CONSTRAINT
&END MOTION

&FORCE_EVAL
  &SUBSYS
    &TOPOLOGY
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
"""


def _make_restart(
    *,
    target_au: float = 10.0,
    target_growth_au: float = 1e-4,
    timestep_fs: float = 1.0,
    step_start: int = 0,
    total_steps: int = 8000,
) -> ColvarRestart:
    """Create a synthetic ColvarRestart for testing."""
    return ColvarRestart(
        project_name="slowgrowth",
        step_start=step_start,
        time_start_fs=0.0,
        timestep_fs=timestep_fs,
        total_steps=total_steps,
        colvars=ColvarInfo(constraints=(
            ConstraintInfo(
                colvar_id=1,
                target_au=target_au,
                target_growth_au=target_growth_au,
                intermolecular=True,
            ),
        )),
        lagrange_filename="constraint_force.dat",
        cell_abc_ang=(10.0, 10.0, 20.0),
        fixed_atom_indices=None,
    )


def _write_test_xyz(path: Path, n_frames: int = 10, step_interval: int = 5) -> None:
    """Write a minimal xyz trajectory with CP2K-style comment lines."""
    symbols = ["Cu", "Cu", "O", "H", "H"]
    positions = [
        "0.0 0.0 0.0",
        "1.8 1.8 0.0",
        "0.9 0.9 2.5",
        "0.9 0.2 3.1",
        "0.9 1.6 3.1",
    ]
    lines: list[str] = []
    for fi in range(n_frames):
        step = fi * step_interval
        time_fs = float(step)
        lines.append(f"{len(symbols)}")
        lines.append(
            f" i = {step:>8d}, time = {time_fs:>12.3f}, E = -1000.0"
        )
        for sym, pos in zip(symbols, positions):
            lines.append(f"{sym}  {pos}")
    path.write_text("\n".join(lines) + "\n")


def _write_restart_file(path: Path, restart: ColvarRestart) -> None:
    """Write a minimal .restart file that ColvarParser can parse."""
    c = restart.colvars.primary
    inter = ".TRUE." if c.intermolecular else ".FALSE."
    path.write_text(f"""\
 &GLOBAL
   PROJECT_NAME {restart.project_name}
   RUN_TYPE MD
 &END GLOBAL
 &MOTION
   &MD
     STEP_START_VAL {restart.step_start}
     TIME_START_VAL {restart.time_start_fs}
     TIMESTEP {restart.timestep_fs}
     STEPS {restart.total_steps}
   &END MD
   &CONSTRAINT
     &COLLECTIVE
       COLVAR {c.colvar_id}
       INTERMOLECULAR {inter}
       TARGET {c.target_au}
       TARGET_GROWTH {c.target_growth_au}
     &END COLLECTIVE
     &LAGRANGE_MULTIPLIERS
       FILENAME {restart.lagrange_filename}
     &END LAGRANGE_MULTIPLIERS
   &END CONSTRAINT
 &END MOTION
 &FORCE_EVAL
   &SUBSYS
     &CELL
       A {restart.cell_abc_ang[0]:.4f} 0.0 0.0
       B 0.0 {restart.cell_abc_ang[1]:.4f} 0.0
       C 0.0 0.0 {restart.cell_abc_ang[2]:.4f}
     &END CELL
   &END SUBSYS
 &END FORCE_EVAL
""")


# ---------------------------------------------------------------------------
# Tests for _modify_inp_for_ti
# ---------------------------------------------------------------------------


class TestModifyInp:
    """Tests for inp file modification logic."""

    def test_project_renamed(self):
        result = _modify_inp_for_ti(_SAMPLE_INP, 5.0, 10000)
        assert "PROJECT cMD" in result
        assert "PROJECT slowgrowth" not in result

    def test_steps_replaced(self):
        result = _modify_inp_for_ti(_SAMPLE_INP, 5.0, 15000)
        # Should appear inside &MD
        assert re.search(r"STEPS\s+15000", result)
        assert not re.search(r"STEPS\s+8000", result)

    def test_target_unit_stripped(self):
        result = _modify_inp_for_ti(_SAMPLE_INP, 5.123456789, 10000)
        # TARGET line should not have [angstrom] — extract COLLECTIVE block
        import re as _re
        coll = _re.search(
            r"&COLLECTIVE.*?&END\s+COLLECTIVE", result, _re.DOTALL | _re.IGNORECASE
        )
        assert coll is not None
        coll_text = coll.group()
        assert "[angstrom]" not in coll_text
        assert "[fs^-1*angstrom]" not in coll_text
        # Should have the bare a.u. value
        assert "5.1234567890E" in coll_text

    def test_target_growth_zeroed(self):
        result = _modify_inp_for_ti(_SAMPLE_INP, 5.0, 10000)
        # TARGET_GROWTH should be 0 without units
        assert re.search(r"TARGET_GROWTH\s+0\s*$", result, re.MULTILINE)
        assert "[fs^-1*angstrom]" not in result

    def test_target_growth_bare_au(self):
        """Handle inp with no [unit] annotation (bare a.u.)."""
        result = _modify_inp_for_ti(_INP_NO_COORD, 5.0, 10000)
        assert re.search(r"TARGET_GROWTH\s+0\s*$", result, re.MULTILINE)

    def test_topology_coord_existing(self):
        """Existing COORD_FILE_NAME/FORMAT are kept as init.xyz/XYZ."""
        result = _modify_inp_for_ti(_SAMPLE_INP, 5.0, 10000)
        assert re.search(r"COORD_FILE_NAME\s+init\.xyz", result)
        assert re.search(r"COORD_FILE_FORMAT\s+XYZ", result)

    def test_topology_coord_inserted_when_missing(self):
        """COORD_FILE_NAME/FORMAT inserted when absent."""
        result = _modify_inp_for_ti(_INP_NO_COORD, 5.0, 10000)
        assert re.search(r"COORD_FILE_NAME\s+init\.xyz", result)
        assert re.search(r"COORD_FILE_FORMAT\s+XYZ", result)

    def test_target_does_not_affect_target_growth(self):
        """TARGET replacement must not accidentally match TARGET_GROWTH."""
        result = _modify_inp_for_ti(_SAMPLE_INP, 9.99, 10000)
        # TARGET_GROWTH should be 0, not 9.99
        lines = [l.strip() for l in result.splitlines()
                 if l.strip().startswith("TARGET_GROWTH")]
        assert len(lines) == 1
        assert lines[0] == "TARGET_GROWTH 0"


class TestModifyCollectiveBlock:
    """Tests for per-block COLLECTIVE modification."""

    def test_target_and_growth_replaced(self):
        block = "       COLVAR  1\n       TARGET [angstrom] 5.64\n       TARGET_GROWTH [fs^-1*angstrom] -0.5E-3\n"
        result = _modify_collective_block(block, None, 3.14)
        assert "TARGET_GROWTH 0" in result
        assert "3.14" in result
        assert "[angstrom]" not in result

    def test_non_matching_colvar_growth_only(self):
        """Non-matching colvar_id: only TARGET_GROWTH zeroed."""
        block = "       COLVAR  2\n       TARGET 7.0\n       TARGET_GROWTH 1E-4\n"
        result = _modify_collective_block(block, 1, 3.14)
        assert "TARGET_GROWTH 0" in result
        assert "TARGET 7.0" in result  # TARGET unchanged


# ---------------------------------------------------------------------------
# Tests for frame loading and snapping
# ---------------------------------------------------------------------------


class TestFrameSnapping:
    """Tests for trajectory CV loading and nearest-frame snapping."""

    def test_load_trajectory_cv(self, tmp_path):
        """Frames loaded with correct step numbers and CV values."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5, step_interval=5)
        restart = _make_restart(target_au=10.0, target_growth_au=1e-4)
        frames = _load_trajectory_cv(xyz, restart)
        assert len(frames) == 5
        # Step numbers: 0, 5, 10, 15, 20
        steps = [s for s, _, _ in frames]
        assert steps == [0, 5, 10, 15, 20]

    def test_snap_to_nearest(self, tmp_path):
        """Snap picks the frame whose CV is closest."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5, step_interval=5)
        restart = _make_restart(
            target_au=10.0,
            target_growth_au=0.1,  # large growth for clear separation
            timestep_fs=1.0,
        )
        frames = _load_trajectory_cv(xyz, restart)
        # CV at step 0: 10.0, step 5: ~10.0 + 5*0.1*dt_au, etc.
        # Pick a target close to the CV at step 10
        cv_at_10 = frames[2][1]
        step, snapped, atoms = _snap_to_nearest_frame(frames, cv_at_10 + 0.001)
        assert step == 10
        assert snapped == pytest.approx(cv_at_10)

    def test_empty_trajectory_raises(self, tmp_path):
        """Empty trajectory raises TIGenError."""
        xyz = tmp_path / "empty.xyz"
        xyz.write_text("")
        restart = _make_restart()
        with pytest.raises(TIGenError, match="No frames found"):
            _load_trajectory_cv(xyz, restart)

    def test_inconsistent_atoms_raises(self, tmp_path):
        """Trajectory with different atom ordering raises TIGenError."""
        xyz = tmp_path / "bad.xyz"
        # Frame 0: Cu Cu O H H, Frame 1: O H H Cu Cu (reordered)
        lines = [
            "5",
            " i =        0, time =        0.000, E = -1000.0",
            "Cu  0.0 0.0 0.0",
            "Cu  1.8 1.8 0.0",
            "O   0.9 0.9 2.5",
            "H   0.9 0.2 3.1",
            "H   0.9 1.6 3.1",
            "5",
            " i =        5, time =        5.000, E = -1000.0",
            "O   0.9 0.9 2.5",
            "H   0.9 0.2 3.1",
            "H   0.9 1.6 3.1",
            "Cu  0.0 0.0 0.0",
            "Cu  1.8 1.8 0.0",
        ]
        xyz.write_text("\n".join(lines) + "\n")
        restart = _make_restart()
        with pytest.raises(TIGenError, match="Atom ordering changed"):
            _load_trajectory_cv(xyz, restart)


# ---------------------------------------------------------------------------
# Tests for generate_ti_workdir
# ---------------------------------------------------------------------------


class TestGenerateTIWorkdir:
    """Tests for single TI work directory generation."""

    def test_creates_directory_with_files(self, tmp_path):
        """Workdir contains cMD.inp and init.xyz."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
        )
        assert workdir.is_dir()
        assert (workdir / "cMD.inp").is_file()
        assert (workdir / "init.xyz").is_file()

    def test_inp_is_modified(self, tmp_path):
        """Generated cMD.inp has PROJECT=cMD and TARGET_GROWTH=0."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
        )
        content = (workdir / "cMD.inp").read_text()
        assert "PROJECT cMD" in content
        assert re.search(r"TARGET_GROWTH\s+0\s*$", content, re.MULTILINE)

    def test_auto_workdir_name(self, tmp_path):
        """Auto-generated workdir name follows ti_target_<value> pattern."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
        )
        assert workdir.name.startswith("ti_target_")

    def test_custom_workdir_name(self, tmp_path):
        """Custom workdir_name is respected."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
            workdir_name="my_ti_point",
        )
        assert workdir.name == "my_ti_point"

    def test_custom_steps(self, tmp_path):
        """Custom steps value appears in generated cMD.inp."""
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
            steps=20000,
        )
        content = (workdir / "cMD.inp").read_text()
        assert re.search(r"STEPS\s+20000", content)

    def test_init_xyz_is_valid(self, tmp_path):
        """init.xyz is a readable XYZ file with correct atom count."""
        from ase.io import read as ase_read

        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=5)
        restart = _make_restart()
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)

        workdir = generate_ti_workdir(
            inp_file, xyz, restart_file, 10.0, tmp_path / "out",
        )
        atoms = ase_read(str(workdir / "init.xyz"), format="xyz")
        assert len(atoms) == 5  # 5 atoms in _write_test_xyz


# ---------------------------------------------------------------------------
# Tests for batch_generate_ti_workdirs
# ---------------------------------------------------------------------------


class TestBatchGenerateTIWorkdirs:
    """Tests for batch TI work directory generation."""

    def _setup_files(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_test_xyz(xyz, n_frames=20, step_interval=5)
        restart = _make_restart(
            target_au=10.0,
            target_growth_au=0.01,
            total_steps=100,
        )
        restart_file = tmp_path / "sg.restart"
        _write_restart_file(restart_file, restart)
        inp_file = tmp_path / "sg.inp"
        inp_file.write_text(_SAMPLE_INP)
        return inp_file, xyz, restart_file

    def test_numeric_mode(self, tmp_path):
        """targets_au produces correct number of directories."""
        inp, xyz, restart = self._setup_files(tmp_path)
        targets = [10.0, 10.1, 10.2]
        dirs = batch_generate_ti_workdirs(
            inp, xyz, restart, tmp_path / "out",
            targets_au=targets,
        )
        assert len(dirs) == 3
        for d in dirs:
            assert (d / "cMD.inp").is_file()
            assert (d / "init.xyz").is_file()

    def test_time_mode(self, tmp_path):
        """Time-based targets produce correct number of directories."""
        inp, xyz, restart = self._setup_files(tmp_path)
        dirs = batch_generate_ti_workdirs(
            inp, xyz, restart, tmp_path / "out",
            time_initial_fs=0.0,
            time_final_fs=50.0,
            n_points=5,
        )
        assert len(dirs) == 5

    def test_both_modes_raises(self, tmp_path):
        """Specifying both modes raises TIGenError."""
        inp, xyz, restart = self._setup_files(tmp_path)
        with pytest.raises(TIGenError, match="Cannot specify both"):
            batch_generate_ti_workdirs(
                inp, xyz, restart, tmp_path / "out",
                targets_au=[10.0],
                time_initial_fs=0.0,
                time_final_fs=50.0,
                n_points=5,
            )

    def test_no_mode_raises(self, tmp_path):
        """No mode specified raises TIGenError."""
        inp, xyz, restart = self._setup_files(tmp_path)
        with pytest.raises(TIGenError, match="Must specify"):
            batch_generate_ti_workdirs(
                inp, xyz, restart, tmp_path / "out",
            )

    def test_incomplete_time_mode_raises(self, tmp_path):
        """Partial time-mode params raise TIGenError."""
        inp, xyz, restart = self._setup_files(tmp_path)
        with pytest.raises(TIGenError, match="Time mode requires"):
            batch_generate_ti_workdirs(
                inp, xyz, restart, tmp_path / "out",
                time_initial_fs=0.0,
                # missing time_final_fs and n_points
            )

    def test_all_dirs_have_project_cmd(self, tmp_path):
        """Every generated cMD.inp has PROJECT cMD."""
        inp, xyz, restart = self._setup_files(tmp_path)
        dirs = batch_generate_ti_workdirs(
            inp, xyz, restart, tmp_path / "out",
            targets_au=[10.0, 10.05],
        )
        for d in dirs:
            content = (d / "cMD.inp").read_text()
            assert "PROJECT cMD" in content


class TestFormatTargetDirname:
    """Tests for directory naming."""

    def test_positive_value(self):
        assert _format_target_dirname(3.14) == "ti_target_3.140000"

    def test_negative_value(self):
        assert _format_target_dirname(-1.5) == "ti_target_-1.500000"

    def test_zero(self):
        assert _format_target_dirname(0.0) == "ti_target_0.000000"
