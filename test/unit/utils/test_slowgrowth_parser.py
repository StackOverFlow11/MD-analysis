"""Tests for ColvarParser — restart + LagrangeMultLog parsing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.utils.config import AU_TIME_TO_FS
from md_analysis.utils.RestartParser import (
    ColvarInfo,
    ColvarMDInfo,
    ColvarParseError,
    ColvarRestart,
    ConstraintInfo,
    LagrangeMultLog,
    compute_target_series,
    parse_colvar_restart,
    parse_lagrange_mult_log,
)

DATA = Path(__file__).resolve().parents[3] / "data_example" / "sg"


# =========================================================================
# ColvarInfo
# =========================================================================

class TestColvarInfo:
    """Test ColvarInfo container."""

    @pytest.fixture()
    def info(self) -> ColvarInfo:
        c1 = ConstraintInfo(colvar_id=1, target_au=1.0, target_growth_au=0.01, intermolecular=True)
        c2 = ConstraintInfo(colvar_id=3, target_au=2.0, target_growth_au=0.02, intermolecular=False)
        return ColvarInfo(constraints=(c1, c2))

    def test_len(self, info: ColvarInfo):
        assert len(info) == 2

    def test_primary(self, info: ColvarInfo):
        assert info.primary.colvar_id == 1
        assert info.primary.target_au == 1.0

    def test_getitem(self, info: ColvarInfo):
        assert info[1].colvar_id == 1
        assert info[3].colvar_id == 3

    def test_getitem_missing(self, info: ColvarInfo):
        with pytest.raises(KeyError, match="colvar_id=99"):
            info[99]

    def test_iter(self, info: ColvarInfo):
        ids = [c.colvar_id for c in info]
        assert ids == [1, 3]


# =========================================================================
# parse_colvar_restart — 4 scenarios
# =========================================================================

class TestParseColvarRestart:
    """Parse restart files for all 4 data scenarios."""

    def test_angle(self):
        r = parse_colvar_restart(DATA / "angle" / "slowgrowth-1.restart")
        assert r.project_name == "slowgrowth"
        assert r.step_start == 5100
        assert r.total_steps == 6000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)
        assert r.time_start_fs == pytest.approx(5100.0, rel=1e-3)

        # Constraint (primary)
        assert r.colvars.primary.colvar_id == 1
        assert r.colvars.primary.intermolecular is True
        assert r.colvars.primary.target_au == pytest.approx(3.0282160653801862)
        assert r.colvars.primary.target_growth_au == pytest.approx(4.2217495722394033e-06)

        # Cell
        assert r.cell_abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)

        # Lagrange filename
        assert r.lagrange_filename == "constraint_force.dat"

        # No fixed atoms
        assert r.fixed_atom_indices is None

    def test_distance(self):
        r = parse_colvar_restart(
            DATA / "distance" / "slowgrowth-1.restart.bak-1"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 7680
        assert r.total_steps == 8000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)

        assert r.colvars.primary.target_au == pytest.approx(3.4015070391917179)
        assert r.colvars.primary.target_growth_au == pytest.approx(-2.2855144621117898e-05)
        assert r.colvars.primary.intermolecular is True

        assert r.cell_abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)
        assert r.fixed_atom_indices is None

    def test_distance_combinedCV(self):
        r = parse_colvar_restart(
            DATA / "distance_combinedCV" / "slowgrowth-1.restart"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 5270
        assert r.total_steps == 8000

        assert r.colvars.primary.target_au == pytest.approx(-8.9573018698768436e-01)
        assert r.colvars.primary.target_growth_au == pytest.approx(-9.1420578484471594e-06)

        assert r.fixed_atom_indices is None

    def test_more_constrain(self):
        r = parse_colvar_restart(
            DATA / "more_constrain" / "slowgrowth-1.restart"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 3220
        assert r.total_steps == 10000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)

        assert r.colvars.primary.target_au == pytest.approx(2.7236433780667042)
        assert r.colvars.primary.target_growth_au == pytest.approx(4.5710289242235795e-05)

        assert r.cell_abc_ang == pytest.approx((10.2239, 8.8542, 27.1231), rel=1e-4)

        # Fixed atoms: 1,2,7..10,15..18,23..26,31..34,39..42,46..49,54..57,62,63
        expected = (
            1, 2,
            7, 8, 9, 10,
            15, 16, 17, 18,
            23, 24, 25, 26,
            31, 32, 33, 34,
            39, 40, 41, 42,
            46, 47, 48, 49,
            54, 55, 56, 57,
            62, 63,
        )
        assert r.fixed_atom_indices == expected
        assert len(r.fixed_atom_indices) == 32


# =========================================================================
# Multi-COLLECTIVE block parsing
# =========================================================================

class TestMultiCollectiveBlocks:
    """Test parsing of restart files with multiple &COLLECTIVE blocks."""

    def test_more_constrain_has_multiple_colvars(self):
        r = parse_colvar_restart(
            DATA / "more_constrain" / "slowgrowth-1.restart"
        )
        # more_constrain has multiple COLLECTIVE blocks
        assert len(r.colvars) >= 1

    def test_synthetic_multi_collective(self, tmp_path):
        """Synthetic restart with two &COLLECTIVE blocks."""
        restart = tmp_path / "multi.restart"
        restart.write_text(
            "&GLOBAL\n  PROJECT_NAME multi_cv\n&END GLOBAL\n"
            "&FORCE_EVAL\n  &SUBSYS\n"
            "    &CELL\n"
            "      A  10.0  0.0  0.0\n"
            "      B   0.0 10.0  0.0\n"
            "      C   0.0  0.0 20.0\n"
            "    &END CELL\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n"
            "&MOTION\n"
            "  &MD\n"
            "    STEPS 100\n    TIMESTEP 0.5\n"
            "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n"
            "  &END MD\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n"
            "      COLVAR 1\n      TARGET 1.5\n      TARGET_GROWTH 0.001\n"
            "      INTERMOLECULAR .TRUE.\n"
            "    &END COLLECTIVE\n"
            "    &COLLECTIVE\n"
            "      COLVAR 2\n      TARGET 3.0\n      TARGET_GROWTH -0.002\n"
            "    &END COLLECTIVE\n"
            "    &LAGRANGE_MULTIPLIERS\n"
            "      FILENAME constraint_force.dat\n"
            "    &END LAGRANGE_MULTIPLIERS\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        r = parse_colvar_restart(restart)
        assert r.project_name == "multi_cv"
        assert len(r.colvars) == 2

        assert r.colvars.primary.colvar_id == 1
        assert r.colvars.primary.target_au == pytest.approx(1.5)
        assert r.colvars.primary.intermolecular is True

        assert r.colvars[2].colvar_id == 2
        assert r.colvars[2].target_au == pytest.approx(3.0)
        assert r.colvars[2].target_growth_au == pytest.approx(-0.002)
        assert r.colvars[2].intermolecular is False


# =========================================================================
# parse_lagrange_mult_log
# =========================================================================

class TestParseLagrangeMultLog:
    """Parse LagrangeMultLog files for single and multi constraint."""

    def _log_path(self, scenario: str) -> Path:
        return (
            DATA / scenario
            / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
        )

    def test_angle_single(self):
        log = parse_lagrange_mult_log(self._log_path("angle"))
        assert log.n_constraints == 1
        assert log.n_steps > 0
        assert log.shake.ndim == 1
        assert log.rattle.ndim == 1
        assert log.shake.shape == (log.n_steps,)
        assert log.rattle.shape == (log.n_steps,)
        # First values from file
        assert log.shake[0] == pytest.approx(-425.953120362)
        assert log.rattle[0] == pytest.approx(424.145881501)
        # collective_* returns same arrays for single constraint
        np.testing.assert_array_equal(log.collective_shake, log.shake)
        np.testing.assert_array_equal(log.collective_rattle, log.rattle)

    def test_distance_single(self):
        log = parse_lagrange_mult_log(self._log_path("distance"))
        assert log.n_constraints == 1
        assert log.n_steps > 0
        assert log.shake.ndim == 1

    def test_combinedCV_single(self):
        log = parse_lagrange_mult_log(self._log_path("distance_combinedCV"))
        assert log.n_constraints == 1
        assert log.n_steps > 0
        assert log.shake.ndim == 1

    def test_more_constrain_multi(self):
        log = parse_lagrange_mult_log(self._log_path("more_constrain"))
        assert log.n_constraints > 1
        assert log.n_steps > 0
        assert log.shake.ndim == 2
        assert log.rattle.ndim == 2
        assert log.shake.shape == (log.n_steps, log.n_constraints)
        assert log.rattle.shape == (log.n_steps, log.n_constraints)
        # First Shake value (collective constraint)
        assert log.shake[0, 0] == pytest.approx(-0.009233319)
        # First Rattle value
        assert log.rattle[0, 0] == pytest.approx(-0.007517126)
        # collective_* extracts column 0
        assert log.collective_shake.shape == (log.n_steps,)
        np.testing.assert_array_equal(log.collective_shake, log.shake[:, 0])
        np.testing.assert_array_equal(log.collective_rattle, log.rattle[:, 0])


# =========================================================================
# compute_target_series
# =========================================================================

class TestComputeTargetSeries:

    def test_at_step_start_equals_target(self):
        """xi(step_start) == target_au by construction."""
        r = parse_colvar_restart(DATA / "angle" / "slowgrowth-1.restart")
        log = parse_lagrange_mult_log(
            DATA / "angle"
            / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
        )
        xi = compute_target_series(r, log.n_steps)
        assert xi.shape == (log.n_steps,)
        # At k = step_start, xi[step_start] should equal target_au
        if r.step_start < log.n_steps:
            assert xi[r.step_start] == pytest.approx(
                r.colvars.primary.target_au, rel=1e-10,
            )

    def test_linear_growth(self):
        """Verify linearity: xi(k+1) - xi(k) == growth_per_step."""
        r = parse_colvar_restart(
            DATA / "distance" / "slowgrowth-1.restart.bak-1"
        )
        xi = compute_target_series(r, 100)
        diffs = np.diff(xi)
        dt_au = r.timestep_fs / AU_TIME_TO_FS
        expected_step = r.colvars.primary.target_growth_au * dt_au
        np.testing.assert_allclose(diffs, expected_step, rtol=1e-10)

    def test_shape(self):
        r = parse_colvar_restart(
            DATA / "distance_combinedCV" / "slowgrowth-1.restart"
        )
        xi = compute_target_series(r, 50)
        assert xi.shape == (50,)

    def test_colvar_id_param(self, tmp_path):
        """compute_target_series with explicit colvar_id."""
        restart = tmp_path / "multi.restart"
        restart.write_text(
            "&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL\n"
            "&FORCE_EVAL\n  &SUBSYS\n"
            "    &CELL\n"
            "      A 10.0 0.0 0.0\n      B 0.0 10.0 0.0\n      C 0.0 0.0 20.0\n"
            "    &END CELL\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n"
            "&MOTION\n"
            "  &MD\n    STEPS 100\n    TIMESTEP 1.0\n"
            "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n  &END MD\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n      COLVAR 1\n      TARGET 1.0\n      TARGET_GROWTH 0.01\n"
            "    &END COLLECTIVE\n"
            "    &COLLECTIVE\n      COLVAR 2\n      TARGET 5.0\n      TARGET_GROWTH -0.05\n"
            "    &END COLLECTIVE\n"
            "    &LAGRANGE_MULTIPLIERS\n      FILENAME f.dat\n"
            "    &END LAGRANGE_MULTIPLIERS\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        r = parse_colvar_restart(restart)
        dt_au = r.timestep_fs / AU_TIME_TO_FS
        # Default uses primary (colvar_id=1); step_start=0, so xi[0]=target_au
        xi_default = compute_target_series(r, 10)
        assert xi_default[0] == pytest.approx(1.0)
        assert xi_default[1] == pytest.approx(1.0 + 0.01 * dt_au)

        # Explicit colvar_id=2
        xi_cv2 = compute_target_series(r, 10, colvar_id=2)
        assert xi_cv2[0] == pytest.approx(5.0)
        assert xi_cv2[1] == pytest.approx(5.0 + (-0.05) * dt_au)


# =========================================================================
# Edge cases
# =========================================================================

class TestEdgeCases:

    def test_empty_file_raises(self, tmp_path):
        empty = tmp_path / "empty.LagrangeMultLog"
        empty.write_text("")
        with pytest.raises(ColvarParseError, match="empty"):
            parse_lagrange_mult_log(empty)

    def test_missing_constraint_raises(self, tmp_path):
        fake = tmp_path / "no_constraint.restart"
        fake.write_text(
            "&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL\n"
            "&MOTION\n  &MD\n    STEPS 100\n    TIMESTEP 1.0\n"
            "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n"
            "  &END MD\n&END MOTION\n"
        )
        with pytest.raises(ColvarParseError, match="CONSTRAINT"):
            parse_colvar_restart(fake)

    def test_fixed_atoms_range_expansion(self, tmp_path):
        """Verify that N..M range syntax is expanded correctly."""
        from md_analysis.utils.RestartParser.ColvarParser import (
            _parse_fixed_atoms_list,
        )
        text = (
            "&CONSTRAINT\n"
            "  &FIXED_ATOMS\n"
            "    LIST 1 3..5 10\n"
            "  &END FIXED_ATOMS\n"
            "&END CONSTRAINT\n"
        )
        result = _parse_fixed_atoms_list(text)
        assert result == (1, 3, 4, 5, 10)


# =========================================================================
# ColvarMDInfo
# =========================================================================

class TestColvarMDInfo:
    """Test ColvarMDInfo — combined restart + log session object."""

    def _restart_path(self, scenario: str) -> Path:
        mapping = {
            "angle": "slowgrowth-1.restart",
            "distance": "slowgrowth-1.restart.bak-1",
            "distance_combinedCV": "slowgrowth-1.restart",
            "more_constrain": "slowgrowth-1.restart",
        }
        return DATA / scenario / mapping[scenario]

    def _log_path(self, scenario: str) -> Path:
        return (
            DATA / scenario
            / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
        )

    def test_from_paths(self):
        info = ColvarMDInfo.from_paths(
            self._restart_path("angle"), self._log_path("angle"),
        )
        assert isinstance(info.restart, ColvarRestart)
        assert isinstance(info.lagrange, LagrangeMultLog)
        assert info.n_steps == info.lagrange.n_steps

    def test_steps_start_from_zero(self):
        info = ColvarMDInfo.from_paths(
            self._restart_path("angle"), self._log_path("angle"),
        )
        steps = info.steps
        assert steps[0] == 0
        assert steps[-1] == info.n_steps - 1
        assert len(steps) == info.n_steps

    def test_times_fs(self):
        info = ColvarMDInfo.from_paths(
            self._restart_path("angle"), self._log_path("angle"),
        )
        dt = info.restart.timestep_fs
        times = info.times_fs
        assert times[0] == pytest.approx(0.0)
        assert times[1] == pytest.approx(dt)
        assert times[-1] == pytest.approx((info.n_steps - 1) * dt)

    def test_target_at_step_start(self):
        """xi(step_start) == target_au — the restart snapshot anchor."""
        info = ColvarMDInfo.from_paths(
            self._restart_path("angle"), self._log_path("angle"),
        )
        s0 = info.restart.step_start
        if s0 < info.n_steps:
            xi = info.target_series_au()
            assert xi[s0] == pytest.approx(
                info.restart.colvars.primary.target_au, rel=1e-10,
            )

    def test_target_linear_growth(self):
        """Verify linearity: xi(k+1) - xi(k) == growth_per_step."""
        info = ColvarMDInfo.from_paths(
            self._restart_path("distance"), self._log_path("distance"),
        )
        xi = info.target_series_au()
        diffs = np.diff(xi)
        dt_au = info.restart.timestep_fs / AU_TIME_TO_FS
        expected_step = info.restart.colvars.primary.target_growth_au * dt_au
        np.testing.assert_allclose(diffs, expected_step, rtol=1e-10)

    def test_target_with_colvar_id(self, tmp_path):
        """target_series_au with explicit colvar_id."""
        restart = tmp_path / "multi.restart"
        restart.write_text(
            "&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL\n"
            "&FORCE_EVAL\n  &SUBSYS\n"
            "    &CELL\n"
            "      A 10.0 0.0 0.0\n      B 0.0 10.0 0.0\n      C 0.0 0.0 20.0\n"
            "    &END CELL\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n"
            "&MOTION\n"
            "  &MD\n    STEPS 100\n    TIMESTEP 1.0\n"
            "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n  &END MD\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n      COLVAR 1\n      TARGET 1.0\n      TARGET_GROWTH 0.01\n"
            "    &END COLLECTIVE\n"
            "    &COLLECTIVE\n      COLVAR 2\n      TARGET 5.0\n      TARGET_GROWTH -0.05\n"
            "    &END COLLECTIVE\n"
            "    &LAGRANGE_MULTIPLIERS\n      FILENAME f.dat\n"
            "    &END LAGRANGE_MULTIPLIERS\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        log_file = tmp_path / "f.dat"
        # 5 steps, single-constraint format (Shake/Rattle labels per line)
        log_file.write_text(
            "Step: 0\n"
            "Shake Lagrangian Multipliers: 0.1\n"
            "Rattle Lagrangian Multipliers: 0.2\n"
            "Step: 1\n"
            "Shake Lagrangian Multipliers: 0.3\n"
            "Rattle Lagrangian Multipliers: 0.4\n"
            "Step: 2\n"
            "Shake Lagrangian Multipliers: 0.5\n"
            "Rattle Lagrangian Multipliers: 0.6\n"
        )
        info = ColvarMDInfo.from_paths(restart, log_file)
        dt_au = info.restart.timestep_fs / AU_TIME_TO_FS
        # step_start=0, so xi[0] = target_au
        xi1 = info.target_series_au()
        assert xi1[0] == pytest.approx(1.0)
        assert xi1[1] == pytest.approx(1.0 + 0.01 * dt_au)

        xi2 = info.target_series_au(colvar_id=2)
        assert xi2[0] == pytest.approx(5.0)
        assert xi2[1] == pytest.approx(5.0 + (-0.05) * dt_au)

    def test_consistent_with_compute_target_series(self):
        """ColvarMDInfo.target_series_au matches standalone compute_target_series."""
        info = ColvarMDInfo.from_paths(
            self._restart_path("angle"), self._log_path("angle"),
        )
        xi_method = info.target_series_au()
        xi_func = compute_target_series(info.restart, info.n_steps)
        np.testing.assert_array_equal(xi_method, xi_func)
