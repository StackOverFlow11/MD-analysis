"""Tests for SlowgrowthParser — restart + LagrangeMultLog parsing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.utils.RestartParser import (
    ColvarDef,
    ConstraintInfo,
    LagrangeMultLog,
    SlowGrowthParseError,
    SlowGrowthRestart,
    compute_target_series,
    parse_lagrange_mult_log,
    parse_slowgrowth_restart,
)

DATA = Path(__file__).resolve().parents[3] / "data_example" / "sg"


# =========================================================================
# parse_slowgrowth_restart — 4 scenarios
# =========================================================================

class TestParseSlowgrowthRestart:
    """Parse restart files for all 4 data scenarios."""

    def test_angle(self):
        r = parse_slowgrowth_restart(DATA / "angle" / "slowgrowth-1.restart")
        assert r.project_name == "slowgrowth"
        assert r.step_start == 5100
        assert r.total_steps == 6000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)
        assert r.time_start_fs == pytest.approx(5100.0, rel=1e-3)

        # Constraint
        assert r.constraint.colvar_id == 1
        assert r.constraint.intermolecular is True
        assert r.constraint.target_au == pytest.approx(3.0282160653801862)
        assert r.constraint.target_growth_au == pytest.approx(4.2217495722394033e-06)

        # COLVAR
        assert r.colvar.cv_type == "ANGLE"
        assert r.colvar.atoms == (267, 119, 268)
        assert r.colvar.function is None
        assert r.colvar.components is None

        # Cell
        assert r.cell_abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)

        # Lagrange filename
        assert r.lagrange_filename == "constraint_force.dat"

        # No fixed atoms
        assert r.fixed_atom_indices is None

    def test_distance(self):
        r = parse_slowgrowth_restart(
            DATA / "distance" / "slowgrowth-1.restart.bak-1"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 7680
        assert r.total_steps == 8000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)

        assert r.constraint.target_au == pytest.approx(3.4015070391917179)
        assert r.constraint.target_growth_au == pytest.approx(-2.2855144621117898e-05)
        assert r.constraint.intermolecular is True

        assert r.colvar.cv_type == "DISTANCE"
        assert r.colvar.atoms == (51, 119)

        assert r.cell_abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)
        assert r.fixed_atom_indices is None

    def test_distance_combinedCV(self):
        r = parse_slowgrowth_restart(
            DATA / "distance_combinedCV" / "slowgrowth-1.restart"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 5270
        assert r.total_steps == 8000

        assert r.constraint.target_au == pytest.approx(-8.9573018698768436e-01)
        assert r.constraint.target_growth_au == pytest.approx(-9.1420578484471594e-06)

        assert r.colvar.cv_type == "COMBINE_COLVAR"
        assert r.colvar.function == "D1-D2"
        assert r.colvar.variables == ("D1", "D2")
        assert r.colvar.atoms is None
        assert len(r.colvar.components) == 2
        assert r.colvar.components[0].cv_type == "DISTANCE"
        assert r.colvar.components[0].atoms == (270, 229)
        assert r.colvar.components[1].cv_type == "DISTANCE"
        assert r.colvar.components[1].atoms == (228, 229)

        assert r.fixed_atom_indices is None

    def test_more_constrain(self):
        r = parse_slowgrowth_restart(
            DATA / "more_constrain" / "slowgrowth-1.restart"
        )
        assert r.project_name == "slowgrowth"
        assert r.step_start == 3220
        assert r.total_steps == 10000
        assert r.timestep_fs == pytest.approx(1.0, rel=1e-6)

        assert r.constraint.target_au == pytest.approx(2.7236433780667042)
        assert r.constraint.target_growth_au == pytest.approx(4.5710289242235795e-05)

        assert r.colvar.cv_type == "COMBINE_COLVAR"
        assert r.colvar.function == "D1-D2"
        assert r.colvar.variables == ("D1", "D2")
        assert r.colvar.components[0].atoms == (82, 80)
        assert r.colvar.components[1].atoms == (82, 221)

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
        r = parse_slowgrowth_restart(DATA / "angle" / "slowgrowth-1.restart")
        log = parse_lagrange_mult_log(
            DATA / "angle"
            / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
        )
        xi = compute_target_series(r, log.n_steps)
        assert xi.shape == (log.n_steps,)
        # At k = step_start, xi should equal target_au
        if r.step_start <= log.n_steps:
            assert xi[r.step_start - 1] == pytest.approx(
                r.constraint.target_au, rel=1e-10,
            )

    def test_linear_growth(self):
        """Verify linearity: xi(k+1) - xi(k) == target_growth_au."""
        r = parse_slowgrowth_restart(
            DATA / "distance" / "slowgrowth-1.restart.bak-1"
        )
        xi = compute_target_series(r, 100)
        diffs = np.diff(xi)
        np.testing.assert_allclose(
            diffs, r.constraint.target_growth_au, rtol=1e-10,
        )

    def test_shape(self):
        r = parse_slowgrowth_restart(
            DATA / "distance_combinedCV" / "slowgrowth-1.restart"
        )
        xi = compute_target_series(r, 50)
        assert xi.shape == (50,)


# =========================================================================
# Edge cases
# =========================================================================

class TestEdgeCases:

    def test_empty_file_raises(self, tmp_path):
        empty = tmp_path / "empty.LagrangeMultLog"
        empty.write_text("")
        with pytest.raises(SlowGrowthParseError, match="empty"):
            parse_lagrange_mult_log(empty)

    def test_missing_constraint_raises(self, tmp_path):
        fake = tmp_path / "no_constraint.restart"
        fake.write_text(
            "&GLOBAL\n  PROJECT_NAME test\n&END GLOBAL\n"
            "&MOTION\n  &MD\n    STEPS 100\n    TIMESTEP 1.0\n"
            "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n"
            "  &END MD\n&END MOTION\n"
        )
        with pytest.raises(SlowGrowthParseError, match="CONSTRAINT"):
            parse_slowgrowth_restart(fake)

    def test_fixed_atoms_range_expansion(self, tmp_path):
        """Verify that N..M range syntax is expanded correctly."""
        from md_analysis.utils.RestartParser.SlowgrowthParser import (
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
