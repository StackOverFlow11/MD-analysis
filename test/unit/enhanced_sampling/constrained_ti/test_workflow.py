"""Tests for workflow: analyze_standalone, analyze_ti, validation."""

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.models import (
    ArctanFitResult,
    BlockAverageResult,
    InsufficientSamplingError,
)
from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    analyze_standalone,
    analyze_ti,
)


def make_iid(n, mu=0.0, sigma=1.0, seed=42):
    rng = np.random.default_rng(seed)
    return rng.normal(mu, sigma, n)


def make_ar1(n, phi, sigma=1.0, seed=42):
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + rng.normal(0, sigma * np.sqrt(1 - phi**2))
    return x


def make_drifting(n, slope=0.01, sigma=1.0, seed=42):
    rng = np.random.default_rng(seed)
    return slope * np.arange(n) + rng.normal(0, sigma, n)


# ---------------------------------------------------------------------------
# analyze_standalone
# ---------------------------------------------------------------------------

class TestAnalyzeStandalone:
    def test_no_sem_target_passed_is_none(self):
        series = make_iid(2000)
        report = analyze_standalone(series, sem_target=None)
        assert report.passed is None
        assert report.sem_max is None
        assert report.point_index is None

    def test_with_sem_target_pass(self):
        series = make_iid(10000, sigma=1.0, seed=0)
        # Very loose target — should pass
        report = analyze_standalone(series, sem_target=10.0)
        assert report.passed == True  # noqa: E712 (np.bool_ safe)

    def test_with_sem_target_fail(self):
        series = make_iid(500, sigma=1.0)
        # Extremely tight target — should fail on SEM
        report = analyze_standalone(series, sem_target=1e-10)
        assert report.passed == False  # noqa: E712
        assert len(report.failure_reasons) > 0

    def test_all_diagnostics_populated(self):
        series = make_iid(1000)
        report = analyze_standalone(series)
        assert report.autocorr is not None
        assert report.block_avg is not None
        assert report.running_avg is not None
        assert report.geweke is not None
        assert report.lambda_mean == pytest.approx(np.mean(series), abs=0.01)

    def test_failure_reasons_empty_when_no_target(self):
        """With no target, may still have warnings but passed is None."""
        series = make_iid(5000)
        report = analyze_standalone(series)
        assert report.passed is None

    def test_equilibration_trim(self):
        series = np.concatenate([np.ones(100) * 1000, make_iid(1000)])
        report = analyze_standalone(series, equilibration=100)
        # After trimming the spike, mean should be near 0
        assert abs(report.lambda_mean) < 1.0

    def test_equilibration_exceeds_length(self):
        series = make_iid(50)
        with pytest.raises(InsufficientSamplingError):
            analyze_standalone(series, equilibration=100)

    def test_nan_handling(self):
        series = make_iid(1000)
        series_with_nan = series.copy()
        series_with_nan[-5:] = np.nan  # 0.5% NaN
        report = analyze_standalone(series_with_nan)
        # Should run without error, with warning
        assert any("NaN" in r for r in report.failure_reasons)

    def test_too_many_nans(self):
        series = np.full(100, np.nan)
        with pytest.raises(InsufficientSamplingError):
            analyze_standalone(series)

    def test_short_series_warning(self):
        series = make_iid(50)
        report = analyze_standalone(series)
        assert any("short" in r.lower() or "unreliable" in r.lower()
                    for r in report.failure_reasons)


# ---------------------------------------------------------------------------
# analyze_ti
# ---------------------------------------------------------------------------

class TestAnalyzeTI:
    def test_clean_ar1_all_pass(self):
        xi = np.array([1.0, 2.0, 3.0])
        # Use longer series and lower phi for robust convergence
        series_list = [make_ar1(10000, 0.3, seed=10 + i) for i in range(3)]
        report = analyze_ti(xi, series_list, dt=1.0, epsilon_tol_kcal=100.0)
        assert report.all_passed == True  # noqa: E712
        assert len(report.failing_indices) == 0
        assert len(report.point_reports) == 3

    def test_drifting_fails(self):
        xi = np.array([1.0, 2.0, 3.0])
        series_list = [
            make_iid(5000, seed=0),
            make_drifting(5000, slope=0.1, seed=1),  # Should fail
            make_iid(5000, seed=2),
        ]
        report = analyze_ti(xi, series_list, dt=1.0, epsilon_tol_kcal=100.0)
        # At least one point should have failure reasons (drift or Geweke)
        any_failures = any(
            len(r.failure_reasons) > 0 for r in report.point_reports
        )
        assert any_failures or not report.all_passed

    def test_unequal_lengths(self):
        xi = np.array([1.0, 2.0])
        series_list = [make_iid(3000, seed=0), make_iid(1000, seed=1)]
        report = analyze_ti(xi, series_list, dt=1.0, epsilon_tol_kcal=100.0)
        assert len(report.point_reports) == 2

    def test_k_less_than_2_raises(self):
        xi = np.array([1.0])
        with pytest.raises(ValueError, match="at least 2"):
            analyze_ti(xi, [make_iid(1000)], dt=1.0)

    def test_non_monotonic_xi_raises(self):
        xi = np.array([3.0, 1.0, 2.0])
        series_list = [make_iid(1000, seed=i) for i in range(3)]
        with pytest.raises(ValueError, match="monotonic"):
            analyze_ti(xi, series_list, dt=1.0)

    def test_epsilon_tol_unit_conversion(self):
        xi = np.array([1.0, 2.0])
        series_list = [make_iid(2000, seed=i) for i in range(2)]
        report = analyze_ti(xi, series_list, dt=1.0, epsilon_tol_kcal=1.0)
        # epsilon_tol_au should be ~1/627.5
        assert report.epsilon_tol_au == pytest.approx(1.0 / 627.509474, rel=1e-6)

    def test_report_fields(self):
        xi = np.array([1.0, 2.0, 3.0])
        series_list = [make_iid(2000, seed=i) for i in range(3)]
        report = analyze_ti(xi, series_list, dt=1.0)
        assert report.xi_values.shape == (3,)
        assert report.weights.shape == (3,)
        assert report.forces.shape == (3,)
        assert report.force_errors.shape == (3,)
        assert isinstance(report.delta_A, float)
        assert isinstance(report.sigma_A, float)


# ---------------------------------------------------------------------------
# Arctan 3-tier SEM selection tests
# ---------------------------------------------------------------------------


class TestArctanIntegration:
    """Tests for arctan integration in the workflow."""

    def test_arctan_field_exists(self):
        """#1: block_avg.arctan field populated for AR(1) data."""
        series = make_ar1(5000, 0.9, seed=42)
        report = analyze_standalone(series)
        assert report.block_avg.arctan is not None

    def test_path1_arctan_reliable_uses_asymptote(self):
        """#2: When arctan is reliable, sem_final == arctan.sem_asymptote."""
        series = make_ar1(5000, 0.9, seed=42)
        report = analyze_standalone(series)
        arctan = report.block_avg.arctan
        if arctan is not None and arctan.reliable:
            assert report.sem_final == pytest.approx(
                arctan.sem_asymptote, rel=1e-10
            )

    def test_path2_arctan_unreliable_plateau(self):
        """#3: Arctan unreliable + plateau → sem_final = sem_plateau.

        Uses mock BlockAverageResult since this path is hard to trigger
        naturally.
        """
        from unittest.mock import patch

        unreliable_arctan = ArctanFitResult(
            sem_asymptote=999.0,
            A=999.0,
            B=999.0,
            r2=0.5,
            reliable=False,
            tau_corr_implied=float("nan"),
            fit_curve=None,
        )

        def patched_analyze_block_average(series, **kwargs):
            real_module = __import__(
                "md_analysis.enhanced_sampling.constrained_ti.analysis.block_average",
                fromlist=["analyze_block_average"],
            )
            # Call real function with dense_sampling=False to get real plateau result
            real_result = real_module.analyze_block_average(
                series, dense_sampling=False, **{
                    k: v for k, v in kwargs.items() if k != "dense_sampling"
                }
            )
            # Construct result with unreliable arctan but real plateau
            return BlockAverageResult(
                block_sizes=real_result.block_sizes,
                sem_curve=real_result.sem_curve,
                sem_plateau=real_result.sem_plateau,
                sem_at_max_B=real_result.sem_at_max_B,
                plateau_rtol=real_result.plateau_rtol,
                plateau_reached=True,  # Force plateau reached
                cross_valid_ok=real_result.cross_valid_ok,
                passed=real_result.passed,
                arctan=unreliable_arctan,
            )

        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.workflow.analyze_block_average",
            side_effect=patched_analyze_block_average,
        ):
            series = make_ar1(5000, 0.3, seed=42)
            report = analyze_standalone(series)

        assert report.sem_final != 999.0  # Not the unreliable arctan value
        assert report.block_avg.arctan is not None
        assert report.block_avg.arctan.reliable is False

    def test_path3_all_fallback_to_acf(self):
        """#4: Arctan=None + no plateau → sem_final = sem_auto."""
        from unittest.mock import patch

        def patched_analyze_block_average(series, **kwargs):
            real_module = __import__(
                "md_analysis.enhanced_sampling.constrained_ti.analysis.block_average",
                fromlist=["analyze_block_average"],
            )
            real_result = real_module.analyze_block_average(
                series, dense_sampling=False, **{
                    k: v for k, v in kwargs.items() if k != "dense_sampling"
                }
            )
            return BlockAverageResult(
                block_sizes=real_result.block_sizes,
                sem_curve=real_result.sem_curve,
                sem_plateau=real_result.sem_plateau,
                sem_at_max_B=real_result.sem_at_max_B,
                plateau_rtol=real_result.plateau_rtol,
                plateau_reached=False,  # Force no plateau
                cross_valid_ok=real_result.cross_valid_ok,
                passed=real_result.passed,
                arctan=None,  # No arctan
            )

        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.workflow.analyze_block_average",
            side_effect=patched_analyze_block_average,
        ):
            series = make_ar1(5000, 0.3, seed=42)
            report = analyze_standalone(series)

        # sem_final should be sem_auto (ACF fallback)
        assert report.sem_final == pytest.approx(
            report.autocorr.sem_auto, rel=1e-10
        )

    def test_failure_reasons_contain_arctan_message(self):
        """#5: When arctan fails, failure_reasons mention it."""
        from unittest.mock import patch

        unreliable_arctan = ArctanFitResult(
            sem_asymptote=0.001,
            A=0.001,
            B=0.01,
            r2=0.5,
            reliable=False,
            tau_corr_implied=float("nan"),
            fit_curve=None,
        )

        def patched_analyze_block_average(series, **kwargs):
            real_module = __import__(
                "md_analysis.enhanced_sampling.constrained_ti.analysis.block_average",
                fromlist=["analyze_block_average"],
            )
            real_result = real_module.analyze_block_average(
                series, dense_sampling=False, **{
                    k: v for k, v in kwargs.items() if k != "dense_sampling"
                }
            )
            return BlockAverageResult(
                block_sizes=real_result.block_sizes,
                sem_curve=real_result.sem_curve,
                sem_plateau=real_result.sem_plateau,
                sem_at_max_B=real_result.sem_at_max_B,
                plateau_rtol=real_result.plateau_rtol,
                plateau_reached=False,
                cross_valid_ok=real_result.cross_valid_ok,
                passed=real_result.passed,
                arctan=unreliable_arctan,
            )

        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.workflow.analyze_block_average",
            side_effect=patched_analyze_block_average,
        ):
            series = make_ar1(5000, 0.3, seed=42)
            report = analyze_standalone(series)

        arctan_msgs = [r for r in report.failure_reasons if "Arctan" in r or "arctan" in r]
        assert len(arctan_msgs) > 0

    def test_clean_ar1_all_pass_regression(self):
        """#6 (P0): Regression — test_clean_ar1_all_pass behavior preserved."""
        xi = np.array([1.0, 2.0, 3.0])
        series_list = [make_ar1(10000, 0.3, seed=10 + i) for i in range(3)]
        report = analyze_ti(xi, series_list, dt=1.0, epsilon_tol_kcal=100.0)
        assert report.all_passed is True
        # SEM values should be in reasonable range
        for r in report.point_reports:
            assert 0 < r.sem_final < 1.0
