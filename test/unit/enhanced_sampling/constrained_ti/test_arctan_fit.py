"""Tests for the arctan SEM extrapolation module (_arctan_fit.py)."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.analysis._arctan_fit import (
    fit_arctan_sem,
)
from md_analysis.enhanced_sampling.constrained_ti.models import ArctanFitResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_synthetic_arctan(
    a: float = 0.001,
    b: float = 0.05,
    n_points: int = 30,
    noise_std: float = 1e-5,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a clean arctan SEM curve with small additive noise."""
    rng = np.random.default_rng(seed)
    block_sizes = np.arange(1, n_points + 1)
    sem_curve = a * np.arctan(b * block_sizes) + rng.normal(0, noise_std, n_points)
    sem_curve = np.maximum(sem_curve, 1e-15)  # ensure positive
    return block_sizes, sem_curve


def _make_ar1(n: int, phi: float, sigma: float = 1.0, seed: int = 42) -> np.ndarray:
    """Generate AR(1) time series: x_t = phi * x_{t-1} + eps."""
    rng = np.random.default_rng(seed)
    eps = rng.normal(0, sigma * np.sqrt(1 - phi**2), n)
    series = np.empty(n)
    series[0] = eps[0]
    for t in range(1, n):
        series[t] = phi * series[t - 1] + eps[t]
    return series


def _block_average_sem_curve(
    series: np.ndarray, block_sizes: np.ndarray
) -> np.ndarray:
    """Compute SEM(B) curve from a time series."""
    n = len(series)
    sem = np.empty(len(block_sizes))
    for i, bs in enumerate(block_sizes):
        n_blocks = n // bs
        if n_blocks < 2:
            sem[i] = np.nan
            continue
        trimmed = series[: n_blocks * bs].reshape(n_blocks, bs)
        block_means = trimmed.mean(axis=1)
        sem[i] = np.std(block_means, ddof=1) / np.sqrt(n_blocks)
    return sem


# ---------------------------------------------------------------------------
# P0 tests
# ---------------------------------------------------------------------------


class TestSyntheticRecovery:
    """#1: Recover A, B from synthetic arctan curve."""

    def test_parameter_recovery(self) -> None:
        a_true, b_true = 0.001, 0.05
        bs, sem = _make_synthetic_arctan(a=a_true, b=b_true, noise_std=1e-6)
        result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.01)
        assert result is not None
        assert result.reliable == True  # noqa: E712 (np.bool_ safe)
        assert result.r2 > 0.99
        assert result.A == pytest.approx(a_true, rel=0.05)
        assert result.B == pytest.approx(b_true, rel=0.05)
        assert result.sem_asymptote == pytest.approx(a_true * np.pi / 2, rel=0.05)


class TestAR1RealSEM:
    """#2: Arctan on real block-average SEM curve from AR(1).

    Note: with dense sampling, the arctan model may not always achieve
    R² ≥ 0.95 for AR(1) data.  We test that it returns a valid result
    and that the asymptote is in the right ballpark when reliable.
    """

    def test_asymptote_vs_theory(self) -> None:
        phi = 0.9
        n = 50000  # Long series for better SEM curve
        series = _make_ar1(n, phi, seed=42)
        sigma = np.std(series, ddof=0)
        tau_theory = (1 + phi) / (2 * (1 - phi))  # = 9.5
        sem_theory = sigma * np.sqrt(2 * tau_theory / n)

        from md_analysis.enhanced_sampling.constrained_ti.analysis.block_average import (
            _generate_block_sizes,
        )

        bs = _generate_block_sizes(n, min_blocks=4)
        sem_curve = _block_average_sem_curve(series, bs)

        result = fit_arctan_sem(bs, sem_curve, n_total=n, sigma_series=sigma)
        assert result is not None
        # Asymptote should be in the right ballpark regardless of R²
        assert result.sem_asymptote == pytest.approx(sem_theory, rel=0.20)


class TestAllZeroSEM:
    """#3: All-zero SEM → ArctanFitResult(reliable=False)."""

    def test_zero_sem(self) -> None:
        bs = np.arange(1, 20)
        sem = np.zeros(19)
        result = fit_arctan_sem(bs, sem, n_total=1000, sigma_series=1.0)
        assert result is not None
        assert isinstance(result, ArctanFitResult)
        assert result.reliable == False  # noqa: E712 (np.bool_ safe)


class TestMonotonicallyDecreasing:
    """#4: Monotonically decreasing SEM → reliable=False."""

    def test_decreasing(self) -> None:
        bs = np.arange(1, 20)
        sem = np.linspace(0.01, 0.001, 19)  # strictly decreasing
        result = fit_arctan_sem(bs, sem, n_total=1000, sigma_series=1.0)
        assert result is not None
        assert isinstance(result, ArctanFitResult)
        assert result.reliable == False  # noqa: E712 (np.bool_ safe)


class TestTooFewPoints:
    """#5: Fewer than min_points → returns None."""

    def test_too_few(self) -> None:
        bs = np.array([1, 2, 4])
        sem = np.array([0.001, 0.0015, 0.0016])
        result = fit_arctan_sem(bs, sem, n_total=100, sigma_series=1.0, min_points=5)
        assert result is None


class TestScipyUnavailable:
    """#6: scipy not importable → returns None."""

    def test_no_scipy(self) -> None:
        import builtins

        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if "scipy" in name:
                raise ImportError("mocked")
            return original_import(name, *args, **kwargs)

        bs, sem = _make_synthetic_arctan()
        with patch("builtins.__import__", side_effect=mock_import):
            result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.01)
        assert result is None


class TestCurveFitRuntimeError:
    """#7: curve_fit raises RuntimeError → ArctanFitResult(reliable=False)."""

    def test_runtime_error(self) -> None:
        bs, sem = _make_synthetic_arctan()
        with patch(
            "scipy.optimize.curve_fit", side_effect=RuntimeError("no convergence")
        ):
            result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.01)
        assert result is not None
        assert isinstance(result, ArctanFitResult)
        assert result.reliable is False
        assert np.isnan(result.sem_asymptote)


# ---------------------------------------------------------------------------
# P1 tests
# ---------------------------------------------------------------------------


class TestNaNInSEMCurve:
    """#8: SEM curve with NaN values."""

    def test_nan_filtered(self) -> None:
        bs, sem = _make_synthetic_arctan(n_points=20, noise_std=1e-6)
        sem[5] = np.nan
        sem[10] = np.nan
        sem[15] = np.nan
        result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.01)
        # 17 valid points > min_points(5), so should still fit
        assert result is not None


class TestSmallSEMValues:
    """#9: Very small SEM values (~1e-10) — numerical stability."""

    def test_tiny_sem(self) -> None:
        bs, sem = _make_synthetic_arctan(a=1e-10, b=0.05, noise_std=1e-12)
        result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=1e-8)
        assert result is not None
        # Should either be reliable or gracefully unreliable, not crash


class TestPcovSingular:
    """#11: pcov diagonal contains inf → reliable=False."""

    def test_inf_pcov(self) -> None:
        bs, sem = _make_synthetic_arctan()
        inf_pcov = np.array([[np.inf, 0], [0, np.inf]])
        with patch(
            "scipy.optimize.curve_fit",
            return_value=(np.array([0.001, 0.05]), inf_pcov),
        ):
            result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.01)
        assert result is not None
        assert result.reliable is False


class TestTauCorrImplied:
    """#12: tau_corr_implied matches known theoretical value."""

    def test_tau_from_synthetic(self) -> None:
        # Known: SEM = sigma * sqrt(2*tau/N) → tau = N * SEM^2 / (2 * sigma^2)
        a_true = 0.001
        sem_asym = a_true * np.pi / 2
        sigma = 0.01
        n_total = 2000
        tau_expected = n_total * sem_asym**2 / (2 * sigma**2)

        bs, sem = _make_synthetic_arctan(a=a_true, b=0.05, noise_std=1e-6)
        result = fit_arctan_sem(bs, sem, n_total=n_total, sigma_series=sigma)
        assert result is not None
        assert result.reliable is True
        assert result.tau_corr_implied == pytest.approx(tau_expected, rel=0.15)


class TestSigmaZero:
    """#13: sigma_series = 0 → tau_corr_implied = nan."""

    def test_sigma_zero(self) -> None:
        bs, sem = _make_synthetic_arctan(noise_std=1e-6)
        result = fit_arctan_sem(bs, sem, n_total=2000, sigma_series=0.0)
        assert result is not None
        assert np.isnan(result.tau_corr_implied)
