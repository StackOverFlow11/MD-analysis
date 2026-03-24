"""Tests for dense block-size sampling and arctan integration in block_average.py."""

from __future__ import annotations

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.analysis.block_average import (
    _generate_block_sizes,
    _generate_block_sizes_pow2,
    analyze_block_average,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ar1(n: int, phi: float, sigma: float = 1.0, seed: int = 42) -> np.ndarray:
    """Generate AR(1) time series."""
    rng = np.random.default_rng(seed)
    eps = rng.normal(0, sigma * np.sqrt(1 - phi**2), n)
    series = np.empty(n)
    series[0] = eps[0]
    for t in range(1, n):
        series[t] = phi * series[t - 1] + eps[t]
    return series


# ---------------------------------------------------------------------------
# P0 tests: _generate_block_sizes
# ---------------------------------------------------------------------------


class TestGenerateBlockSizes:
    """Block-size generation strategy tests."""

    def test_n2000_basic_properties(self) -> None:
        """#1: N=2000 → includes 1..20, ~28-38 total, strictly increasing."""
        sizes = _generate_block_sizes(2000)
        assert all(s in sizes for s in range(1, 21))
        assert 28 <= len(sizes) <= 38
        # Strictly increasing
        assert np.all(np.diff(sizes) > 0)

    def test_n40_small(self) -> None:
        """#2: N=40, max_b=10."""
        sizes = _generate_block_sizes(40)
        assert sizes[0] == 1
        assert sizes[-1] == 10
        assert np.array_equal(sizes, np.arange(1, 11))

    def test_n4_extreme(self) -> None:
        """#3: N=4, max_b=1 → [1]."""
        sizes = _generate_block_sizes(4)
        assert np.array_equal(sizes, np.array([1]))

    def test_n0_edge(self) -> None:
        """#4: N=0 → [1], no crash."""
        sizes = _generate_block_sizes(0)
        assert np.array_equal(sizes, np.array([1]))

    def test_n3_edge(self) -> None:
        """#5: N=3, max_b=0 → [1]."""
        sizes = _generate_block_sizes(3)
        assert np.array_equal(sizes, np.array([1]))


# ---------------------------------------------------------------------------
# P0 tests: analyze_block_average with arctan
# ---------------------------------------------------------------------------


class TestDenseSamplingArctan:
    """Arctan integration in block average analysis."""

    def test_ar1_arctan_present(self) -> None:
        """#6: AR(1) N=5000, φ=0.9 → arctan computed, asymptote reasonable."""
        phi = 0.9
        n = 5000
        series = _make_ar1(n, phi, seed=42)
        sigma = np.std(series, ddof=0)
        tau_theory = (1 + phi) / (2 * (1 - phi))
        sem_theory = sigma * np.sqrt(2 * tau_theory / n)

        result = analyze_block_average(series, dense_sampling=True)
        assert result.arctan is not None
        # Asymptote should be in the right ballpark regardless of reliability
        assert result.arctan.sem_asymptote == pytest.approx(sem_theory, rel=0.20)

    def test_legacy_pow2_still_computes_arctan(self) -> None:
        """#7: dense_sampling=False → powers-of-2 but arctan still computed."""
        series = _make_ar1(5000, 0.9, seed=42)
        result = analyze_block_average(series, dense_sampling=False)
        # block_sizes should be powers of 2
        bs = result.block_sizes
        for b in bs:
            assert b & (b - 1) == 0 or b == 0  # power-of-2 check
        # arctan should still be attempted
        assert result.arctan is not None


# ---------------------------------------------------------------------------
# P1 tests
# ---------------------------------------------------------------------------


class TestTauCorrImplied:
    """Comparison of tau_corr_implied with ACF tau."""

    def test_tau_vs_acf(self) -> None:
        """#8: tau_corr_implied vs ACF tau_corr — within 25%."""
        from md_analysis.enhanced_sampling.constrained_ti.analysis.autocorrelation import (
            analyze_autocorrelation,
        )

        # Use long series for more reliable fit
        series = _make_ar1(50000, 0.9, seed=42)
        ba_result = analyze_block_average(series, dense_sampling=True)
        acf_result = analyze_autocorrelation(series)

        assert ba_result.arctan is not None
        if ba_result.arctan.reliable:
            assert ba_result.arctan.tau_corr_implied == pytest.approx(
                acf_result.tau_corr, rel=0.25
            )


class TestConstantSeries:
    """Edge case: near-constant series."""

    def test_constant_no_crash(self) -> None:
        """#9: Constant series → arctan None or tau=nan, no crash."""
        series = np.ones(1000)
        result = analyze_block_average(series, dense_sampling=True)
        # SEM curve will be all zeros; arctan should handle gracefully
        if result.arctan is not None:
            assert np.isnan(result.arctan.tau_corr_implied) or not result.arctan.reliable
