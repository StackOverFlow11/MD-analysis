"""Tests for Flyvbjerg-Petersen block averaging."""

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.analysis.block_average import (
    _generate_block_sizes,
    analyze_block_average,
)


def make_ar1(n, phi, sigma=1.0, seed=42):
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + rng.normal(0, sigma * np.sqrt(1 - phi**2))
    return x


class TestGenerateBlockSizes:
    def test_all_powers_of_2(self):
        bs = _generate_block_sizes(2000, min_blocks=4)
        for b in bs:
            assert b & (b - 1) == 0  # power of 2

    def test_ascending(self):
        bs = _generate_block_sizes(2000, min_blocks=4)
        assert np.all(np.diff(bs) > 0)

    def test_max_respects_min_blocks(self):
        bs = _generate_block_sizes(2000, min_blocks=4)
        assert 2000 // bs[-1] >= 4

    def test_small_n(self):
        bs = _generate_block_sizes(5, min_blocks=4)
        assert bs[0] == 1
        assert len(bs) >= 1

    def test_n0_edge(self):
        bs = _generate_block_sizes(0, min_blocks=4)
        assert len(bs) == 1


class TestFPPlateau:
    def test_ar1_plateau_detected(self):
        series = make_ar1(20000, 0.9, seed=42)
        result = analyze_block_average(series)
        assert result.plateau_reached
        assert result.plateau_index is not None
        assert result.plateau_block_size is not None
        assert result.plateau_sem > 0

    def test_iid_plateau_at_small_b(self):
        rng = np.random.default_rng(42)
        series = rng.normal(0, 1, 5000)
        result = analyze_block_average(series)
        assert result.plateau_reached
        # IID: plateau should be at small B (no autocorrelation)
        assert result.plateau_block_size <= 4

    def test_delta_sem_formula(self):
        rng = np.random.default_rng(42)
        series = rng.normal(0, 1, 1000)
        result = analyze_block_average(series)
        for i, bs in enumerate(result.block_sizes):
            if np.isnan(result.sem_curve[i]):
                continue
            nb = result.n_total // bs
            expected = result.sem_curve[i] / np.sqrt(2.0 * (nb - 1))
            assert result.delta_sem[i] == pytest.approx(expected, rel=1e-10)

    def test_n_total_stored(self):
        rng = np.random.default_rng(42)
        series = rng.normal(0, 1, 1234)
        result = analyze_block_average(series)
        assert result.n_total == 1234

    def test_nb_lt_2_gives_nan(self):
        """Block sizes where n_b < 2 should produce NaN."""
        series = np.ones(7)  # N=7, B=4 → n_b=1
        result = analyze_block_average(series, min_blocks=1)
        # B=4 gives n_b=1, should be NaN
        for i, bs in enumerate(result.block_sizes):
            nb = 7 // bs
            if nb < 2:
                assert np.isnan(result.sem_curve[i])

    def test_constant_series(self):
        series = np.ones(1000)
        result = analyze_block_average(series)
        # All SEMs should be 0
        valid = ~np.isnan(result.sem_curve)
        assert np.all(result.sem_curve[valid] == 0.0)

    def test_pass_fail_with_sem_max(self):
        series = make_ar1(5000, 0.5, seed=42)
        result = analyze_block_average(series, sem_max=100.0)
        assert result.passed == True  # noqa: E712

        result2 = analyze_block_average(series, sem_max=1e-20)
        assert result2.passed == False  # noqa: E712

    def test_n_consecutive_affects_detection(self):
        series = make_ar1(5000, 0.9, seed=42)
        r1 = analyze_block_average(series, n_consecutive=1)
        r3 = analyze_block_average(series, n_consecutive=3)
        # n_consecutive=1 should detect plateau at same or earlier B
        if r1.plateau_reached and r3.plateau_reached:
            assert r1.plateau_block_size <= r3.plateau_block_size
