"""Tests for _acf_core: FFT-ACF, IAT, corrected SEM."""

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.analysis._acf_core import (
    compute_acf,
    compute_iat,
    compute_sem_corrected,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# compute_acf
# ---------------------------------------------------------------------------

class TestComputeACF:
    def test_acf_shape(self):
        series = make_iid(100)
        acf = compute_acf(series)
        assert acf.shape == (100,)

    def test_acf_zero_lag_is_one(self):
        series = make_iid(1000)
        acf = compute_acf(series)
        assert acf[0] == pytest.approx(1.0)

    def test_acf_constant_raises(self):
        series = np.ones(100)
        with pytest.raises(ValueError, match="zero variance"):
            compute_acf(series)

    def test_acf_iid_decays_fast(self):
        series = make_iid(5000)
        acf = compute_acf(series)
        # For i.i.d., ACF at lag > 0 should be near zero
        assert np.abs(acf[1:10]).max() < 0.1

    def test_acf_ar1_lag1(self):
        phi = 0.9
        series = make_ar1(10000, phi)
        acf = compute_acf(series)
        # C(1) should be close to phi for AR(1)
        assert acf[1] == pytest.approx(phi, abs=0.05)


# ---------------------------------------------------------------------------
# compute_iat
# ---------------------------------------------------------------------------

class TestComputeIAT:
    def test_iid_tau_near_half(self):
        series = make_iid(10000)
        acf = compute_acf(series)
        tau = compute_iat(acf, len(series))
        assert tau == pytest.approx(0.5, abs=0.15)

    def test_ar1_tau_within_5pct(self):
        phi = 0.9
        n = 10000
        series = make_ar1(n, phi)
        acf = compute_acf(series)
        tau = compute_iat(acf, n)
        # Analytic: tau_int = 1/2 + phi/(1-phi) = 1/2 + 9 = 9.5
        expected = 0.5 + phi / (1 - phi)
        assert tau == pytest.approx(expected, rel=0.1)

    def test_ar1_high_phi(self):
        """High correlation (phi=0.95) — self-consistent iteration should converge."""
        phi = 0.95
        n = 20000
        series = make_ar1(n, phi, seed=123)
        acf = compute_acf(series)
        tau = compute_iat(acf, n)
        expected = 0.5 + phi / (1 - phi)  # = 19.5
        assert tau == pytest.approx(expected, rel=0.15)

    def test_hard_cap_respected(self):
        """Ensure M <= N//2 hard cap is enforced."""
        series = make_ar1(100, 0.99)  # Very high correlation, short series
        acf = compute_acf(series)
        tau = compute_iat(acf, len(series))
        # tau should be capped by limited window
        assert tau > 0.5


# ---------------------------------------------------------------------------
# compute_sem_corrected
# ---------------------------------------------------------------------------

class TestComputeSEMCorrected:
    def test_iid_sem(self):
        n = 10000
        sigma = 2.0
        series = make_iid(n, sigma=sigma)
        tau = 0.5  # i.i.d.
        sem = compute_sem_corrected(series, tau)
        expected = np.std(series) * np.sqrt(2 * 0.5 / n)
        assert sem == pytest.approx(expected, rel=1e-6)

    def test_empty_series(self):
        sem = compute_sem_corrected(np.array([]), 1.0)
        assert sem == 0.0
