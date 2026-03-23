"""Step 4: Geweke stationarity test."""

from __future__ import annotations

import numpy as np

from ..config import (
    DEFAULT_ACF_ALPHA,
    DEFAULT_GEWEKE_FA,
    DEFAULT_GEWEKE_FB,
    DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES,
    DEFAULT_GEWEKE_ZCRIT,
)
from ..models import GewekeResult
from ._acf_core import compute_acf, compute_iat, compute_sem_corrected


def analyze_geweke(
    series: np.ndarray,
    *,
    f_a: float = DEFAULT_GEWEKE_FA,
    f_b: float = DEFAULT_GEWEKE_FB,
    z_crit: float = DEFAULT_GEWEKE_ZCRIT,
    alpha: int = DEFAULT_ACF_ALPHA,
    min_neff_subseries: int = DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES,
) -> GewekeResult:
    """Run Geweke stationarity test on a Lagrange multiplier time series.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Lagrange multiplier time series.
    f_a : float
        Fraction of series for front segment.
    f_b : float
        Fraction of series for rear segment.
    z_crit : float
        Critical z-value (two-sided).
    alpha : int
        ACF truncation multiplier for spectral variance estimation.
    min_neff_subseries : int
        Minimum N_eff in front segment for reliable result.

    Returns
    -------
    GewekeResult
    """
    n = len(series)
    n_a = max(int(n * f_a), 1)
    n_b = max(int(n * f_b), 1)

    seg_a = series[:n_a]
    seg_b = series[n - n_b :]

    mean_a = float(np.mean(seg_a))
    mean_b = float(np.mean(seg_b))

    # Spectral variance for each sub-series via ACF
    def _spectral_var(seg: np.ndarray) -> tuple[float, float]:
        """Return (spectral_variance, n_eff) for a sub-series."""
        nn = len(seg)
        if nn < 2 or np.var(seg) == 0.0:
            # Degenerate case: no variance
            return 0.0, float(nn)
        acf = compute_acf(seg)
        tau = compute_iat(acf, nn, alpha=alpha)
        n_eff = nn / (2.0 * tau)
        sem = compute_sem_corrected(seg, tau)
        return sem**2, n_eff

    var_a, n_eff_a_val = _spectral_var(seg_a)
    var_b, _ = _spectral_var(seg_b)

    # z statistic
    denom = np.sqrt(var_a + var_b)
    if denom > 0:
        z = (mean_a - mean_b) / denom
    else:
        z = 0.0

    # Reliability check
    reliable = n_eff_a_val >= min_neff_subseries

    return GewekeResult(
        z=float(z),
        mean_a=mean_a,
        mean_b=mean_b,
        n_eff_a=float(n_eff_a_val) if reliable else None,
        passed=abs(z) < z_crit,
        reliable=reliable,
    )
