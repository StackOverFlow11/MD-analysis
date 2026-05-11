"""Shared ACF computation used by autocorrelation and Geweke engines.

Public (internal to analysis package)
--------------------------------------
- ``compute_acf``  — FFT-accelerated normalized ACF
- ``compute_iat``  — self-consistent truncation IAT (Sokal 1997)
- ``compute_sem_corrected`` — sigma * sqrt(2 * tau / N)

tau_corr definition (Sokal 1997 convention)
-------------------------------------------
tau_int = 1/2 + sum_{j=1}^{M} C(j), where C(j) is the normalized ACF.
N_eff = N / (2 * tau_int).
"""

from __future__ import annotations

import numpy as np

from ..config import DEFAULT_ACF_ALPHA, DEFAULT_ACF_M_MAX_FRACTION


def compute_acf(series: np.ndarray) -> np.ndarray:
    """FFT-accelerated normalized autocorrelation function.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Input time series.

    Returns
    -------
    np.ndarray, shape (N,)
        C(j) for j = 0, ..., N-1 with C(0) = 1.

    Raises
    ------
    ValueError
        If series has zero variance (constant input).
    """
    n = len(series)
    if n < 2:
        return np.ones(1)

    mean = np.mean(series)
    var = np.var(series)
    if var == 0.0:
        raise ValueError(
            "Cannot compute ACF for a constant series (zero variance)."
        )

    # Zero-pad to next power of 2 for FFT efficiency
    x = series - mean
    fft_size = 1
    while fft_size < 2 * n:
        fft_size *= 2

    fft_x = np.fft.rfft(x, n=fft_size)
    acf_raw = np.fft.irfft(fft_x * np.conj(fft_x), n=fft_size)[:n]

    # Normalize: C(0) = 1
    acf_raw /= acf_raw[0]
    return acf_raw


def compute_iat(
    acf: np.ndarray,
    n: int,
    *,
    alpha: int = DEFAULT_ACF_ALPHA,
    m_max_fraction: float = DEFAULT_ACF_M_MAX_FRACTION,
) -> float:
    """Self-consistent truncation integrated autocorrelation time.

    Iteratively sums C(j) until M >= alpha * tau_corr, with a hard upper
    limit M <= n * m_max_fraction to prevent noise-dominated estimates.

    Parameters
    ----------
    acf : np.ndarray
        Normalized ACF from ``compute_acf``.
    n : int
        Length of original series (for hard cap computation).
    alpha : int
        Self-consistent truncation multiplier.
    m_max_fraction : float
        Hard upper limit as fraction of N.

    Returns
    -------
    float
        tau_int (in frames), including the 1/2 base term.
    """
    m_max = int(n * m_max_fraction)
    m_max = min(m_max, len(acf) - 1)

    # Start with tau = 0.5 (the base term)
    tau = 0.5
    for j in range(1, m_max + 1):
        tau += acf[j]
        # Self-consistent check: stop when window is large enough
        if j >= alpha * tau:
            break

    # Ensure tau >= 0.5 (minimum for uncorrelated data)
    return max(tau, 0.5)


def compute_sem_corrected(
    series: np.ndarray,
    tau_corr: float,
) -> float:
    """Standard error of the mean corrected for autocorrelation.

    SEM = sigma_lambda * sqrt(2 * tau_corr / N)

    Parameters
    ----------
    series : np.ndarray
        Input time series.
    tau_corr : float
        Integrated autocorrelation time (frames).

    Returns
    -------
    float
        Corrected SEM.
    """
    n = len(series)
    if n == 0:
        return 0.0
    sigma = np.std(series, ddof=0)
    return sigma * np.sqrt(2.0 * tau_corr / n)
