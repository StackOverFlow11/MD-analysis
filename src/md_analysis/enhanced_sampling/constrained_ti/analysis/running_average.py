"""Step 1: Running average drift check — cumulative-mean stability."""

from __future__ import annotations

import numpy as np

from ..config import DEFAULT_DRIFT_FACTOR
from ..models import RunningAverageResult


def analyze_running_average(
    series: np.ndarray,
    *,
    sem: float,
    drift_factor: float = DEFAULT_DRIFT_FACTOR,
) -> RunningAverageResult:
    """Check cumulative-mean stability of a Lagrange multiplier time series.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Lagrange multiplier time series.
    sem : float
        Standard error estimate (from autocorrelation or block averaging).
    drift_factor : float
        Multiplicative factor for the drift threshold.

    Returns
    -------
    RunningAverageResult
    """
    n = len(series)
    running_mean = np.cumsum(series) / np.arange(1, n + 1)

    # Drift = max - min over [N/2, N]
    half = n // 2
    second_half = running_mean[half:]
    drift_D = float(np.max(second_half) - np.min(second_half))

    drift_limit = drift_factor * sem
    passed = drift_D < drift_limit

    return RunningAverageResult(
        running_mean=running_mean,
        drift_D=drift_D,
        drift_limit=drift_limit,
        passed=passed,
    )
