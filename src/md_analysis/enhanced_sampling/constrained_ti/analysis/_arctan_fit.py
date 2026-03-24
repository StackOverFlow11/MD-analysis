"""Arctan extrapolation for block-average SEM curves.

Pure-function module, no dataclass definitions, no logging.
Follows the same conventions as ``_acf_core.py``.

The model ``SEM(B) = A · arctan(B · x)`` is an empirical fit — not the
exact theoretical SEM(B) curve (which depends on the ACF shape).  The
asymptotic SEM is ``A · π/2``.  See TIreconstruction.md §0.4 for details.
"""

from __future__ import annotations

import numpy as np

from ..config import DEFAULT_ARCTAN_MIN_POINTS, DEFAULT_ARCTAN_R2_MIN
from ..models import ArctanFitResult


def fit_arctan_sem(
    block_sizes: np.ndarray,
    sem_curve: np.ndarray,
    *,
    n_total: int,
    sigma_series: float,
    r2_min: float = DEFAULT_ARCTAN_R2_MIN,
    min_points: int = DEFAULT_ARCTAN_MIN_POINTS,
) -> ArctanFitResult | None:
    """Fit SEM(B) = A · arctan(B · x) via WLS and extrapolate to A · π/2.

    Parameters
    ----------
    block_sizes : np.ndarray
        Block sizes tested.
    sem_curve : np.ndarray
        SEM(B) values for each block size.
    n_total : int
        Total series length (for tau_corr_implied calculation).
    sigma_series : float
        Series standard deviation with ddof=0 (for tau_corr_implied).
    r2_min : float
        Minimum R² for the fit to be considered reliable.
    min_points : int
        Minimum number of non-NaN data points to attempt the fit.

    Returns
    -------
    ArctanFitResult or None
        None if scipy is unavailable or fewer than *min_points* valid
        data points exist.  ``ArctanFitResult(reliable=False)`` if the
        fit was attempted but failed quality checks.

    Notes
    -----
    WLS weighting: ``sigma_i = sqrt(block_sizes[i] / n_total)`` passed to
    ``scipy.optimize.curve_fit`` with ``absolute_sigma=False``.  This is a
    simplified approximation — the exact weight should also depend on
    ``SEM(B)`` — but in relative-sigma mode only the weight *ratios*
    matter, so the practical difference is small.
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        return None

    # Filter NaN values
    mask = ~np.isnan(sem_curve)
    bs_valid = block_sizes[mask].astype(float)
    sem_valid = sem_curve[mask]

    if len(bs_valid) < min_points:
        return None

    # Arctan model
    def _arctan_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
        return a * np.arctan(b * x)

    # WLS weights: sigma_i = sqrt(block_size_i / n_total) = 1/sqrt(n_blocks_i)
    sigma_weights = np.sqrt(bs_valid / n_total)

    # Initial guesses (clamp to bounds)
    a0 = max(2.0 * np.max(sem_valid) / np.pi, 1e-14)
    b0 = max(1.0 / np.median(bs_valid), 1e-14)

    # Attempt fit
    _FAIL = ArctanFitResult(
        sem_asymptote=float("nan"),
        A=float("nan"),
        B=float("nan"),
        r2=float("nan"),
        reliable=False,
        tau_corr_implied=float("nan"),
        fit_curve=None,
    )
    try:
        popt, pcov = curve_fit(
            _arctan_model,
            bs_valid,
            sem_valid,
            p0=[a0, b0],
            bounds=([1e-15, 1e-15], [np.inf, np.inf]),
            sigma=sigma_weights,
            absolute_sigma=False,
            maxfev=10000,
        )
    except (RuntimeError, ValueError):
        # curve_fit did not converge or initial guess outside bounds
        return _FAIL

    a_fit, b_fit = popt
    sem_asymptote = a_fit * np.pi / 2.0

    # Compute R²
    y_pred = _arctan_model(bs_valid, a_fit, b_fit)
    ss_res = np.sum((sem_valid - y_pred) ** 2)
    ss_tot = np.sum((sem_valid - np.mean(sem_valid)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Check pcov singularity
    pcov_singular = np.any(np.isinf(np.diag(pcov)))

    # Reliability check
    reliable = (
        r2 >= r2_min
        and a_fit > 1e-15
        and b_fit > 1e-15
        and not pcov_singular
    )

    # Compute tau_corr_implied
    if reliable and sigma_series > 1e-30:
        tau_corr_implied = n_total * sem_asymptote**2 / (2.0 * sigma_series**2)
    else:
        tau_corr_implied = float("nan")

    # Compute fit curve on ALL block sizes (including NaN positions → None there)
    fit_curve_full = np.full(len(block_sizes), np.nan)
    fit_curve_full[mask] = _arctan_model(bs_valid, a_fit, b_fit)

    return ArctanFitResult(
        sem_asymptote=float(sem_asymptote),
        A=float(a_fit),
        B=float(b_fit),
        r2=float(r2),
        reliable=reliable,
        tau_corr_implied=float(tau_corr_implied),
        fit_curve=fit_curve_full,
    )
