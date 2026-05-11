"""Step 2: Autocorrelation analysis — tau_corr, N_eff, SEM_auto."""

from __future__ import annotations

import numpy as np

from ..config import DEFAULT_ACF_ALPHA, DEFAULT_NEFF_MIN
from ..models import AutocorrResult
from ._acf_core import compute_acf, compute_iat, compute_sem_corrected


def analyze_autocorrelation(
    series: np.ndarray,
    *,
    sem_max: float | None = None,
    alpha: int = DEFAULT_ACF_ALPHA,
    neff_min: int = DEFAULT_NEFF_MIN,
) -> AutocorrResult:
    """Run autocorrelation analysis on a Lagrange multiplier time series.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Lagrange multiplier time series (post-equilibration).
    sem_max : float or None
        Precision target. None skips the SEM pass/fail check.
    alpha : int
        Self-consistent truncation multiplier for IAT.
    neff_min : int
        Minimum effective independent samples.

    Returns
    -------
    AutocorrResult
    """
    n = len(series)

    acf = compute_acf(series)
    tau_corr = compute_iat(acf, n, alpha=alpha)
    n_eff = n / (2.0 * tau_corr)
    sem_auto = compute_sem_corrected(series, tau_corr)

    passed_neff = n_eff >= neff_min
    passed_sem = (sem_auto <= sem_max) if sem_max is not None else None

    # If N_eff insufficient, compute minimum frames needed:
    # N_min = 2 * neff_min * tau_corr
    t_min_frames: int | None = None
    if not passed_neff:
        t_min_frames = int(np.ceil(2.0 * neff_min * tau_corr))

    return AutocorrResult(
        acf=acf,
        tau_corr=tau_corr,
        n_eff=n_eff,
        sem_auto=sem_auto,
        passed_neff=passed_neff,
        passed_sem=passed_sem,
        t_min_frames=t_min_frames,
    )
