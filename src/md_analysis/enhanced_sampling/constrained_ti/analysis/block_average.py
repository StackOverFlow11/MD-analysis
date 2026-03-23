"""Step 3: Block averaging (Flyvbjerg-Petersen) — SEM(B) plateau detection."""

from __future__ import annotations

import numpy as np

from ..config import (
    DEFAULT_CROSS_VALID_RTOL,
    DEFAULT_MIN_BLOCKS,
    DEFAULT_PLATEAU_RTOL,
    DEFAULT_PLATEAU_WINDOW,
)
from ..models import BlockAverageResult


def analyze_block_average(
    series: np.ndarray,
    *,
    sem_max: float | None = None,
    sem_auto: float | None = None,
    min_blocks: int = DEFAULT_MIN_BLOCKS,
    plateau_window: int = DEFAULT_PLATEAU_WINDOW,
    plateau_rtol: float = DEFAULT_PLATEAU_RTOL,
    cross_valid_rtol: float = DEFAULT_CROSS_VALID_RTOL,
) -> BlockAverageResult:
    """Run block-averaging analysis on a Lagrange multiplier time series.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Lagrange multiplier time series.
    sem_max : float or None
        Precision target. None skips the SEM pass/fail check.
    sem_auto : float or None
        SEM from autocorrelation analysis, for cross-validation.
    min_blocks : int
        Minimum number of blocks at largest B.
    plateau_window : int
        Number of largest-B points used in plateau check.
    plateau_rtol : float
        Maximum relative spread for plateau detection.
    cross_valid_rtol : float
        Cross-validation tolerance between SEM_block and SEM_auto.

    Returns
    -------
    BlockAverageResult
    """
    n = len(series)

    # Generate block sizes B = 1, 2, 4, 8, ..., up to N // min_blocks
    max_b = n // min_blocks if min_blocks > 0 else n
    block_sizes = []
    b = 1
    while b <= max_b:
        block_sizes.append(b)
        b *= 2
    block_sizes = np.array(block_sizes, dtype=int)

    # Compute SEM(B) for each block size
    sem_curve = np.empty(len(block_sizes))
    for i, bs in enumerate(block_sizes):
        n_blocks = n // bs
        if n_blocks < 1:
            sem_curve[i] = np.nan
            continue
        # Compute block means
        trimmed = series[: n_blocks * bs].reshape(n_blocks, bs)
        block_means = trimmed.mean(axis=1)
        # SEM of block means
        sem_curve[i] = np.std(block_means, ddof=1) / np.sqrt(n_blocks)

    # SEM at largest B (for fallback transparency)
    sem_at_max_B = float(sem_curve[-1]) if len(sem_curve) > 0 else 0.0

    # Plateau detection: last `plateau_window` valid entries
    valid_mask = ~np.isnan(sem_curve)
    valid_indices = np.where(valid_mask)[0]

    if len(valid_indices) >= plateau_window:
        plateau_indices = valid_indices[-plateau_window:]
        plateau_values = sem_curve[plateau_indices]
        s_bar = np.mean(plateau_values)
        delta_s = np.max(plateau_values) - np.min(plateau_values)
        rtol = delta_s / s_bar if s_bar > 0 else np.inf
        plateau_reached = rtol < plateau_rtol
        sem_plateau_val = float(s_bar)
    else:
        # Not enough points for plateau check
        sem_plateau_val = float(sem_at_max_B)
        rtol = np.inf
        plateau_reached = False

    # Cross-validation with SEM_auto
    cross_valid_ok = True
    if sem_auto is not None and sem_plateau_val > 0:
        relative_diff = abs(sem_plateau_val - sem_auto) / sem_plateau_val
        cross_valid_ok = relative_diff < cross_valid_rtol

    # Pass/fail
    if sem_max is not None:
        passed: bool | None = plateau_reached and sem_plateau_val <= sem_max
    else:
        passed = None

    return BlockAverageResult(
        block_sizes=block_sizes,
        sem_curve=sem_curve,
        sem_plateau=sem_plateau_val,
        sem_at_max_B=sem_at_max_B,
        plateau_rtol=float(rtol),
        plateau_reached=plateau_reached,
        cross_valid_ok=cross_valid_ok,
        passed=passed,
    )
