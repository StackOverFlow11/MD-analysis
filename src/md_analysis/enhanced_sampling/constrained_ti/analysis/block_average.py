"""Step 3: Flyvbjerg-Petersen block averaging — SEM(B) plateau detection.

Reference: Flyvbjerg & Petersen, J. Chem. Phys. 91, 461 (1989).
"""

from __future__ import annotations

import numpy as np

from ..config import DEFAULT_FP_CONSECUTIVE, DEFAULT_FP_MIN_BLOCKS
from ..models import BlockAverageResult


def _generate_block_sizes(n: int, min_blocks: int = DEFAULT_FP_MIN_BLOCKS) -> np.ndarray:
    """Powers-of-2 block sizes up to n // min_blocks."""
    max_b = n // min_blocks if min_blocks > 0 else n
    block_sizes = []
    b = 1
    while b <= max_b:
        block_sizes.append(b)
        b *= 2
    return np.array(block_sizes, dtype=int) if block_sizes else np.array([1], dtype=int)


def analyze_block_average(
    series: np.ndarray,
    *,
    sem_max: float | None = None,
    min_blocks: int = DEFAULT_FP_MIN_BLOCKS,
    n_consecutive: int = DEFAULT_FP_CONSECUTIVE,
) -> BlockAverageResult:
    """Run Flyvbjerg-Petersen block averaging on a time series.

    Parameters
    ----------
    series : np.ndarray, shape (N,)
        Time series (e.g. Lagrange multiplier).
    sem_max : float or None
        Precision target.  None skips the pass/fail check.
    min_blocks : int
        Minimum number of blocks at largest B.
    n_consecutive : int
        Consecutive pow2 levels with increase < δSEM to declare plateau.

    Returns
    -------
    BlockAverageResult
    """
    n = len(series)
    block_sizes = _generate_block_sizes(n, min_blocks)

    sem_curve = np.empty(len(block_sizes))
    delta_sem = np.empty(len(block_sizes))

    for i, bs in enumerate(block_sizes):
        nb = n // bs
        if nb < 2:
            sem_curve[i] = np.nan
            delta_sem[i] = np.nan
            continue
        block_means = series[: nb * bs].reshape(nb, bs).mean(axis=1)
        sem_curve[i] = np.std(block_means, ddof=1) / np.sqrt(nb)
        delta_sem[i] = sem_curve[i] / np.sqrt(2.0 * (nb - 1))

    # Plateau detection: consecutive levels where increase < δSEM
    plateau_index: int | None = None
    consec = 0
    for i in range(len(block_sizes) - 1):
        if np.isnan(sem_curve[i]) or np.isnan(sem_curve[i + 1]):
            consec = 0
            continue
        increase = sem_curve[i + 1] - sem_curve[i]
        if increase < delta_sem[i]:
            consec += 1
            if consec >= n_consecutive and plateau_index is None:
                # plateau starts at the first point of the qualifying run
                plateau_index = i - n_consecutive + 2
        else:
            consec = 0

    plateau_reached = plateau_index is not None

    if plateau_reached:
        plateau_sem = float(sem_curve[plateau_index])
        plateau_delta = float(delta_sem[plateau_index])
        plateau_block_size: int | None = int(block_sizes[plateau_index])
    else:
        # Fallback: use the last valid point
        valid = np.where(~np.isnan(sem_curve))[0]
        last = valid[-1] if len(valid) > 0 else 0
        plateau_sem = float(sem_curve[last])
        plateau_delta = float(delta_sem[last])
        plateau_block_size = None

    # Pass/fail
    if sem_max is not None:
        passed: bool | None = plateau_reached and plateau_sem <= sem_max
    else:
        passed = None

    return BlockAverageResult(
        block_sizes=block_sizes,
        sem_curve=sem_curve,
        delta_sem=delta_sem,
        n_total=n,
        plateau_index=plateau_index,
        plateau_sem=plateau_sem,
        plateau_delta=plateau_delta,
        plateau_block_size=plateau_block_size,
        plateau_reached=plateau_reached,
        passed=passed,
    )
