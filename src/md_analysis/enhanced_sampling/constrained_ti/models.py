"""Frozen dataclasses and exception hierarchy for constrained-TI analysis.

All dataclasses are ``frozen=True``, matching the ``Slowgrowth`` convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ...exceptions import MDAnalysisError

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class ConvergenceError(MDAnalysisError):
    """Base for constrained-TI convergence errors."""


class InsufficientSamplingError(ConvergenceError):
    """Raised when sampling is too short for any reliable estimate."""


# ---------------------------------------------------------------------------
# Per-engine result dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RunningAverageResult:
    """Output of Step 1: cumulative-mean drift analysis."""

    running_mean: np.ndarray  # shape (N,) — cumulative mean curve
    drift_D: float  # max - min of running_mean over [N/2, N]
    drift_limit: float  # drift_factor * SEM (the comparison target)
    passed: bool


@dataclass(frozen=True)
class AutocorrResult:
    """Output of Step 2: autocorrelation analysis.

    ``tau_corr`` uses the Sokal (1997) convention:
    tau_int = 1/2 + sum_{j=1}^{M} C(j), so N_eff = N / (2 * tau_int).
    """

    acf: np.ndarray  # C(j), truncated at self-consistent cutoff M
    tau_corr: float  # integrated autocorrelation time (frames)
    n_eff: float  # effective independent sample count
    sem_auto: float  # sigma_lambda * sqrt(2 * tau / N)
    passed_neff: bool  # N_eff >= N_EFF_MIN
    passed_sem: bool | None  # SEM_auto <= SEM_max; None if sem_max not provided
    t_min_frames: int | None  # minimum frames needed if N_eff < threshold; else None


@dataclass(frozen=True)
class BlockAverageResult:
    """Output of Step 3: Flyvbjerg-Petersen block averaging.

    Block sizes are powers of 2.  Plateau is detected when the SEM
    increase between consecutive levels falls below δSEM for
    *n_consecutive* levels.  δSEM(B) = SEM(B) / √(2(n_b − 1)).
    """

    block_sizes: np.ndarray  # pow2 block sizes tested
    sem_curve: np.ndarray  # SEM(B) at each level
    delta_sem: np.ndarray  # δSEM(B) = SEM(B) / √(2(n_b − 1))
    n_total: int  # series length (n_b = n_total // B)
    plateau_index: int | None  # index where plateau first detected; None if not
    plateau_sem: float  # SEM at plateau (or SEM at max B if no plateau)
    plateau_delta: float  # δSEM at plateau point
    plateau_block_size: int | None  # B at plateau; None if not detected
    plateau_reached: bool  # True if plateau was detected
    passed: bool | None  # plateau_reached AND SEM <= SEM_max; None if no target


@dataclass(frozen=True)
class GewekeResult:
    """Output of Step 4: Geweke stationarity test."""

    z: float  # test statistic
    mean_a: float  # front-segment mean
    mean_b: float  # rear-segment mean
    n_eff_a: float | None  # effective samples in front segment; None if unreliable
    passed: bool  # |z| < z_crit
    reliable: bool  # False if front sub-series too short for spectral variance


# ---------------------------------------------------------------------------
# Per-point composite structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstraintPointInput:
    """Everything needed to diagnose one constraint point.

    Supports two modes:

    - **TI context** (multi-point): ``weight`` and ``sem_max`` are computed
      from trapezoid quadrature and epsilon_tol; pass/fail is determined.
    - **Standalone** (single-point): ``weight`` and ``sem_max`` are None;
      all four diagnostics run but SEM pass/fail is skipped.
    """

    xi: float  # CV constraint value (or user label)
    lambda_series: np.ndarray  # shape (N,) — post-equilibration
    dt: float  # frame interval (fs)
    time_start_fs: float  # absolute start time of analyzed window (fs)
    weight: float | None  # trapezoid weight w_k (None in standalone)
    sem_max: float | None  # precision target (None in standalone)
    point_index: int | None  # 0-based index; None in standalone


@dataclass(frozen=True)
class ConstraintPointReport:
    """Complete diagnostic output for one constraint point."""

    xi: float
    point_index: int | None
    n_analyzed: int  # number of frames used in analysis (post-equilibration)
    time_start_fs: float  # start time of analyzed window (fs)
    time_end_fs: float  # end time of analyzed window (fs)
    lambda_mean: float  # mean constraint force at this point
    sigma_lambda: float  # sample standard deviation of lambda

    # Engine results (carried verbatim for downstream inspection)
    autocorr: AutocorrResult
    block_avg: BlockAverageResult
    running_avg: RunningAverageResult
    geweke: GewekeResult

    # Summary
    sem_final: float  # SEM_block (primary) or SEM_auto (fallback)
    sem_max: float | None  # precision target; None in standalone
    passed: bool | None  # all four passed; None if sem_max not set
    failure_reasons: tuple[str, ...]  # empty if passed


# ---------------------------------------------------------------------------
# Multi-point report
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TIReport:
    """Full thermodynamic-integration convergence report."""

    point_reports: tuple[ConstraintPointReport, ...]
    xi_values: np.ndarray  # shape (K,)
    weights: np.ndarray  # shape (K,) — trapezoid weights
    forces: np.ndarray  # shape (K,) — lambda_mean at each point
    force_errors: np.ndarray  # shape (K,) — SEM_final at each point

    delta_A: float  # integrated free-energy difference
    sigma_A: float  # propagated statistical error
    epsilon_tol_au: float  # tolerance in Hartree (internal unit)

    all_passed: bool  # every point passed
    failing_indices: tuple[int, ...]  # 0-based indices into point_reports

    # Optional: optimal allocation suggestion
    suggested_time_ratios: np.ndarray | None  # shape (K,) or None


# ---------------------------------------------------------------------------
# I/O data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TIPointDefinition:
    """Discovered constraint point file paths."""

    xi: float  # CV constraint value parsed from directory name
    restart_path: Path  # path to *.restart file
    log_path: Path  # path to *.LagrangeMultLog file
