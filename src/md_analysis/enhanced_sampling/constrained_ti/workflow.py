"""Orchestrator for constrained-TI convergence diagnostics.

Hosts:
- ``analyze_single_point``    — single point in TI context
- ``analyze_standalone``      — standalone single-point (no TI context)
- ``analyze_ti``              — full multi-point TI analysis
- ``standalone_diagnostics``  — unified entry: parse + analyze + plot + CSV
- ``run_ti_full_from_root``   — agent-facing end-to-end wrapper (discover +
                                 load + analyze + plot + CSV + metrics)
- ``TIFullAnalysisReport``    — JSON-serializable result of the wrapper
- CSV export helpers
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .analysis.autocorrelation import analyze_autocorrelation
from .analysis.block_average import analyze_block_average
from .analysis.geweke import analyze_geweke
from .analysis.running_average import analyze_running_average
from .config import (
    DEFAULT_AUTO_EQUIL_MIN_FRAMES,
    DEFAULT_CROSS_CHECK_RTOL,
    DEFAULT_EPSILON_TOL_EV,
    DEFAULT_N_MIN,
    DEFAULT_N_WARN_SHORT,
    DEFAULT_N_WARN_UNRELIABLE,
    DEFAULT_NAN_FRACTION_MAX,
    DEFAULT_STANDALONE_CSV_NAME,
    DEFAULT_STANDALONE_PNG_NAME,
    EV_TO_HARTREE,
)
from ...utils.constants import HA_TO_EV
from .integration import (
    _compute_sem_targets,
    _integrate_free_energy,
    _suggest_time_allocation,
    compute_trapezoid_weights,
)
from .models import (
    ConstraintPointInput,
    ConstraintPointReport,
    InsufficientSamplingError,
    TIReport,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pre-analysis validation
# ---------------------------------------------------------------------------


def _validate_and_trim(
    series: np.ndarray,
    equilibration: int = 0,
    nan_fraction_max: float = DEFAULT_NAN_FRACTION_MAX,
) -> tuple[np.ndarray, list[str]]:
    """Validate series, trim equilibration, handle NaN.

    Returns
    -------
    tuple[np.ndarray, list[str]]
        Cleaned series and list of warning messages.
    """
    warnings: list[str] = []

    # Equilibration trim
    if equilibration > 0:
        if equilibration >= len(series):
            raise InsufficientSamplingError(
                f"Equilibration ({equilibration}) >= series length ({len(series)}). "
                "Nothing left after trimming."
            )
        series = series[equilibration:]

    # NaN handling
    nan_mask = np.isnan(series)
    nan_count = int(np.sum(nan_mask))
    if nan_count > 0:
        nan_frac = nan_count / len(series)
        if nan_frac > nan_fraction_max:
            raise InsufficientSamplingError(
                f"NaN fraction ({nan_frac:.1%}) exceeds maximum "
                f"({nan_fraction_max:.0%}). Series too corrupted for analysis."
            )
        # Truncate to last contiguous non-NaN segment
        last_valid = len(series) - 1
        while last_valid >= 0 and np.isnan(series[last_valid]):
            last_valid -= 1
        if last_valid < 0:
            raise InsufficientSamplingError("Series is entirely NaN.")
        # Find the start of the last contiguous non-NaN block
        first_valid = 0
        while first_valid <= last_valid and np.isnan(series[first_valid]):
            first_valid += 1
        series = series[first_valid : last_valid + 1]
        # Remove any remaining interior NaN by forward-fill (minimal impact)
        interior_nan = np.isnan(series)
        if np.any(interior_nan):
            # Use last valid value for fill
            for i in range(len(series)):
                if np.isnan(series[i]) and i > 0:
                    series[i] = series[i - 1]
        warnings.append(
            f"Removed {nan_count} NaN frames; using {len(series)} valid frames."
        )

    # Length checks
    n = len(series)
    if n < DEFAULT_N_MIN:
        raise InsufficientSamplingError(
            f"Series length ({n}) < minimum ({DEFAULT_N_MIN}) after trimming."
        )
    if n < DEFAULT_N_WARN_UNRELIABLE:
        warnings.append(
            f"Extremely short series (N={n}) — ACF/block results may be unreliable."
        )
    elif n < DEFAULT_N_WARN_SHORT:
        warnings.append(
            f"Short series (N={n}) — results should be interpreted with caution."
        )

    return series, warnings


# ---------------------------------------------------------------------------
# Engine override dispatch
# ---------------------------------------------------------------------------

_ACF_KEYS = {"alpha", "neff_min"}
_BLOCK_KEYS = {"min_blocks", "n_consecutive"}
_RUNNING_KEYS = {"drift_factor"}
_GEWEKE_KEYS = {"f_a", "f_b", "z_crit", "alpha", "min_neff_subseries"}


def _dispatch_overrides(overrides: dict) -> dict:
    """Split engine overrides into per-engine dicts."""
    acf_kw = {k: v for k, v in overrides.items() if k in _ACF_KEYS}
    block_kw = {k: v for k, v in overrides.items() if k in _BLOCK_KEYS}
    running_kw = {k: v for k, v in overrides.items() if k in _RUNNING_KEYS}
    geweke_kw = {k: v for k, v in overrides.items() if k in _GEWEKE_KEYS}
    return {
        "acf": acf_kw,
        "block": block_kw,
        "running": running_kw,
        "geweke": geweke_kw,
    }


# ---------------------------------------------------------------------------
# Single-point analysis
# ---------------------------------------------------------------------------


def analyze_single_point(
    inp: ConstraintPointInput,
    **overrides: object,
) -> ConstraintPointReport:
    """Run all four diagnostics on one constraint point.

    Parameters
    ----------
    inp : ConstraintPointInput
    **overrides
        Engine-specific parameter overrides.

    Returns
    -------
    ConstraintPointReport
    """
    series = inp.lambda_series
    dispatch = _dispatch_overrides(overrides)
    failure_reasons: list[str] = []

    # Step 1 (executed as step 2 internally): Autocorrelation
    autocorr = analyze_autocorrelation(
        series, sem_max=inp.sem_max, **dispatch["acf"]
    )

    # Step 2 (executed as step 3): Block averaging (F&P)
    block_avg = analyze_block_average(
        series,
        sem_max=inp.sem_max,
        **dispatch["block"],
    )

    # Determine sem_final: always prefer block-average SEM.
    # When plateau is detected, plateau_sem is from the plateau start;
    # otherwise block_average.py falls back to the largest valid block size.
    sem_final = block_avg.plateau_sem

    if not block_avg.plateau_reached:
        failure_reasons.append(
            "F&P block-average plateau not reached; using largest-block SEM."
        )

    # Cross-check: SEM_block vs SEM_auto
    _max_sem = max(block_avg.plateau_sem, autocorr.sem_auto)
    if _max_sem > 0:
        _rel_diff = abs(block_avg.plateau_sem - autocorr.sem_auto) / _max_sem
        if _rel_diff > DEFAULT_CROSS_CHECK_RTOL:
            failure_reasons.append(
                f"SEM_block ({block_avg.plateau_sem:.2e}) and SEM_auto "
                f"({autocorr.sem_auto:.2e}) disagree by {_rel_diff:.0%}."
            )

    # Step 3 (executed as step 4): Running average
    running_avg = analyze_running_average(
        series, sem=sem_final, **dispatch["running"]
    )

    # Step 4 (executed last): Geweke
    geweke = analyze_geweke(series, **dispatch["geweke"])

    # Aggregate pass/fail
    lambda_mean = float(np.mean(series))
    sigma_lambda = float(np.std(series, ddof=0))

    # Check tau_corr vs N/10 warning
    n = len(series)
    if autocorr.tau_corr > n / 10:
        failure_reasons.append(
            f"tau_corr ({autocorr.tau_corr:.1f}) > N/10 ({n/10:.0f}); "
            "IAT estimate may be unreliable."
        )

    if not autocorr.passed_neff:
        failure_reasons.append(
            f"N_eff ({autocorr.n_eff:.1f}) < minimum; "
            f"need ~{autocorr.t_min_frames} frames."
        )

    if autocorr.passed_sem is False:
        failure_reasons.append(
            f"SEM_auto ({autocorr.sem_auto:.6f}) > SEM_max ({inp.sem_max})."
        )

    if block_avg.passed is False and block_avg.plateau_reached:
        failure_reasons.append(
            f"SEM_block ({block_avg.plateau_sem:.6f}) > SEM_max ({inp.sem_max})."
        )

    if not running_avg.passed:
        failure_reasons.append(
            f"Running-average drift D={running_avg.drift_D:.6f} >= "
            f"limit={running_avg.drift_limit:.6f}."
        )

    if not geweke.passed:
        failure_reasons.append(
            f"Geweke |z|={abs(geweke.z):.3f} >= 1.96; "
            "series may be non-stationary."
        )
    if not geweke.reliable:
        failure_reasons.append(
            "Geweke front sub-series too short for reliable spectral variance."
        )

    # Overall passed (use bool() to handle np.bool_ safely)
    # The overall pass/fail is based on:
    #   1. N_eff >= threshold
    #   2. sem_final <= sem_max (the already-selected best SEM estimate)
    #   3. Running average drift check
    #   4. Geweke stationarity
    # Note: block_avg.passed and autocorr.passed_sem are per-engine
    # diagnostics; the overall SEM check uses sem_final directly.
    if inp.sem_max is not None:
        sem_ok = sem_final <= inp.sem_max
        all_engine_pass = (
            bool(autocorr.passed_neff)
            and sem_ok
            and bool(running_avg.passed)
            and bool(geweke.passed)
        )
        passed: bool | None = all_engine_pass
    else:
        # Standalone without SEM target: can still fail on non-SEM criteria
        passed = None

    n = len(inp.lambda_series)
    return ConstraintPointReport(
        xi=inp.xi,
        point_index=inp.point_index,
        n_analyzed=n,
        time_start_fs=inp.time_start_fs,
        time_end_fs=inp.time_start_fs + (n - 1) * inp.dt,
        lambda_mean=lambda_mean,
        sigma_lambda=sigma_lambda,
        autocorr=autocorr,
        block_avg=block_avg,
        running_avg=running_avg,
        geweke=geweke,
        sem_final=sem_final,
        sem_max=inp.sem_max,
        passed=passed,
        failure_reasons=tuple(failure_reasons),
    )


# ---------------------------------------------------------------------------
# Auto-equilibration (binary halving)
# ---------------------------------------------------------------------------


def _is_converged(report: ConstraintPointReport) -> bool:
    """Check if a report indicates convergence.

    For TI context (sem_max set): uses the composite ``passed`` flag.
    For standalone (sem_max is None): checks Geweke, running-avg, and N_eff.
    """
    if report.passed is not None:
        return bool(report.passed)
    return (
        bool(report.geweke.passed)
        and bool(report.running_avg.passed)
        and bool(report.autocorr.passed_neff)
    )


def _auto_equilibrate(
    series: np.ndarray,
    *,
    dt: float,
    xi: float,
    time_start_fs: float,
    weight: float | None,
    sem_max: float | None,
    point_index: int | None,
    min_frames: int = DEFAULT_AUTO_EQUIL_MIN_FRAMES,
    **engine_overrides: object,
) -> ConstraintPointReport:
    """Iteratively discard the first half until diagnostics pass.

    At each iteration the series is halved (keep second half).  Stops when
    the diagnostics pass or the remaining series is shorter than
    *min_frames*.

    The returned report's ``failure_reasons`` records how many frames were
    used and how many iterations were performed.
    """
    n_original = len(series)
    remaining = series
    offset = 0
    n_iter = 0

    while True:
        inp = ConstraintPointInput(
            xi=xi,
            lambda_series=remaining,
            dt=dt,
            time_start_fs=time_start_fs + offset * dt,
            weight=weight,
            sem_max=sem_max,
            point_index=point_index,
        )
        report = analyze_single_point(inp, **engine_overrides)
        n_iter += 1

        if _is_converged(report):
            break

        half = len(remaining) // 2
        if half < min_frames:
            break

        offset += half
        remaining = remaining[half:]

    # Annotate the report with auto-equilibration info
    n_used = len(remaining)
    pct = n_used / n_original * 100
    info = (
        f"Auto-equilibration: used last {n_used}/{n_original} frames "
        f"({pct:.0f}%), {n_iter} iteration(s)."
    )
    return ConstraintPointReport(
        xi=report.xi,
        point_index=report.point_index,
        n_analyzed=report.n_analyzed,
        time_start_fs=report.time_start_fs,
        time_end_fs=report.time_end_fs,
        lambda_mean=report.lambda_mean,
        sigma_lambda=report.sigma_lambda,
        autocorr=report.autocorr,
        block_avg=report.block_avg,
        running_avg=report.running_avg,
        geweke=report.geweke,
        sem_final=report.sem_final,
        sem_max=report.sem_max,
        passed=report.passed,
        failure_reasons=(info,) + report.failure_reasons,
    )


# ---------------------------------------------------------------------------
# Standalone single-point
# ---------------------------------------------------------------------------


def analyze_standalone(
    lambda_series: np.ndarray,
    *,
    dt: float = 1.0,
    xi: float = 0.0,
    sem_target: float | None = None,
    equilibration: int = 0,
    time_start_fs: float = 0.0,
    auto_equilibration: bool = False,
    **engine_overrides: object,
) -> ConstraintPointReport:
    """Diagnose one constraint point independently (no TI context).

    Parameters
    ----------
    lambda_series : np.ndarray, shape (N,)
    dt : float
        Frame interval in fs.
    xi : float
        CV label for reporting.
    sem_target : float or None
        Optional precision target (same unit as lambda).
    equilibration : int
        Frames to discard from start.
    auto_equilibration : bool
        If True, iteratively halve the series (keep second half) until
        convergence diagnostics pass or data is exhausted.
    **engine_overrides
        Forwarded to analysis engines.

    Returns
    -------
    ConstraintPointReport
    """
    series = lambda_series.copy()
    series, pre_warnings = _validate_and_trim(series, equilibration=equilibration)

    t_start = time_start_fs + equilibration * dt

    if auto_equilibration:
        report = _auto_equilibrate(
            series,
            dt=dt,
            xi=xi,
            time_start_fs=t_start,
            weight=None,
            sem_max=sem_target,
            point_index=None,
            **engine_overrides,
        )
    else:
        inp = ConstraintPointInput(
            xi=xi,
            lambda_series=series,
            dt=dt,
            time_start_fs=t_start,
            weight=None,
            sem_max=sem_target,
            point_index=None,
        )
        report = analyze_single_point(inp, **engine_overrides)

    # Prepend pre-analysis warnings
    if pre_warnings:
        report = ConstraintPointReport(
            xi=report.xi,
            point_index=report.point_index,
            n_analyzed=report.n_analyzed,
            time_start_fs=report.time_start_fs,
            time_end_fs=report.time_end_fs,
            lambda_mean=report.lambda_mean,
            sigma_lambda=report.sigma_lambda,
            autocorr=report.autocorr,
            block_avg=report.block_avg,
            running_avg=report.running_avg,
            geweke=report.geweke,
            sem_final=report.sem_final,
            sem_max=report.sem_max,
            passed=report.passed,
            failure_reasons=tuple(pre_warnings) + report.failure_reasons,
        )
    return report


# ---------------------------------------------------------------------------
# Full TI analysis
# ---------------------------------------------------------------------------


def analyze_ti(
    xi_values: np.ndarray,
    lambda_series_list: list[np.ndarray],
    dt: float,
    *,
    epsilon_tol_ev: float = DEFAULT_EPSILON_TOL_EV,
    equilibration: int | list[int] = 0,
    time_starts: list[float] | None = None,
    auto_equilibration: bool = False,
    **engine_overrides: object,
) -> TIReport:
    """Run full constrained-TI convergence analysis.

    Parameters
    ----------
    xi_values : np.ndarray, shape (K,)
        Constraint point CV values.
    lambda_series_list : list[np.ndarray]
        K arrays, each shape (N_k,).
    dt : float
        Frame interval in fs.
    epsilon_tol_ev : float
        Free-energy tolerance in eV.
    equilibration : int or list[int]
        Frames to discard per point.
    time_starts : list[float] or None
        Absolute start time (fs) per point from restart files.
        If None, defaults to 0.0 for all points.
    auto_equilibration : bool
        If True, iteratively halve each point's series until convergence
        diagnostics pass or data is exhausted.
    **engine_overrides
        Forwarded to analysis engines.

    Returns
    -------
    TIReport
    """
    k = len(xi_values)
    epsilon_tol_au = epsilon_tol_ev * EV_TO_HARTREE

    # Compute weights and SEM targets
    weights = compute_trapezoid_weights(xi_values)
    sem_targets = _compute_sem_targets(weights, epsilon_tol_au)

    # Normalize equilibration
    if isinstance(equilibration, int):
        equil_list = [equilibration] * k
    else:
        equil_list = list(equilibration)
        if len(equil_list) != k:
            raise ValueError(
                f"equilibration list length ({len(equil_list)}) != K ({k})."
            )

    # Normalize time_starts
    if time_starts is None:
        ts_list = [0.0] * k
    else:
        ts_list = list(time_starts)

    # Analyze each point
    reports: list[ConstraintPointReport] = []
    for i in range(k):
        series = lambda_series_list[i].copy()
        series, pre_warnings = _validate_and_trim(
            series, equilibration=equil_list[i]
        )

        # Analyzed window starts after equilibration
        t_start_analyzed = ts_list[i] + equil_list[i] * dt

        if auto_equilibration:
            report = _auto_equilibrate(
                series,
                dt=dt,
                xi=float(xi_values[i]),
                time_start_fs=t_start_analyzed,
                weight=float(weights[i]),
                sem_max=float(sem_targets[i]),
                point_index=i,
                **engine_overrides,
            )
        else:
            inp = ConstraintPointInput(
                xi=float(xi_values[i]),
                lambda_series=series,
                dt=dt,
                time_start_fs=t_start_analyzed,
                weight=float(weights[i]),
                sem_max=float(sem_targets[i]),
                point_index=i,
            )
            report = analyze_single_point(inp, **engine_overrides)
        if pre_warnings:
            report = ConstraintPointReport(
                xi=report.xi,
                point_index=report.point_index,
                n_analyzed=report.n_analyzed,
                time_start_fs=report.time_start_fs,
                time_end_fs=report.time_end_fs,
                lambda_mean=report.lambda_mean,
                sigma_lambda=report.sigma_lambda,
                autocorr=report.autocorr,
                block_avg=report.block_avg,
                running_avg=report.running_avg,
                geweke=report.geweke,
                sem_final=report.sem_final,
                sem_max=report.sem_max,
                passed=report.passed,
                failure_reasons=tuple(pre_warnings) + report.failure_reasons,
            )
        reports.append(report)

    # Integrate
    # dA/dξ = -⟨λ_shake⟩  (standard constrained-MD / Blue Moon formula)
    forces = -np.array([r.lambda_mean for r in reports])
    force_errors = np.array([r.sem_final for r in reports])
    delta_A, sigma_A = _integrate_free_energy(forces, weights, force_errors)

    # Failing indices (use == to handle np.bool_ types)
    failing = tuple(
        i for i, r in enumerate(reports) if r.passed is not None and not r.passed
    )

    # Suggest time allocation if any failed
    suggested = None
    if failing:
        sigmas = np.array([r.sigma_lambda for r in reports])
        tau_corrs = np.array([r.autocorr.tau_corr for r in reports])
        suggested = _suggest_time_allocation(weights, sigmas, tau_corrs)

    return TIReport(
        point_reports=tuple(reports),
        xi_values=xi_values.copy(),
        weights=weights,
        forces=forces,
        force_errors=force_errors,
        delta_A=float(delta_A),
        sigma_A=float(sigma_A),
        epsilon_tol_au=epsilon_tol_au,
        all_passed=len(failing) == 0,
        failing_indices=failing,
        suggested_time_ratios=suggested,
    )


# ---------------------------------------------------------------------------
# Unified standalone entry point
# ---------------------------------------------------------------------------


def standalone_diagnostics(
    restart_path: str,
    log_path: str,
    *,
    equilibration: int = 0,
    sem_target: float | None = None,
    colvar_id: int | None = None,
    output_dir: Path | None = None,
    auto_equilibration: bool = False,
) -> dict[str, Path | ConstraintPointReport]:
    """Parse + analyze + plot + CSV for one constraint point.

    Parameters
    ----------
    restart_path : str
        Path to .restart file.
    log_path : str
        Path to .LagrangeMultLog file.
    equilibration : int
        Frames to discard.
    sem_target : float or None
        Optional precision target.
    colvar_id : int or None
        Which CV to use.
    output_dir : Path or None
        Output directory (created if needed).

    Returns
    -------
    dict with keys "report", "diagnostics_png", "csv".
    """
    from ...engines.cp2k import read_constraint_run_from_files

    md_info = read_constraint_run_from_files(restart_path, log_path)
    lambda_series = md_info.lambda_series.collective_shake

    constraint = (
        md_info.metadata.colvars[colvar_id]
        if colvar_id is not None
        else md_info.metadata.colvars.primary
    )
    xi = float(constraint.target_au)
    dt = float(md_info.metadata.timestep_fs)

    t0 = float(md_info.metadata.time_start_fs)
    report = analyze_standalone(
        lambda_series,
        dt=dt,
        xi=xi,
        sem_target=sem_target,
        equilibration=equilibration,
        time_start_fs=t0,
        auto_equilibration=auto_equilibration,
    )

    n_total = len(lambda_series)
    result: dict = {
        "report": report,
        "n_total": n_total,
        "n_analyzed": n_total - equilibration,
        "equilibration": equilibration,
        "dt_fs": dt,
        "time_start_fs": t0 + equilibration * dt,
        "time_end_fs": t0 + (n_total - 1) * dt,
    }

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # Plot
        from .plot import plot_point_diagnostics

        png_path = plot_point_diagnostics(report, output_dir=out)
        result["diagnostics_png"] = png_path

        # CSV
        csv_path = write_single_point_csv(report, output_dir=out)
        result["csv"] = csv_path

    # Console summary
    r = report
    logger.info(
        "Standalone diagnostics for xi=%.6f: tau_corr=%.1f, N_eff=%.1f, "
        "SEM_final=%.6f, Geweke z=%.3f, passed=%s",
        r.xi,
        r.autocorr.tau_corr,
        r.autocorr.n_eff,
        r.sem_final,
        r.geweke.z,
        r.passed,
    )
    if r.failure_reasons:
        for reason in r.failure_reasons:
            logger.info("  - %s", reason)

    return result


# ---------------------------------------------------------------------------
# CSV exports
# ---------------------------------------------------------------------------

_POINT_CSV_COLUMNS = [
    "xi",
    "n_analyzed",
    "time_start_fs",
    "time_end_fs",
    "lambda_mean",
    "sigma_lambda",
    "tau_corr",
    "n_eff",
    "sem_auto",
    "sem_block",
    "delta_sem_block",
    "plateau_B",
    "plateau_reached",
    "sem_final",
    "sem_final_method",
    "sem_max",
    "geweke_z",
    "geweke_reliable",
    "drift_D",
    "passed",
    "failure_reasons",
]


def _point_to_row(r: ConstraintPointReport) -> dict:
    """Convert a ConstraintPointReport to a CSV row dict.

    λ = dA/dξ is NOT pure energy; its unit is Hartree/(ξ_unit).
    All λ-related quantities are kept in a.u. (matching CP2K output).
    """
    ba = r.block_avg
    if ba.plateau_reached:
        method = "plateau"
    else:
        method = "largest_block"
    return {
        "xi": f"{r.xi:.6f}",
        "n_analyzed": str(r.n_analyzed),
        "time_start_fs": f"{r.time_start_fs:.1f}",
        "time_end_fs": f"{r.time_end_fs:.1f}",
        "lambda_mean": f"{r.lambda_mean:.8f}",
        "sigma_lambda": f"{r.sigma_lambda:.8f}",
        "tau_corr": f"{r.autocorr.tau_corr:.2f}",
        "n_eff": f"{r.autocorr.n_eff:.1f}",
        "sem_auto": f"{r.autocorr.sem_auto:.8f}",
        "sem_block": f"{ba.plateau_sem:.8f}",
        "delta_sem_block": f"{ba.plateau_delta:.8f}",
        "plateau_B": (
            str(ba.plateau_block_size) if ba.plateau_block_size is not None
            else "N/A"
        ),
        "plateau_reached": str(ba.plateau_reached),
        "sem_final": f"{r.sem_final:.8f}",
        "sem_final_method": method,
        "sem_max": f"{r.sem_max:.8f}" if r.sem_max is not None else "N/A",
        "geweke_z": f"{r.geweke.z:.4f}",
        "geweke_reliable": str(r.geweke.reliable),
        "drift_D": f"{r.running_avg.drift_D:.8f}",
        "passed": str(r.passed) if r.passed is not None else "N/A",
        "failure_reasons": "; ".join(r.failure_reasons),
    }


def write_single_point_csv(
    report: ConstraintPointReport,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Write single-point diagnostic CSV."""
    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    path = out / DEFAULT_STANDALONE_CSV_NAME

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_POINT_CSV_COLUMNS)
        writer.writeheader()
        writer.writerow(_point_to_row(report))

    return path


def write_convergence_csv(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Write multi-point convergence report CSV."""
    from .config import DEFAULT_REPORT_CSV_NAME

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    path = out / DEFAULT_REPORT_CSV_NAME

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_POINT_CSV_COLUMNS)
        writer.writeheader()
        for r in ti_report.point_reports:
            writer.writerow(_point_to_row(r))

    return path


def write_free_energy_csv(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Write free-energy profile CSV."""
    from .config import DEFAULT_FE_CSV_NAME

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    path = out / DEFAULT_FE_CSV_NAME

    columns = ["xi", "weight", "dA_dxi", "sem", "A_integrated_eV", "sigma_A_cumulative_eV"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()

        # dA/dxi and sem in a.u.; integrated A in eV
        cumul_A = np.cumsum(ti_report.weights * ti_report.forces) * HA_TO_EV
        cumul_sigma = np.sqrt(
            np.cumsum(ti_report.weights**2 * ti_report.force_errors**2)
        ) * HA_TO_EV

        for i in range(len(ti_report.xi_values)):
            writer.writerow(
                {
                    "xi": f"{ti_report.xi_values[i]:.6f}",
                    "weight": f"{ti_report.weights[i]:.8f}",
                    "dA_dxi": f"{ti_report.forces[i]:.8f}",
                    "sem": f"{ti_report.force_errors[i]:.8f}",
                    "A_integrated_eV": f"{cumul_A[i]:.8f}",
                    "sigma_A_cumulative_eV": f"{cumul_sigma[i]:.8f}",
                }
            )

    return path


# ---------------------------------------------------------------------------
# Agent-facing end-to-end wrapper
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TIFullAnalysisReport:
    """Return value of :func:`run_ti_full_from_root`.

    Artifacts (files) and metrics (JSON-serializable values) are carried
    separately so the agent dispatch layer can route them to
    ``TaskResult.outputs`` and ``TaskResult.summary`` respectively.

    ``ti_report`` carries the full :class:`TIReport` object (with embedded
    engine results).  It is *not* JSON-serializable; future Resources-layer
    code may reference it for failure-mode detection.
    """

    # Artifacts
    convergence_csv: Path
    free_energy_csv: Path
    free_energy_png: Path
    diagnostics_pngs: tuple[Path, ...]
    # Metrics (all JSON-serializable)
    n_points: int
    delta_A_eV: float
    sigma_A_eV: float
    all_passed: bool
    failing_indices: tuple[int, ...]
    per_point: tuple[dict[str, Any], ...]
    # Raw model (not MCP-returned; for Resources-layer inspection)
    ti_report: Any   # TIReport — typed as Any to avoid hard import ordering


def _parse_point_slice(spec: str) -> slice:
    """Parse a Python-slice string such as ``"0:2"``, ``":4"``, ``"::2"``.

    Accepts 2–3 colon-separated parts; empty parts become ``None``.
    Rejects forms without a colon (e.g. ``"2"``) and forms with more than
    two colons.  Raises :class:`ValueError` with a user-readable message.
    """
    if ":" not in spec:
        raise ValueError(
            f"Invalid point_slice {spec!r}: must contain at least one ':' "
            "(e.g. '0:2', ':4', '::2')"
        )
    parts = spec.split(":")
    if not (2 <= len(parts) <= 3):
        raise ValueError(
            f"Invalid point_slice {spec!r}: expected 2 or 3 colon-separated "
            f"parts, got {len(parts)}"
        )
    args: list[int | None] = []
    for p in parts:
        stripped = p.strip()
        if stripped == "":
            args.append(None)
        else:
            try:
                args.append(int(stripped))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid point_slice {spec!r}: part {p!r} is not an integer"
                ) from exc
    return slice(*args)


def run_ti_full_from_root(
    root_dir: str | Path = ".",
    output_dir: str | Path = "analysis",
    *,
    parser: str = "auto",
    dir_filter: str | None = None,
    reverse: bool = False,
    equilibration: int | list[int] = 0,
    epsilon_tol_ev: float = DEFAULT_EPSILON_TOL_EV,
    auto_equilibration: bool = False,
    point_slice: str | None = None,
) -> TIFullAnalysisReport:
    """End-to-end constrained-TI analysis from a root directory.

    This is the wrapper behind the ``ti_full_analysis`` agent task.  It
    orchestrates discover → optional slice → load → analyze → CSV / PNG
    → metrics in a single call; the signature matches the agent contract
    exactly (contract-first, see ``agent/_contracts.py``).

    Parameters
    ----------
    root_dir : str or Path
        TI root directory containing constraint-point subdirectories.
    output_dir : str or Path
        Output directory for CSV and PNG files (created if missing).
        Existing files are overwritten.
    parser : str, default ``"auto"``
        Engine parser to use.  ``"auto"`` sniffs registered parsers
        against the first matching subdirectory (currently CP2K is the
        only registered engine).  Pass a registered parser name (e.g.
        ``"cp2k"``) to skip sniffing.
    dir_filter : str or None, default ``None``
        Optional glob pattern to restrict which subdirectories are
        treated as constraint points (e.g. ``"ti_target_*"``).
        ``None`` means: any subdirectory the parser recognises.
    reverse : bool
        If True, treat the max-xi point as the initial state.
    equilibration : int or list[int]
        Frames to discard from the start of each series.  Scalar is broadcast
        to all points; list must have one value per loaded point.
    epsilon_tol_ev : float
        Free-energy tolerance in eV (must be > 0).
    auto_equilibration : bool
        Iteratively discard the front half until convergence (or bottom out).
    point_slice : str or None
        Python-slice syntax string to select a subset of discovered points
        (e.g. ``"0:2"``, ``":4"``, ``"::2"``).  Must be a valid slice with
        2–3 colon-separated parts.  After slicing, at least 2 points must
        remain.

    Returns
    -------
    TIFullAnalysisReport

    Raises
    ------
    FileNotFoundError
        ``root_dir`` does not exist, or a discovered point is missing
        required files.
    ParserInferenceError
        ``parser="auto"`` requested but no registered parser recognises
        any candidate directory under *root_dir*.
    ValueError
        Invalid ``point_slice``; fewer than 2 points after slicing;
        inconsistent ``dt`` across points.
    InsufficientSamplingError
        ``auto_equilibration`` bisected below the min-frames threshold for
        some point.
    """
    # Lazy imports (keep top-level light and avoid import cycles).
    from .io import discover_ti_points, load_ti_series
    from .plot import plot_free_energy_profile, plot_point_diagnostics

    # ── Validate inputs ──────────────────────────────────────────────
    root_path = Path(root_dir)
    if not root_path.is_dir():
        raise FileNotFoundError(
            f"root_dir not found or not a directory: {root_path}"
        )

    # ── 1. Discover constraint points (strict: agent-facing) ─────────
    # strict=True so a matched directory missing required files
    # surfaces as FileNotFoundError rather than being silently skipped
    # (the latter is acceptable for the CLI menu path).
    point_defs = discover_ti_points(
        root_path,
        parser=parser,
        dir_filter=dir_filter,
        reverse=reverse,
        strict=True,
    )

    # ── 2. Optional slice selection (hard-validated) ─────────────────
    if point_slice is not None and point_slice != "":
        sl = _parse_point_slice(point_slice)   # raises ValueError on bad input
        point_defs = point_defs[sl]

    if len(point_defs) < 2:
        raise ValueError(
            f"After slicing, at least 2 TI points are required, got "
            f"{len(point_defs)}. Check root_dir / pattern / point_slice."
        )

    # ── 3. Load series + parse time_starts ───────────────────────────
    series_data = load_ti_series(point_defs)
    xi_values = np.array([x for x, _, _ in series_data])
    lambda_list = [s for _, s, _ in series_data]
    dts = [d for _, _, d in series_data]
    dt = dts[0]
    if not all(abs(d - dt) < 1e-9 for d in dts):
        raise ValueError(
            f"Inconsistent dt across TI points: {dts}. All points must share "
            "the same frame interval."
        )
    time_starts = [float(p.metadata.time_start_fs) for p in point_defs]

    # ── 4. Analyze ───────────────────────────────────────────────────
    ti_report = analyze_ti(
        xi_values,
        lambda_list,
        dt,
        epsilon_tol_ev=epsilon_tol_ev,
        equilibration=equilibration,
        time_starts=time_starts,
        auto_equilibration=auto_equilibration,
    )

    # ── 5. Write outputs ─────────────────────────────────────────────
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    csv_conv = write_convergence_csv(ti_report, output_dir=output_path)
    csv_fe = write_free_energy_csv(ti_report, output_dir=output_path)
    png_fe = plot_free_energy_profile(ti_report, output_dir=output_path)
    diag_pngs = tuple(
        plot_point_diagnostics(r, output_dir=output_path)
        for r in ti_report.point_reports
    )

    # ── 6. Build per-point metrics (JSON-serializable) ───────────────
    per_point: list[dict[str, Any]] = []
    for idx, rep in enumerate(ti_report.point_reports):
        per_point.append({
            "point_index": rep.point_index if rep.point_index is not None else idx,
            "xi": float(rep.xi),
            "n_analyzed": int(rep.n_analyzed),
            "time_start_fs": float(rep.time_start_fs),
            "time_end_fs": float(rep.time_end_fs),
            "time_total_fs": float(rep.n_analyzed * dt),
            "tau_corr": float(rep.autocorr.tau_corr),
            "n_eff": float(rep.autocorr.n_eff),
            "sem_final_au": float(rep.sem_final),
            "sem_max_au": (float(rep.sem_max) if rep.sem_max is not None else None),
            "geweke_z": float(rep.geweke.z),
            "drift_D": float(rep.running_avg.drift_D),
            "passed": (bool(rep.passed) if rep.passed is not None else None),
            "failure_reasons": list(rep.failure_reasons),
        })

    # ── 7. Assemble report ───────────────────────────────────────────
    return TIFullAnalysisReport(
        convergence_csv=csv_conv,
        free_energy_csv=csv_fe,
        free_energy_png=png_fe,
        diagnostics_pngs=diag_pngs,
        n_points=len(ti_report.point_reports),
        delta_A_eV=round(ti_report.delta_A * HA_TO_EV, 6),
        sigma_A_eV=round(ti_report.sigma_A * HA_TO_EV, 6),
        all_passed=bool(ti_report.all_passed),
        failing_indices=tuple(int(i) for i in ti_report.failing_indices),
        per_point=tuple(per_point),
        ti_report=ti_report,
    )
