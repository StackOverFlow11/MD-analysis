"""Orchestrator for constrained-TI convergence diagnostics.

Hosts:
- ``analyze_single_point``  — single point in TI context
- ``analyze_standalone``    — standalone single-point (no TI context)
- ``analyze_ti``            — full multi-point TI analysis
- ``standalone_diagnostics`` — unified entry: parse + analyze + plot + CSV
- CSV export helpers
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np

from .analysis.autocorrelation import analyze_autocorrelation
from .analysis.block_average import analyze_block_average
from .analysis.geweke import analyze_geweke
from .analysis.running_average import analyze_running_average
from .config import (
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
from ...utils.config import HA_TO_EV
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

    # Determine sem_final: 2-tier hierarchy (F&P plateau → ACF fallback)
    if block_avg.plateau_reached:
        sem_final = block_avg.plateau_sem
    else:
        sem_final = autocorr.sem_auto
        failure_reasons.append(
            "F&P block-average plateau not reached; falling back to SEM_auto."
        )

    # Cross-check: SEM_block vs SEM_auto
    if block_avg.plateau_reached:
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

    return ConstraintPointReport(
        xi=inp.xi,
        point_index=inp.point_index,
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
# Standalone single-point
# ---------------------------------------------------------------------------


def analyze_standalone(
    lambda_series: np.ndarray,
    *,
    dt: float = 1.0,
    xi: float = 0.0,
    sem_target: float | None = None,
    equilibration: int = 0,
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
    **engine_overrides
        Forwarded to analysis engines.

    Returns
    -------
    ConstraintPointReport
    """
    series = lambda_series.copy()
    series, pre_warnings = _validate_and_trim(series, equilibration=equilibration)

    inp = ConstraintPointInput(
        xi=xi,
        lambda_series=series,
        dt=dt,
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

    # Analyze each point
    reports: list[ConstraintPointReport] = []
    for i in range(k):
        series = lambda_series_list[i].copy()
        series, pre_warnings = _validate_and_trim(
            series, equilibration=equil_list[i]
        )

        inp = ConstraintPointInput(
            xi=float(xi_values[i]),
            lambda_series=series,
            dt=dt,
            weight=float(weights[i]),
            sem_max=float(sem_targets[i]),
            point_index=i,
        )
        report = analyze_single_point(inp, **engine_overrides)
        if pre_warnings:
            report = ConstraintPointReport(
                xi=report.xi,
                point_index=report.point_index,
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
    forces = np.array([r.lambda_mean for r in reports])
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
    from ...utils.RestartParser.ColvarParser import ColvarMDInfo

    md_info = ColvarMDInfo.from_paths(restart_path, log_path)
    lambda_series = md_info.lagrange.collective_shake

    constraint = (
        md_info.restart.colvars[colvar_id]
        if colvar_id is not None
        else md_info.restart.colvars.primary
    )
    xi = float(constraint.target_au)
    dt = float(md_info.restart.timestep_fs)

    report = analyze_standalone(
        lambda_series,
        dt=dt,
        xi=xi,
        sem_target=sem_target,
        equilibration=equilibration,
    )

    result: dict[str, Path | ConstraintPointReport] = {"report": report}

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
        method = "acf"
    return {
        "xi": f"{r.xi:.6f}",
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
