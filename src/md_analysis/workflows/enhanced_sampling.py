"""Enhanced-sampling workflow facade.

Exposes five programmatic entry points covering the slow-growth (SG)
and constrained-thermodynamic-integration (TI) analysis surfaces:

- :func:`run_slowgrowth_quick_plot`       — SG quick diagnostic plot + CSV
- :func:`run_slowgrowth_publication_plot` — SG publication-style plot + CSV
- :func:`run_ti_single_diagnostics`       — single-point TI 2x2 diagnostics
- :func:`run_ti_full_analysis`            — multi-point TI end-to-end
- :func:`run_ti_constant_potential_correction` — TI full + Nørskov correction

Business logic lives under
:mod:`md_analysis.enhanced_sampling.{slowgrowth,constrained_ti}`; this
module only organises parameters, output paths, and the
:class:`WorkflowResult` contract.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..enhanced_sampling.constrained_ti.config import DEFAULT_EPSILON_TOL_EV
from .models import WorkflowResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Slow-growth workflows
# ---------------------------------------------------------------------------


def _run_slowgrowth_with_style(
    *,
    restart_path: Path | str,
    log_path: Path | str,
    plot_style: str,
    initial_step: int,
    final_step: int | None,
    output_dir: Path | str,
    colvar_id: int | None,
    workflow_name: str,
) -> WorkflowResult:
    """Shared SG workflow body used by the quick / publication facades."""
    from ..enhanced_sampling.slowgrowth.SlowGrowthPlot import (
        slowgrowth_analysis_with_report,
    )

    restart_p = Path(restart_path)
    log_p = Path(log_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting slow-growth analysis: plot_style=%s, output_dir=%s",
        plot_style,
        out_dir,
    )

    report = slowgrowth_analysis_with_report(
        str(restart_p),
        str(log_p),
        initial_step=initial_step,
        final_step=final_step,
        output_dir=out_dir,
        plot_style=plot_style,
        colvar_id=colvar_id,
    )

    artifacts: dict[str, Path] = {
        key: Path(path) for key, path in report.artifacts.items()
    }
    metadata: dict[str, Any] = {
        "plot_style": plot_style,
        "n_steps": report.n_steps,
        "target_start_au": report.target_start_au,
        "target_end_au": report.target_end_au,
        "delta_F_eV": report.delta_F_eV,
        "delta_F_barrier_eV": report.delta_F_barrier_eV,
        "barrier_step": report.barrier_step,
        "is_reversed": report.is_reversed,
        "initial_step": initial_step,
        "final_step": final_step,
        "colvar_id": colvar_id,
    }
    return WorkflowResult(
        name=workflow_name,
        output_dir=out_dir,
        artifacts=artifacts,
        metadata=metadata,
        extra=report,
    )


def run_slowgrowth_quick_plot(
    *,
    restart_path: Path | str,
    log_path: Path | str,
    output_dir: Path | str,
    initial_step: int = 0,
    final_step: int | None = None,
    colvar_id: int | None = None,
) -> WorkflowResult:
    """SG quick plot + CSV (dual-axis: FE / Lagrange).

    Artifacts: ``csv``, ``quick_png``. Convergence / barrier metrics
    are reported via ``metadata`` and the strongly-typed
    :class:`SlowgrowthAnalysisReport` on ``extra``.
    """
    return _run_slowgrowth_with_style(
        restart_path=restart_path,
        log_path=log_path,
        plot_style="quick",
        initial_step=initial_step,
        final_step=final_step,
        output_dir=output_dir,
        colvar_id=colvar_id,
        workflow_name="slowgrowth_quick_plot",
    )


def run_slowgrowth_publication_plot(
    *,
    restart_path: Path | str,
    log_path: Path | str,
    output_dir: Path | str,
    initial_step: int = 0,
    final_step: int | None = None,
    colvar_id: int | None = None,
) -> WorkflowResult:
    """SG publication-style plot + CSV.

    Artifacts: ``csv``, ``publication_png``. Same metadata surface as
    :func:`run_slowgrowth_quick_plot`.
    """
    return _run_slowgrowth_with_style(
        restart_path=restart_path,
        log_path=log_path,
        plot_style="publication",
        initial_step=initial_step,
        final_step=final_step,
        output_dir=output_dir,
        colvar_id=colvar_id,
        workflow_name="slowgrowth_publication_plot",
    )


# ---------------------------------------------------------------------------
# Constrained-TI workflows
# ---------------------------------------------------------------------------


def run_ti_single_diagnostics(
    *,
    restart_path: Path | str,
    log_path: Path | str,
    output_dir: Path | str,
    equilibration: int = 0,
    sem_target: float | None = None,
    colvar_id: int | None = None,
    auto_equilibration: bool = False,
) -> WorkflowResult:
    """Run the 2x2 convergence diagnostics for a single constraint point.

    Artifacts: ``csv``, ``diagnostics_png``. The full
    :class:`ConstraintPointReport` is carried on ``extra``; a small
    subset of its scalar fields (n_eff, tau_corr, geweke_z, passed,
    etc.) is mirrored to ``metadata`` for quick scanning.
    """
    from ..enhanced_sampling.constrained_ti.workflow import (
        standalone_diagnostics,
    )

    restart_p = Path(restart_path)
    log_p = Path(log_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting TI single-point diagnostics: restart=%s, output_dir=%s",
        restart_p,
        out_dir,
    )

    result = standalone_diagnostics(
        str(restart_p),
        str(log_p),
        equilibration=equilibration,
        sem_target=sem_target,
        colvar_id=colvar_id,
        output_dir=out_dir,
        auto_equilibration=auto_equilibration,
    )

    point_report = result["report"]
    artifacts: dict[str, Path] = {}
    for key in ("csv", "diagnostics_png"):
        value = result.get(key)
        if value is not None:
            artifacts[key] = Path(value)

    metadata: dict[str, Any] = {
        "n_total": result["n_total"],
        "n_analyzed": result["n_analyzed"],
        "equilibration": result["equilibration"],
        "dt_fs": result["dt_fs"],
        "time_start_fs": result["time_start_fs"],
        "time_end_fs": result["time_end_fs"],
        "xi": float(point_report.xi),
        "tau_corr": float(point_report.autocorr.tau_corr),
        "n_eff": float(point_report.autocorr.n_eff),
        "geweke_z": float(point_report.geweke.z),
        "passed": (
            bool(point_report.passed) if point_report.passed is not None else None
        ),
        "sem_target": sem_target,
        "auto_equilibration": auto_equilibration,
    }
    return WorkflowResult(
        name="ti_single_diagnostics",
        output_dir=out_dir,
        artifacts=artifacts,
        metadata=metadata,
        extra=point_report,
    )


def run_ti_full_analysis(
    *,
    root_dir: Path | str,
    output_dir: Path | str,
    parser: str = "auto",
    dir_filter: str | None = None,
    reverse: bool = False,
    equilibration: int | list[int] = 0,
    epsilon_tol_ev: float = DEFAULT_EPSILON_TOL_EV,
    auto_equilibration: bool = False,
    point_slice: str | None = None,
) -> WorkflowResult:
    """End-to-end multi-point constrained TI analysis.

    Discovers constraint points under ``root_dir``, loads each
    series, applies the four-step convergence diagnostics, integrates
    ⟨λ⟩ over ξ for the free-energy profile, and writes the
    convergence CSV / free-energy CSV / free-energy PNG / per-point
    diagnostic PNGs.

    Per-point diagnostic PNGs are flattened into the artifact map as
    ``diagnostics_png_0``, ``diagnostics_png_1``, ... (the
    :class:`WorkflowResult` artifact contract is a flat
    ``dict[str, Path]``; tuples are flattened by integer index). The
    full :class:`TIFullAnalysisReport` is preserved on ``extra``.
    """
    from ..enhanced_sampling.constrained_ti.workflow import (
        run_ti_full_from_root,
    )

    root_p = Path(root_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting TI full analysis: root=%s, output_dir=%s, point_slice=%s",
        root_p,
        out_dir,
        point_slice,
    )

    report = run_ti_full_from_root(
        root_dir=root_p,
        output_dir=out_dir,
        parser=parser,
        dir_filter=dir_filter,
        reverse=reverse,
        equilibration=equilibration,
        epsilon_tol_ev=epsilon_tol_ev,
        auto_equilibration=auto_equilibration,
        point_slice=point_slice,
    )

    artifacts: dict[str, Path] = {
        "convergence_csv": Path(report.convergence_csv),
        "free_energy_csv": Path(report.free_energy_csv),
        "free_energy_png": Path(report.free_energy_png),
    }
    for i, png in enumerate(report.diagnostics_pngs):
        artifacts[f"diagnostics_png_{i}"] = Path(png)

    metadata: dict[str, Any] = {
        "n_points": report.n_points,
        "delta_A_eV": report.delta_A_eV,
        "sigma_A_eV": report.sigma_A_eV,
        "all_passed": report.all_passed,
        "failing_indices": list(report.failing_indices),
        "parser": parser,
        "dir_filter": dir_filter,
        "reverse": reverse,
        "point_slice": point_slice,
        "auto_equilibration": auto_equilibration,
        "epsilon_tol_ev": epsilon_tol_ev,
    }
    return WorkflowResult(
        name="ti_full_analysis",
        output_dir=out_dir,
        artifacts=artifacts,
        metadata=metadata,
        extra=report,
    )


def run_ti_constant_potential_correction(
    *,
    root_dir: Path | str,
    output_dir: Path | str,
    calibration_json_path: Path | str,
    target_side: str = "aligned",
    method: str = "counterion",
    normal: str = "c",
    parser: str = "auto",
    dir_filter: str | None = None,
    reverse: bool = False,
    equilibration: int | list[int] = 0,
    epsilon_tol_ev: float = DEFAULT_EPSILON_TOL_EV,
    auto_equilibration: bool = False,
    point_slice: str | None = None,
) -> WorkflowResult:
    """Run TI full analysis and apply the Nørskov constant-potential correction.

    Pipeline:

    1. Full constrained-TI analysis (delegated to
       :func:`run_ti_full_analysis`).
    2. Re-discover constraint-point directories (lightweight) and apply
       the same ``point_slice`` so per-point Bader data alignment
       matches the analysed report.
    3. Load the calibration JSON, reconstruct the σ→φ mapper, and call
       :func:`compute_constant_potential_correction`.
    4. Write the corrected free-energy CSV + PNG into ``output_dir``.

    Artifacts include all TI baseline files plus
    ``corrected_free_energy_csv`` and ``corrected_free_energy_png``.
    Both the TI baseline metrics and the correction summary are
    surfaced via ``metadata``; the strongly-typed
    :class:`ConstantPotentialResult` (which embeds the TIReport) is
    placed on ``extra``.
    """
    from ..electrochemical.calibration._data import load_calibration_json
    from ..electrochemical.calibration._mapper import mapper_from_dict
    from ..enhanced_sampling.constrained_ti.correction import (
        compute_constant_potential_correction,
        plot_corrected_free_energy_profile,
        write_corrected_free_energy_csv,
    )
    from ..enhanced_sampling.constrained_ti.io import discover_ti_points
    from ..enhanced_sampling.constrained_ti.workflow import (
        _parse_point_slice,
    )

    root_p = Path(root_dir)
    out_dir = Path(output_dir)
    cal_path = Path(calibration_json_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: TI baseline analysis (writes ti_*.csv, ti_*.png, diag pngs)
    ti_result = run_ti_full_analysis(
        root_dir=root_p,
        output_dir=out_dir,
        parser=parser,
        dir_filter=dir_filter,
        reverse=reverse,
        equilibration=equilibration,
        epsilon_tol_ev=epsilon_tol_ev,
        auto_equilibration=auto_equilibration,
        point_slice=point_slice,
    )
    ti_report = ti_result.extra.ti_report  # underlying TIReport

    # Phase 2: re-discover point_defs (cheap; no series reload) and slice
    # the same way run_ti_full_from_root did so indices align with
    # ti_report.point_reports.
    point_defs = discover_ti_points(
        root_p, parser=parser, dir_filter=dir_filter, reverse=reverse,
    )
    if point_slice is not None:
        point_defs = point_defs[_parse_point_slice(point_slice)]

    # Phase 3: load mapper from calibration JSON.
    _cal_data, fit_params = load_calibration_json(cal_path)
    mapper = mapper_from_dict(fit_params)

    # Phase 4: correction.
    logger.info(
        "Applying constant-potential correction: target_side=%s, method=%s",
        target_side,
        method,
    )
    correction_result = compute_constant_potential_correction(
        ti_report,
        point_defs,
        mapper,
        target_side=target_side,
        method=method,
        normal=normal,
    )
    csv_corr = write_corrected_free_energy_csv(
        correction_result, output_dir=out_dir
    )
    png_corr = plot_corrected_free_energy_profile(
        correction_result, output_dir=out_dir
    )

    artifacts = dict(ti_result.artifacts)
    artifacts["corrected_free_energy_csv"] = Path(csv_corr)
    artifacts["corrected_free_energy_png"] = Path(png_corr)

    corr = correction_result.correction
    metadata: dict[str, Any] = {
        **ti_result.metadata,
        "target_side": target_side,
        "method": method,
        "normal": normal,
        "calibration_json": str(cal_path),
        "area_A2": corr.area_A2,
        "delta_A_const_phi_eV": correction_result.delta_A_const_phi_eV,
        "delta_A_const_q_eV": float(correction_result.A_const_q_eV[-1]),
        "total_correction_eV": (
            correction_result.delta_A_const_phi_eV
            - float(correction_result.A_const_q_eV[-1])
        ),
    }
    return WorkflowResult(
        name="ti_constant_potential_correction",
        output_dir=out_dir,
        artifacts=artifacts,
        metadata=metadata,
        extra=correction_result,
    )
