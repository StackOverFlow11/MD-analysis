"""Task handler factory and Phase 1 task registrations.

Convention: all analysis modules MUST be imported inside function bodies
(lazy import), never at module top level.  The ``_make_handler`` factory
enforces this structurally via ``importlib.import_module`` in the closure.
"""

from __future__ import annotations

import logging
from importlib import import_module
from pathlib import Path
from typing import Any, Callable

from ._core import TaskDef, TaskResult, register

logger = logging.getLogger(__name__)

# Type alias for summary extraction functions
SummaryExtractor = Callable[[Any, dict[str, Any]], dict[str, Any]]


# ── Handler factory ──────────────────────────────────────────────────


def _make_handler(
    task_name: str,
    target_fn_path: str,
    summary_extractor: SummaryExtractor | None = None,
) -> Callable[[dict[str, Any]], TaskResult]:
    """Generate a pass-through handler that calls the target function.

    Exception handling is done by the ``dispatch()`` layer — handlers
    produced here do NOT wrap calls in try-except.
    """

    def handler(params: dict[str, Any]) -> TaskResult:
        module_path, fn_name = target_fn_path.rsplit(":", 1)
        fn = getattr(import_module(module_path), fn_name)

        result = fn(**params)

        outputs = _normalize_outputs(result)
        summary = summary_extractor(result, params) if summary_extractor else {}

        return TaskResult(
            success=True,
            task=task_name,
            outputs=outputs,
            summary=summary,
        )

    return handler


def _normalize_outputs(result: Any) -> dict[str, str]:
    """Normalize diverse return types to ``{name: path_string}``."""
    if isinstance(result, dict):
        return {k: str(v) for k, v in result.items()}
    if isinstance(result, Path):
        return {"output": str(result)}
    if isinstance(result, list):
        # list[Path] from batch scripts → enumerate as workdir_0, workdir_1, ...
        return {f"workdir_{i}": str(p) for i, p in enumerate(result)}
    # SurfaceChargeResult and similar dataclasses with csv_path
    if hasattr(result, "csv_path"):
        out: dict[str, str] = {"csv": str(result.csv_path)}
        png = result.csv_path.parent / (result.csv_path.stem + ".png")
        if png.exists():
            out["png"] = str(png)
        return out
    logger.debug(
        "Unrecognized return type %s, outputs will be empty",
        type(result).__name__,
    )
    return {}


# ── Phase 1 task registrations (10 high-frequency tasks) ────────────

# 1. water_three_panel (CLI 105)
register(TaskDef(
    name="water_three_panel",
    category="water",
    description="Water density, orientation, adsorbed layer three-panel analysis",
    handler=_make_handler(
        "water_three_panel", "md_analysis.main:run_water_analysis",
    ),
    target_fn="md_analysis.main:run_water_analysis",
    cli_codes=("105",),
    param_descriptions={
        "xyz_path": "XYZ trajectory file",
        "md_inp_path": "CP2K md.inp (for cell_abc)",
        "cell_abc": "Cell dimensions [a, b, c] in Angstrom",
        "output_dir": "Output directory",
        "layer_tol_A": "Layer detection tolerance (Angstrom)",
    },
))

# 2. potential_full (CLI 216)
register(TaskDef(
    name="potential_full",
    category="potential",
    description="Full potential analysis (center + fermi + electrode + phi_z + sensitivity)",
    handler=_make_handler(
        "potential_full", "md_analysis.main:run_potential_analysis",
    ),
    target_fn="md_analysis.main:run_potential_analysis",
    cli_codes=("216",),
    param_descriptions={
        "output_dir": "Output directory",
        "cube_pattern": "Glob pattern for cube files",
        "md_out_path": "CP2K md.out for Fermi energy",
        "thickness_ang": "Slab thickness for averaging (Angstrom)",
        "input_mode": "continuous (MD) or distributed (SP subdirs)",
    },
    param_choices={
        "input_mode": ["continuous", "distributed"],
        "center_mode": ["interface", "cell"],
        "fermi_unit": ["au", "ev"],
    },
))

# 3. charge_surface (CLI 221-223)
register(TaskDef(
    name="charge_surface",
    category="charge",
    description="Bader surface charge density time series",
    handler=_make_handler(
        "charge_surface", "md_analysis.main:run_charge_analysis",
    ),
    target_fn="md_analysis.main:run_charge_analysis",
    cli_codes=("221", "222", "223"),
    param_descriptions={
        "output_dir": "Output directory",
        "method": "Charge partitioning method (defaults to counterion)",
        "root_dir": "Bader calculation root directory",
    },
    param_choices={
        "method": ["counterion", "layer"],
        "normal": ["a", "b", "c"],
    },
))

# 4. charge_tracked (CLI 225)
register(TaskDef(
    name="charge_tracked",
    category="charge",
    description="Track Bader charges for specified atoms",
    handler=_make_handler(
        "charge_tracked", "md_analysis.main:run_tracked_charge_analysis",
    ),
    target_fn="md_analysis.main:run_tracked_charge_analysis",
    cli_codes=("225",),
    param_descriptions={
        "output_dir": "Output directory",
        "atom_indices_xyz": "0-based XYZ atom indices to track (list of int)",
    },
))

# 5. charge_counterion (CLI 226)
register(TaskDef(
    name="charge_counterion",
    category="charge",
    description="Per-frame counterion detection and charge tracking",
    handler=_make_handler(
        "charge_counterion", "md_analysis.main:run_counterion_charge_analysis",
    ),
    target_fn="md_analysis.main:run_counterion_charge_analysis",
    cli_codes=("226",),
    param_descriptions={
        "output_dir": "Output directory",
        "root_dir": "Bader calculation root directory",
    },
    param_choices={
        "normal": ["a", "b", "c"],
    },
))

# 6. run_all (composite)
register(TaskDef(
    name="run_all",
    category="composite",
    description="Run water + potential analysis with standard directory layout",
    handler=_make_handler("run_all", "md_analysis.main:run_all"),
    target_fn="md_analysis.main:run_all",
    param_descriptions={
        "xyz_path": "XYZ trajectory file",
        "output_dir": "Output root directory",
        "cube_pattern": "Glob pattern for cube files",
        "md_out_path": "CP2K md.out for Fermi energy",
    },
))

# 7. calibration_fit_csv (CLI 231)
register(TaskDef(
    name="calibration_fit_csv",
    category="calibration",
    description="Fit sigma-phi calibration curve from CSV data",
    handler=_make_handler(
        "calibration_fit_csv",
        "md_analysis.electrochemical.calibration.CalibrationWorkflow:calibrate",
    ),
    target_fn="md_analysis.electrochemical.calibration.CalibrationWorkflow:calibrate",
    cli_codes=("231",),
    param_descriptions={
        "csv_path": "CSV file with sigma and phi columns",
        "output_dir": "Output directory for calibration results",
    },
    param_choices={
        "fitting_method": ["linear", "polynomial", "spline"],
    },
))

# 8. calibration_predict (CLI 233)
register(TaskDef(
    name="calibration_predict",
    category="calibration",
    description="Predict electrode potential from sigma using calibration",
    handler=_make_handler(
        "calibration_predict",
        "md_analysis.electrochemical.calibration.CalibrationWorkflow:predict_potential",
    ),
    target_fn="md_analysis.electrochemical.calibration.CalibrationWorkflow:predict_potential",
    cli_codes=("233",),
    param_descriptions={
        "calibration_json": "Path to calibration.json",
        "sigma_value": "Surface charge density (uC/cm2)",
    },
))

# 9. slowgrowth_quick (CLI 301)
register(TaskDef(
    name="slowgrowth_quick",
    category="enhanced_sampling",
    description="Quick slow-growth free energy integration plot",
    handler=_make_handler(
        "slowgrowth_quick",
        "md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot:slowgrowth_analysis",
    ),
    target_fn="md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot:slowgrowth_analysis",
    cli_codes=("301",),
    param_descriptions={
        "log_path": "LagrangeMultLog file path",
        "restart_path": "CP2K restart file",
        "output_dir": "Output directory",
    },
))

# 10. config_show (CLI 900)
def _handle_config_show(params: dict[str, Any]) -> TaskResult:
    """Show current user configuration."""
    from ..config import load_config
    config = load_config()
    return TaskResult(
        success=True,
        task="config_show",
        outputs={},
        summary={"config": config},
    )


register(TaskDef(
    name="config_show",
    category="meta",
    description="Show current user configuration",
    handler=_handle_config_show,
    target_fn="md_analysis.config:load_config",
    cli_codes=("900",),
))


# 11. ti_full_analysis (CLI 312)
def _handle_ti_full_analysis(params: dict[str, Any]) -> TaskResult:
    """Full constrained-TI convergence analysis + free-energy integration."""
    import numpy as np

    from ..enhanced_sampling.constrained_ti.io import (
        discover_ti_points,
        load_ti_series,
    )
    from ..enhanced_sampling.constrained_ti.workflow import (
        analyze_ti,
        write_convergence_csv,
        write_free_energy_csv,
    )
    from ..enhanced_sampling.constrained_ti.plot import (
        plot_free_energy_profile,
        plot_point_diagnostics,
    )
    from ..utils.RestartParser.ColvarParser import parse_colvar_restart
    from ..utils.constants import HA_TO_EV

    root_dir = Path(params.get("root_dir", "."))
    output_dir = Path(params.get("output_dir", "analysis"))
    pattern = params.get("pattern", "auto")
    reverse = params.get("reverse", False)
    equilibration = params.get("equilibration", 0)
    epsilon_tol_ev = params.get("epsilon_tol_ev", 0.05)
    auto_equilibration = params.get("auto_equilibration", False)
    point_slice = params.get("point_slice", None)

    # 1. Discover constraint points
    point_defs = discover_ti_points(root_dir, pattern=pattern, reverse=reverse)

    # 2. Optional slice selection (e.g. "3:8", "::2")
    if point_slice:
        parts = point_slice.split(":")
        args = [int(x) if x.strip() else None for x in parts]
        point_defs = point_defs[slice(*args)]

    # 3. Load series + parse time_starts
    series_data = load_ti_series(point_defs)
    xi_values = np.array([x for x, _, _ in series_data])
    lambda_list = [s for _, s, _ in series_data]
    dts = [d for _, _, d in series_data]
    dt = dts[0]
    time_starts = [
        parse_colvar_restart(str(p.restart_path)).time_start_fs
        for p in point_defs
    ]

    # 4. Analyze
    ti_report = analyze_ti(
        xi_values,
        lambda_list,
        dt,
        epsilon_tol_ev=epsilon_tol_ev,
        equilibration=equilibration,
        time_starts=time_starts,
        auto_equilibration=auto_equilibration,
    )

    # 5. Write outputs (plot_free_energy_profile shows CV ticks on x-axis)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_conv = write_convergence_csv(ti_report, output_dir=output_dir)
    csv_fe = write_free_energy_csv(ti_report, output_dir=output_dir)
    png_fe = plot_free_energy_profile(ti_report, output_dir=output_dir)
    diag_pngs = []
    for r in ti_report.point_reports:
        diag_pngs.append(plot_point_diagnostics(r, output_dir=output_dir))

    # 6. Build outputs dict
    outputs: dict[str, str] = {
        "convergence_csv": str(csv_conv),
        "free_energy_csv": str(csv_fe),
        "free_energy_png": str(png_fe),
    }
    for i, p in enumerate(diag_pngs):
        outputs[f"diagnostics_png_{i}"] = str(p)

    # 7. Summary
    delta_A_eV = ti_report.delta_A * HA_TO_EV
    sigma_A_eV = ti_report.sigma_A * HA_TO_EV
    summary = {
        "n_points": len(ti_report.point_reports),
        "delta_A_eV": round(delta_A_eV, 6),
        "sigma_A_eV": round(sigma_A_eV, 6),
        "all_passed": ti_report.all_passed,
        "failing_indices": list(ti_report.failing_indices),
    }

    return TaskResult(
        success=True,
        task="ti_full_analysis",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="ti_full_analysis",
    category="enhanced_sampling",
    description="Full constrained-TI convergence analysis + free-energy integration",
    handler=_handle_ti_full_analysis,
    target_fn="md_analysis.enhanced_sampling.constrained_ti.workflow:analyze_ti",
    cli_codes=("312",),
    param_descriptions={
        "root_dir": "TI root directory containing constraint-point subdirs",
        "output_dir": "Output directory for CSV and PNG files",
        "pattern": "Directory discovery pattern (ti_target/xi/auto)",
        "reverse": "Reverse integration direction (initial state = max xi)",
        "equilibration": "Frames to discard from start (int or list[int])",
        "epsilon_tol_ev": "Free-energy tolerance in eV",
        "auto_equilibration": "Iteratively discard first half until converged",
        "point_slice": "Python slice string to select points (e.g. 3:8, ::2)",
    },
    param_choices={
        "pattern": ["ti_target", "xi", "auto"],
    },
))


# 12. sp_gen_batch (CLI 442) — Batch-generate CP2K SP work directories for DeePMD training
register(TaskDef(
    name="sp_gen_batch",
    category="scripts",
    description=(
        "Batch-generate CP2K SP work directories for DeePMD training. "
        "Supports two frame-selection modes: mode='index' uses frame_start/end/step "
        "(0-based indices); mode='time' uses time_start_fs/end_fs/step_fs (all three "
        "required together, greedy matching on atoms.info['time'])."
    ),
    handler=_make_handler(
        "sp_gen_batch",
        "md_analysis.scripts.SpGen:batch_generate_sp_workdirs",
    ),
    target_fn="md_analysis.scripts.SpGen:batch_generate_sp_workdirs",
    cli_codes=("442",),
    param_descriptions={
        "xyz_path": "CP2K XYZ trajectory file",
        "cell_abc": "Cell dimensions [a, b, c] in Angstrom",
        "output_dir": "Parent output directory (sub-dirs sp_t{time}_i{step}/ created here)",
        "inp_template_path": "DP SP inp template (fallback: KEY_DP_SP_INP_TEMPLATE_PATH)",
        "mode": "Frame selection mode: 'index' (use frame_*) or 'time' (use time_*_fs)",
        "frame_start": "Index mode: first frame index (0-based)",
        "frame_end": "Index mode: exclusive upper bound (None = all)",
        "frame_step": "Index mode: step between frames",
        "time_start_fs": "Time mode: start time in fs (inclusive). Required when mode='time'",
        "time_end_fs": "Time mode: end time in fs (inclusive). Required when mode='time'",
        "time_step_fs": "Time mode: time step in fs (greedy). Required when mode='time'",
        "script_path": "Submission script to copy (fallback: KEY_CP2K_SCRIPT_PATH)",
    },
    param_choices={
        "mode": ["index", "time"],
    },
))
