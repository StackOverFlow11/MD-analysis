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

from ._contracts import ExceptionMapping, FieldSpec, TaskContract
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
    # Dataclass reports with a `workdirs` tuple field (e.g. TIGenBatchReport)
    wds = getattr(result, "workdirs", None)
    if wds is not None and isinstance(wds, (tuple, list)):
        return {f"workdir_{i}": str(p) for i, p in enumerate(wds)}
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


# 11. ti_full_analysis (CLI 312) — full convergence analysis + ΔA integration.
#     Second task to ship with a full TaskContract (composite handler backed by
#     the real wrapper run_ti_full_from_root).


_TI_FULL_CONTRACT = TaskContract(
    inputs={
        "root_dir": FieldSpec(
            description="TI root directory containing constraint-point subdirs",
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
            required=False, default=".",
        ),
        "output_dir": FieldSpec(
            description="Output directory for CSV and PNG files (created if missing; existing files overwritten)",
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
            required=False, default="analysis",
        ),
        "pattern": FieldSpec(
            description="Directory discovery pattern",
            json_schema={"type": "string", "enum": ["auto", "ti_target", "xi"]},
            type="str",
            choices=("auto", "ti_target", "xi"),
            required=False, default="auto",
        ),
        "reverse": FieldSpec(
            description="Reverse integration direction (initial state = max xi)",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
        "equilibration": FieldSpec(
            description="Frames to discard from start (scalar or per-point list)",
            json_schema={
                "oneOf": [
                    {"type": "integer", "minimum": 0},
                    {"type": "array", "items": {"type": "integer", "minimum": 0}},
                ],
            },
            type="int | list[int]",
            required=False, default=0,
        ),
        "epsilon_tol_ev": FieldSpec(
            description="Free-energy tolerance in eV (must be > 0)",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="eV",
            required=False, default=0.05,
        ),
        "auto_equilibration": FieldSpec(
            description="Iteratively discard first half until converged",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
        "point_slice": FieldSpec(
            description=(
                "Python slice string to select subset of discovered points "
                "(e.g. '0:2', ':4', '::2'). After slicing, at least 2 points "
                "must remain."
            ),
            json_schema={"type": ["string", "null"]},
            type="str | None",
            required=False, default=None,
        ),
    },
    outputs_artifacts={
        "convergence_csv": FieldSpec(
            description="Per-point diagnostics CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "free_energy_csv": FieldSpec(
            description="Integrated ΔA vs xi CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "free_energy_png": FieldSpec(
            description="Free-energy curve figure",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "diagnostics_pngs": FieldSpec(
            description="2x2 diagnostic plots, one per constraint point",
            json_schema={"type": "array", "items": {"type": "string"}},
            type="list[Path]", path_kind="file", shape="(N,)",
            category="artifact",
        ),
    },
    outputs_metrics={
        "n_points": FieldSpec(
            description="Number of constraint points analyzed",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "delta_A_eV": FieldSpec(
            description="Total free-energy change",
            json_schema={"type": "number"},
            type="float", unit="eV", category="metric",
        ),
        "sigma_A_eV": FieldSpec(
            description="Propagated SEM on ΔA",
            json_schema={"type": "number", "minimum": 0},
            type="float", unit="eV", category="metric",
        ),
        "all_passed": FieldSpec(
            description="All constraint points passed the 4-step diagnostics",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "failing_indices": FieldSpec(
            description="Indices (in discovered order) of points that failed",
            json_schema={"type": "array", "items": {"type": "integer"}},
            type="list[int]", shape="(K,)", category="metric",
        ),
        "per_point": FieldSpec(
            description=(
                "Per-point summary supporting future Resources-layer "
                "failure-mode classification. Fields: point_index, xi, "
                "n_analyzed, time_start_fs, time_end_fs, time_total_fs, "
                "tau_corr, n_eff, sem_final_au, sem_max_au, geweke_z, "
                "drift_D, passed, failure_reasons."
            ),
            json_schema={
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "point_index": {"type": "integer"},
                        "xi": {"type": "number"},
                        "n_analyzed": {"type": "integer"},
                        "time_start_fs": {"type": "number"},
                        "time_end_fs": {"type": "number"},
                        "time_total_fs": {"type": "number"},
                        "tau_corr": {"type": "number"},
                        "n_eff": {"type": "number"},
                        "sem_final_au": {"type": "number"},
                        "sem_max_au": {"type": ["number", "null"]},
                        "geweke_z": {"type": "number"},
                        "drift_D": {"type": "number"},
                        "passed": {"type": ["boolean", "null"]},
                        "failure_reasons": {
                            "type": "array", "items": {"type": "string"},
                        },
                    },
                },
            },
            type="list[dict]", shape="(N,)", category="metric",
        ),
    },
    outputs_raw_model={
        "ti_report": FieldSpec(
            description=(
                "Full TIReport object (models.py). Not JSON-serializable; "
                "MCP does NOT return this. Reserved for future Resources-layer "
                "field-path references (e.g. ti_report.point_reports[*])."
            ),
            json_schema={"$ref": "#/definitions/TIReport"},
            type="TIReport", category="raw_model",
        ),
    },
    preconditions=(
        "root_dir exists and is a directory",
        "Each discovered ti_target_*/ subdir contains a LagrangeMultLog file "
        "and a .restart (not a _N.restart checkpoint)",
        "dt is consistent across all loaded constraint points",
        "pattern ∈ {auto, ti_target, xi}",
        "If point_slice is provided: 2–3 colon-separated parts; integer "
        "components parseable; at least 2 points remain after slicing",
        (
            "If auto_equilibration=True, each point has >= "
            "DEFAULT_AUTO_EQUIL_MIN_FRAMES frames (currently 100, see "
            "constrained_ti/config.py)"
        ),
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Writes ti_convergence_report.csv, ti_free_energy.csv, "
        "ti_free_energy.png into output_dir (overwrites if present)",
        "Writes one diagnostics PNG per analyzed point into output_dir",
        "Does NOT mutate input TI directories",
    ),
    # NOTE: "API 成功但业务失败" (NOT_CONVERGED_* / CONVERGED_BUT_SHORT_TRAJ
    # 等) 不在 Tools 层。Resources 层基于 outputs_metrics.per_point 的信号
    # (passed / failure_reasons / time_total_fs / n_eff / drift_D / geweke_z)
    # 做判定。
    exceptions=(
        # Ordered: specific subclasses first.
        ExceptionMapping(
            exception_fqn=(
                "md_analysis.enhanced_sampling.constrained_ti.models"
                ".InsufficientSamplingError"
            ),
            triggered_by=(
                "auto_equilibration bisected below DEFAULT_AUTO_EQUIL_MIN_FRAMES"
            ),
            # Agent视角：输入数据不够 → validation，不是代码内部 analysis
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "root_dir missing, or a discovered ti_target_*/ subdir is "
                "missing required files"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid pattern / point_slice / inconsistent dt / fewer than "
                "2 points after slicing"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _ti_full_analysis_handler(params: dict[str, Any]) -> TaskResult:
    """Thin handler: delegate to run_ti_full_from_root; route artifacts + metrics."""
    from ..enhanced_sampling.constrained_ti.workflow import (
        run_ti_full_from_root,
    )

    report = run_ti_full_from_root(**params)

    # Artifacts → TaskResult.outputs
    outputs: dict[str, str] = {
        "convergence_csv": str(report.convergence_csv),
        "free_energy_csv": str(report.free_energy_csv),
        "free_energy_png": str(report.free_energy_png),
    }
    for i, p in enumerate(report.diagnostics_pngs):
        outputs[f"diagnostics_png_{i}"] = str(p)

    # Metrics → TaskResult.summary (JSON-serializable only)
    summary: dict[str, Any] = {
        "n_points": report.n_points,
        "delta_A_eV": report.delta_A_eV,
        "sigma_A_eV": report.sigma_A_eV,
        "all_passed": report.all_passed,
        "failing_indices": list(report.failing_indices),
        "per_point": list(report.per_point),
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
    handler=_ti_full_analysis_handler,
    target_fn=(
        "md_analysis.enhanced_sampling.constrained_ti.workflow"
        ":run_ti_full_from_root"
    ),
    cli_codes=("312",),
    contract=_TI_FULL_CONTRACT,
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


# 13. ti_gen_batch (CLI 422) — Batch-generate CP2K constrained-MD workdirs
#     for TI sampling points.  First task to ship with a full TaskContract;
#     contract is authoritative for schema + coercion + exception mapping.


def _ti_gen_batch_summary(
    result: Any, params: dict[str, Any],
) -> dict[str, Any]:
    """Summary extractor for TIGenBatchReport: promote metrics into summary."""
    if not hasattr(result, "workdirs"):
        return {}
    return {
        "n_targets": len(result.workdirs),
        "requested_targets_au": list(result.requested_targets_au),
        "snapped_targets_au": list(result.snapped_targets_au),
        "snap_deltas_au": list(result.snap_deltas_au),
        "steps": int(result.steps),
    }


_TI_GEN_BATCH_CONTRACT = TaskContract(
    inputs={
        "inp_path": FieldSpec(
            description="SG CP2K input file (filename may not end with .inp)",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "xyz_path": FieldSpec(
            description="SG trajectory XYZ file",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "restart_path": FieldSpec(
            description="SG .restart file (NOT a _N.restart checkpoint)",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "output_dir": FieldSpec(
            description="Parent directory to create ti_target_<cv>/ under",
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "targets_au": FieldSpec(
            description=(
                "CV target values in atomic units (numeric mode; mutually "
                "exclusive with time_range)"
            ),
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]",
            unit="a.u.", shape="(N,)",
            required=False, default=None,
        ),
        "time_range": FieldSpec(
            description=(
                "Time-mode target specification (mutually exclusive with "
                "targets_au): {time_initial_fs, time_final_fs, n_points}"
            ),
            json_schema={
                "type": "object",
                "properties": {
                    "time_initial_fs": {"type": "number"},
                    "time_final_fs": {"type": "number"},
                    "n_points": {"type": "integer", "minimum": 2},
                },
                "required": ["time_initial_fs", "time_final_fs", "n_points"],
                "additionalProperties": False,
            },
            type="dict",
            required=False, default=None,
        ),
        "steps": FieldSpec(
            description="MD steps written into each cMD.inp",
            json_schema={"type": "integer", "minimum": 1},
            type="int",
            default=10000, required=False,
        ),
        "script_path": FieldSpec(
            description=(
                "Submission script to copy as script.sh (falls back to "
                "KEY_CP2K_SCRIPT_PATH user config)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
    },
    outputs_artifacts={
        "workdirs": FieldSpec(
            description="Created ti_target_<cv>/ directories",
            json_schema={"type": "array", "items": {"type": "string"}},
            type="list[Path]", path_kind="dir", shape="(N,)",
            category="artifact",
        ),
    },
    outputs_metrics={
        "n_targets": FieldSpec(
            description="Number of TI points generated",
            json_schema={"type": "integer"}, type="int", category="metric",
        ),
        "requested_targets_au": FieldSpec(
            description="CV values originally requested (pre-snap)",
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="a.u.", shape="(N,)", category="metric",
        ),
        "snapped_targets_au": FieldSpec(
            description="Actual CV values used (snapped to nearest SG frame)",
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="a.u.", shape="(N,)", category="metric",
        ),
        "snap_deltas_au": FieldSpec(
            description="|requested - snapped| per target (informational)",
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="a.u.", shape="(N,)", category="metric",
        ),
        "steps": FieldSpec(
            description="MD steps used per TI point",
            json_schema={"type": "integer"}, type="int", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "inp_path is a valid CP2K SG input (contains PROJECT / TARGET / &TOPOLOGY)",
        "xyz_path frames carry atoms.info['i'] step index",
        "restart_path is a .restart file (not a _N.restart checkpoint)",
        "Exactly one of targets_au / time_range is provided",
        "If time_range: time_initial_fs <= time_final_fs AND n_points >= 2",
        "output_dir does not already contain ti_target_<cv>/ dirs matching "
        "any of the planned (snapped) targets",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Creates one ti_target_<cv>/ subdir per snapped target, containing "
        "init.xyz + cMD.inp [+ script.sh if script_path]",
        "Raises ValueError BEFORE any filesystem write if a planned target "
        "dirname collides with an existing ti_target_<cv>/ directory",
        "Does NOT mutate the input SG directory",
    ),
    exceptions=(
        # Order: specific subclasses first; then FileNotFoundError as generic
        # file-missing; then ValueError (incl. collision); then MDAnalysisError
        # fallback.
        ExceptionMapping(
            exception_fqn="md_analysis.scripts.TIGen.TIGenError",
            triggered_by=(
                "Invalid target specification (both modes or neither, bad "
                "time_range) or inp/restart parsing failure"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "inp_path / xyz_path / restart_path / script_path not found"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Collision: a planned target dirname already exists under "
                "output_dir (no writes performed)"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _make_ti_gen_batch_handler() -> Callable[[dict[str, Any]], TaskResult]:
    """Dedicated handler so we can plug in the summary extractor cleanly."""

    def handler(params: dict[str, Any]) -> TaskResult:
        from ..scripts.TIGen import generate_ti_batch_with_report

        result = generate_ti_batch_with_report(**params)
        return TaskResult(
            success=True,
            task="ti_gen_batch",
            outputs=_normalize_outputs(result),
            summary=_ti_gen_batch_summary(result, params),
        )

    return handler


register(TaskDef(
    name="ti_gen_batch",
    category="scripts",
    description=(
        "Batch-generate CP2K constrained-MD work directories for TI "
        "sampling points. Two target modes: numeric (explicit CV list in "
        "atomic units) or time (time_range object). Each target snaps to "
        "the nearest SG trajectory frame. Raises a validation error if "
        "any planned target directory already exists."
    ),
    handler=_make_ti_gen_batch_handler(),
    target_fn="md_analysis.scripts.TIGen:generate_ti_batch_with_report",
    cli_codes=("422",),
    contract=_TI_GEN_BATCH_CONTRACT,
))


# 14. bader_gen_batch (CLI 412) — Batch-prepare VASP Bader work directories.
#     Pass-through script-preparation contract example.
#     Does NOT submit jobs, parse Bader output, or compute surface charge.


def _bader_gen_batch_summary(
    result: Any, params: dict[str, Any],
) -> dict[str, Any]:
    """Summary extractor for BaderGenBatchReport."""
    if not hasattr(result, "workdirs"):
        return {}
    return {
        "n_frames": int(result.n_frames),
        "frame_indices": list(result.frame_indices),
        "steps": list(result.steps),
        "times_fs": list(result.times_fs),
        "generate_potcar": bool(result.generate_potcar),
    }


_BADER_GEN_BATCH_CONTRACT = TaskContract(
    inputs={
        "xyz_path": FieldSpec(
            description="CP2K XYZ trajectory file",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "cell_abc": FieldSpec(
            description="Orthogonal cell lengths [a, b, c]",
            json_schema={
                "type": "array",
                "items": {"type": "number", "exclusiveMinimum": 0},
                "minItems": 3, "maxItems": 3,
            },
            type="list[float]",
            unit="Angstrom", shape="(3,)",
        ),
        "output_dir": FieldSpec(
            description=(
                "Parent directory under which bader_t{time}_i{step}/ "
                "subdirs are created"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "mode": FieldSpec(
            description="Frame selection mode",
            json_schema={"type": "string", "enum": ["index", "time"]},
            type="str",
            choices=("index", "time"),
            required=False, default="index",
        ),
        "frame_start": FieldSpec(
            description="Index mode: first frame index (0-based, inclusive)",
            json_schema={"type": "integer", "minimum": 0},
            type="int",
            required=False, default=0,
        ),
        "frame_end": FieldSpec(
            description=(
                "Index mode: exclusive upper bound (None → all frames)"
            ),
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_step": FieldSpec(
            description="Index mode: step between frames",
            json_schema={"type": "integer", "minimum": 1},
            type="int",
            required=False, default=1,
        ),
        "time_start_fs": FieldSpec(
            description=(
                "Time mode: start time in fs (inclusive). Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "time_end_fs": FieldSpec(
            description=(
                "Time mode: end time in fs (inclusive). Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "time_step_fs": FieldSpec(
            description=(
                "Time mode: time step in fs (greedy matching). Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"], "exclusiveMinimum": 0},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "script_path": FieldSpec(
            description=(
                "Submission script to copy as script.sh (falls back to "
                "KEY_VASP_SCRIPT_PATH user config if not provided)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "element_order": FieldSpec(
            description=(
                "Optional element grouping order for POSCAR (e.g. "
                "['Cu','Ag','O','H']). Affects IndexMap on the comment line."
            ),
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "string"},
            },
            type="list[str] | None",
            required=False, default=None,
        ),
        "generate_potcar": FieldSpec(
            description=(
                "If True, invoke external `vaspkit 103` in each workdir to "
                "generate POTCAR. This is a local side effect; may fail if "
                "vaspkit is not on PATH."
            ),
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=True,
        ),
        "direct": FieldSpec(
            description="If True, POSCAR uses fractional (Direct) coordinates",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=True,
        ),
        "verbose": FieldSpec(
            description="If True, show tqdm progress bar during generation",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
    },
    outputs_artifacts={
        "workdirs": FieldSpec(
            description="Created bader_t{time}_i{step}/ directories",
            json_schema={"type": "array", "items": {"type": "string"}},
            type="list[Path]", path_kind="dir", shape="(N,)",
            category="artifact",
        ),
    },
    outputs_metrics={
        "n_frames": FieldSpec(
            description="Number of work directories created",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "frame_indices": FieldSpec(
            description="0-based trajectory frame indices that were used",
            json_schema={"type": "array", "items": {"type": "integer"}},
            type="list[int]", shape="(N,)", category="metric",
        ),
        "steps": FieldSpec(
            description="MD step values from atoms.info['i'] per frame",
            json_schema={"type": "array", "items": {"type": "integer"}},
            type="list[int]", shape="(N,)", category="metric",
        ),
        "times_fs": FieldSpec(
            description="MD time values from atoms.info['time'] per frame",
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="fs", shape="(N,)",
            category="metric",
        ),
        "generate_potcar": FieldSpec(
            description="Whether POTCAR generation (vaspkit 103) was invoked",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "xyz_path exists and is a CP2K XYZ trajectory file",
        "xyz_path frames carry atoms.info['i'] (step) and atoms.info['time'] (fs)",
        "cell_abc has length 3 with positive components",
        "mode ∈ {'index', 'time'}",
        "If mode='time': time_start_fs / time_end_fs / time_step_fs all provided",
        "If script_path is provided: exists and is a file",
        "If generate_potcar=True: `vaspkit` is on PATH",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Creates one bader_t{time}_i{step}/ subdir per selected frame",
        "Writes POSCAR, INCAR, KPOINTS per workdir",
        "Copies script.sh if script_path (or KEY_VASP_SCRIPT_PATH) is set",
        "If generate_potcar=True: invokes `vaspkit 103` in each workdir "
        "(local side effect; may be slow or fail if vaspkit is unavailable)",
        "Does NOT submit VASP jobs",
        "Does NOT parse Bader output or compute surface charge",
    ),
    exceptions=(
        # Ordered: specific subclasses first.
        ExceptionMapping(
            exception_fqn="md_analysis.scripts.BaderGen.BaderGenError",
            triggered_by=(
                "vaspkit missing from PATH, or vaspkit non-zero exit / no "
                "POTCAR generated"
            ),
            # Classified as ``validation`` because the dominant failure
            # mode is an unmet environment precondition (``vaspkit`` not on
            # PATH); the agent's retry policy for that is fix-inputs,
            # which matches ``validation`` rather than ``analysis``.
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn=(
                "md_analysis.scripts._frame_selector.FrameSelectionError"
            ),
            triggered_by=(
                "Time-mode missing any of time_start_fs / time_end_fs / "
                "time_step_fs; atoms.info missing 'time' metadata; "
                "invalid index-mode slice"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by="xyz_path or script_path missing on disk",
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by="cell_abc wrong length, or malformed numeric input",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Type mismatch (e.g. list where scalar expected)",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for other MDAnalysisError subclasses (e.g. "
                "IndexMapper parse issues)"
            ),
            error_type="analysis",
        ),
    ),
)


def _make_bader_gen_batch_handler() -> Callable[[dict[str, Any]], TaskResult]:
    """Handler for bader_gen_batch: route artifacts + metrics."""

    def handler(params: dict[str, Any]) -> TaskResult:
        from ..scripts.BaderGen import generate_bader_batch_with_report

        result = generate_bader_batch_with_report(**params)
        return TaskResult(
            success=True,
            task="bader_gen_batch",
            outputs=_normalize_outputs(result),
            summary=_bader_gen_batch_summary(result, params),
        )

    return handler


register(TaskDef(
    name="bader_gen_batch",
    category="scripts",
    description=(
        "Batch-prepare VASP Bader work directories from a CP2K MD "
        "trajectory. Selects frames by index or time, writes POSCAR/INCAR/"
        "KPOINTS per workdir, optionally copies a submission script, and "
        "optionally invokes vaspkit 103 for POTCAR. Does NOT submit jobs "
        "or parse Bader output; that stage is external."
    ),
    handler=_make_bader_gen_batch_handler(),
    target_fn=(
        "md_analysis.scripts.BaderGen:generate_bader_batch_with_report"
    ),
    cli_codes=("412",),
    contract=_BADER_GEN_BATCH_CONTRACT,
))
