"""Legacy (schema-derived) task registrations.

Most tasks here rely on ``target_fn`` signature introspection rather than
an explicit :class:`TaskContract`.  Registration order below is
load-bearing for ``list_tasks()`` output; contract-backed tasks from
Batch 1 onward remain in this module **only when** their CLI position
anchors the order (e.g. ``slowgrowth_quick`` #9, ``config_show`` #10);
larger scientific tasks live in sibling modules.
"""

from __future__ import annotations

from typing import Any

from ._contracts import ExceptionMapping, FieldSpec, TaskContract
from ._core import TaskDef, TaskResult, register
from ._handler_utils import _make_handler

# 7. calibration_fit_csv (CLI 231) — contract-backed.  Task name kept for
#    stability; the function also accepts manual ``data_points`` input.


_CALIBRATION_FIT_CONTRACT = TaskContract(
    inputs={
        "csv_path": FieldSpec(
            description=(
                "Two-column calibration CSV; column 1 = potential "
                "V vs SHE, column 2 = surface charge density uC/cm^2. "
                "Mutually exclusive with data_points."
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "data_points": FieldSpec(
            description=(
                "Manual (phi_V_vs_SHE, sigma_uC_cm2) points; mutually "
                "exclusive with csv_path"
            ),
            json_schema={
                "type": ["array", "null"],
                "items": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "maxItems": 2,
                },
            },
            type="list[tuple[float, float]] | None", shape="(N, 2)",
            required=False, default=None,
        ),
        "method": FieldSpec(
            description="Fitting method",
            json_schema={
                "type": "string",
                "enum": [
                    "linear", "polynomial", "spline",
                    "differential_capacitance",
                ],
            },
            type="str",
            choices=(
                "linear", "polynomial", "spline",
                "differential_capacitance",
            ),
            required=False, default="linear",
        ),
        "poly_degree": FieldSpec(
            description=(
                "Polynomial degree; used only when method='polynomial'"
            ),
            json_schema={"type": "integer", "minimum": 1},
            type="int",
            required=False, default=2,
        ),
        "output_dir": FieldSpec(
            description=(
                "Optional directory for calibration_data.csv and "
                "calibration_fit.png; skipped when None"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="dir",
            required=False, default=None,
        ),
        "calibration_json_path": FieldSpec(
            description=(
                "Destination JSON path for the calibration.  The "
                "underlying non-agent workflow defaults to "
                "~/.config/md_analysis/calibration.json; for agent "
                "calls, pass that path explicitly if the global "
                "default is desired."
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
    },
    outputs_artifacts={
        "calibration_json": FieldSpec(
            description="Saved calibration JSON (always present)",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "calibration_csv": FieldSpec(
            description=(
                "calibration_data.csv in output_dir (only when "
                "output_dir is provided and file exists)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "calibration_png": FieldSpec(
            description=(
                "calibration_fit.png in output_dir (only when "
                "output_dir is provided and file exists)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
    },
    outputs_metrics={
        "n_points": FieldSpec(
            description="Number of calibration data points used",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "reference": FieldSpec(
            description="Reference scale of stored data (currently 'SHE')",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "method": FieldSpec(
            description="Fitting method actually used",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "r_squared": FieldSpec(
            description=(
                "Coefficient of determination (reported fact only; the "
                "Tools layer does not judge calibration quality)"
            ),
            json_schema={"type": "number"},
            type="float", category="metric",
        ),
        "rmse": FieldSpec(
            description="Root-mean-square error (reported fact only)",
            json_schema={"type": "number"},
            type="float", category="metric",
        ),
        "equation": FieldSpec(
            description="Human-readable equation string",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "fit_params": FieldSpec(
            description=(
                "Full persisted fit-parameter dict (method-specific "
                "coefficients + r_squared + rmse + equation)"
            ),
            json_schema={"type": "object"},
            type="dict", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "Exactly one of csv_path / data_points is provided",
        "calibration_json_path is explicitly provided",
        "At least 2 calibration points available",
        "Each data point is (phi_V_vs_SHE, sigma_uC_cm2)",
        "method ∈ {'linear', 'polynomial', 'spline', 'differential_capacitance'}",
        "poly_degree >= 1",
        "For method='spline': scipy is importable",
        (
            "For method='differential_capacitance': σ is monotonically "
            "increasing with φ after sorting by φ"
        ),
    ),
    side_effects=(
        "Reads csv_path if provided",
        "Creates parent directory for calibration_json_path",
        "Writes calibration JSON (overwrites if present)",
        "Creates output_dir if provided",
        "Writes calibration_data.csv and calibration_fit.png in output_dir "
        "if provided (overwrites if present)",
        "Does NOT modify input CSV",
        "Does NOT perform charge analysis or constant-potential correction",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by="Missing csv_path on disk",
            error_type="file_not_found",
        ),
        # json.JSONDecodeError is a subclass of ValueError — list first.
        ExceptionMapping(
            exception_fqn="json.decoder.JSONDecodeError",
            triggered_by=(
                "Saved or loaded calibration JSON is not valid JSON"
            ),
            error_type="validation",
        ),
        # numpy.linalg.LinAlgError is a subclass of ValueError — list first.
        ExceptionMapping(
            exception_fqn="numpy.linalg.LinAlgError",
            triggered_by="Fit cannot be solved for the supplied data",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Both or neither of csv_path / data_points provided; "
                "missing calibration_json_path; malformed CSV rows; "
                "too few points; unknown method; invalid differential-"
                "capacitance data"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Malformed data_points shape or type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ImportError",
            triggered_by="method='spline' but scipy is unavailable",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not "
                "matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _handle_calibration_fit_csv(params: dict[str, Any]) -> TaskResult:
    """Dedicated handler — route artifacts + metrics from the WorkflowResult.

    Phase 6.4: dispatches through ``workflows.calibration.run_calibration_fit``
    instead of the business layer directly.  The TaskResult contract
    (outputs / summary) is byte-equal to the pre-migration baseline so
    agent callers cannot observe the routing change.
    """
    from ..workflows.calibration import run_calibration_fit

    result = run_calibration_fit(**params)
    report = result.extra  # CalibrationFitReport

    outputs: dict[str, str] = {
        "calibration_json": str(result.artifacts["calibration_json"]),
    }
    if "calibration_csv" in result.artifacts:
        outputs["calibration_csv"] = str(result.artifacts["calibration_csv"])
    if "calibration_png" in result.artifacts:
        outputs["calibration_png"] = str(result.artifacts["calibration_png"])

    summary: dict[str, Any] = {
        "n_points": int(report.n_points),
        "reference": str(report.reference),
        "method": str(report.method),
        "r_squared": float(report.r_squared),
        "rmse": float(report.rmse),
        "equation": str(report.equation),
        "fit_params": dict(report.fit_params),
    }
    return TaskResult(
        success=True,
        task="calibration_fit_csv",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="calibration_fit_csv",
    category="calibration",
    description=(
        "Fit a sigma→phi calibration curve from a CSV file or manual "
        "data_points (exactly one); writes a calibration JSON plus "
        "optional CSV/PNG artifacts when output_dir is given.  Task "
        "name retained for stability even though manual data_points "
        "input is also supported."
    ),
    handler=_handle_calibration_fit_csv,
    target_fn="md_analysis.workflows.calibration:run_calibration_fit",
    cli_codes=("231",),
    contract=_CALIBRATION_FIT_CONTRACT,
))


# 8. calibration_predict (CLI 233) — contract-backed.


_CALIBRATION_PREDICT_CONTRACT = TaskContract(
    inputs={
        "sigma": FieldSpec(
            description=(
                "Surface charge density; scalar or array (uC/cm^2)"
            ),
            json_schema={
                "oneOf": [
                    {"type": "number"},
                    {"type": "array", "items": {"type": "number"}},
                ],
            },
            type="float | list[float]",
            unit="uC/cm^2",
        ),
        "calibration_json_path": FieldSpec(
            description=(
                "Calibration JSON path.  The underlying non-agent "
                "workflow defaults to "
                "~/.config/md_analysis/calibration.json; for agent "
                "calls, pass that path explicitly if the global "
                "default is desired."
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "target_reference": FieldSpec(
            description=(
                "Output potential reference; None uses the reference "
                "stored in the calibration JSON"
            ),
            json_schema={
                "type": ["string", "null"],
                "enum": ["SHE", "RHE", "PZC", None],
            },
            type="str | None",
            choices=("SHE", "RHE", "PZC"),
            required=False, default=None,
        ),
        "temperature_K": FieldSpec(
            description="Temperature for RHE conversion",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="K",
            required=False, default=298.15,
        ),
        "pH": FieldSpec(
            description="pH for RHE conversion",
            json_schema={"type": "number"},
            type="float",
            required=False, default=0.0,
        ),
        "phi_pzc": FieldSpec(
            description=(
                "Potential of zero charge in V vs SHE; required for "
                "PZC conversions"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="V vs SHE",
            required=False, default=None,
        ),
    },
    outputs_artifacts={},
    outputs_metrics={
        "sigma_uC_cm2": FieldSpec(
            description=(
                "Input charge densities echoed as a list of floats "
                "(always a list, even for scalar input)"
            ),
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="uC/cm^2", shape="(N,)",
            category="metric",
        ),
        "potential_V": FieldSpec(
            description=(
                "Predicted electrode potentials (always a list, even "
                "for scalar input)"
            ),
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="V",
            shape="(N,)", category="metric",
        ),
        "n_values": FieldSpec(
            description="Number of predictions",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "is_scalar_input": FieldSpec(
            description="Whether the original sigma was a scalar",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "stored_reference": FieldSpec(
            description="Reference scale stored in the calibration JSON",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "target_reference": FieldSpec(
            description=(
                "Resolved reference actually used for the prediction "
                "(= stored_reference when caller passed None)"
            ),
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "temperature_K": FieldSpec(
            description="Temperature used for RHE conversion",
            json_schema={"type": "number"},
            type="float", unit="K", category="metric",
        ),
        "pH": FieldSpec(
            description="pH used for RHE conversion",
            json_schema={"type": "number"},
            type="float", category="metric",
        ),
        "phi_pzc": FieldSpec(
            description="Potential of zero charge, if supplied",
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="V vs SHE", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "calibration_json_path exists with loadable data + fit sections",
        "target_reference ∈ {None, 'SHE', 'RHE', 'PZC'}",
        "For PZC conversions (from or to): phi_pzc is provided",
        "temperature_K > 0 (documented; not runtime-checked in this batch)",
        "sigma contains numeric values",
    ),
    side_effects=(
        "Reads the calibration JSON",
        "Does NOT write files",
        "Does NOT modify calibration data",
        "Does NOT run charge analysis or constant-potential correction",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by="Missing calibration JSON",
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="json.decoder.JSONDecodeError",
            triggered_by="Malformed calibration JSON syntax",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Missing calibration_json_path; unknown method or "
                "reference; missing phi_pzc for PZC conversion; non-"
                "numeric sigma"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.KeyError",
            triggered_by=(
                "Missing required JSON sections / fit parameters"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ImportError",
            triggered_by=(
                "Loading a spline calibration but scipy is unavailable"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not "
                "matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _handle_calibration_predict(params: dict[str, Any]) -> TaskResult:
    """Dedicated handler — route JSON-serializable metrics from the WorkflowResult.

    Phase 6.4: dispatches through ``workflows.calibration.run_calibration_predict``
    instead of the business layer directly.  ``sigma`` is forwarded via
    ``**params`` because the workflow signature is positional-or-keyword
    on that argument; this preserves the existing splat style and avoids
    mutating the caller's params dict.
    """
    from ..workflows.calibration import run_calibration_predict

    result = run_calibration_predict(**params)
    report = result.extra  # CalibrationPredictionReport

    summary: dict[str, Any] = {
        "sigma_uC_cm2": [float(x) for x in report.sigma_uC_cm2],
        "potential_V": [float(x) for x in report.potential_V],
        "n_values": int(report.n_values),
        "is_scalar_input": bool(report.is_scalar_input),
        "stored_reference": str(report.stored_reference),
        "target_reference": str(report.target_reference),
        "temperature_K": float(report.temperature_K),
        "pH": float(report.pH),
        "phi_pzc": (
            float(report.phi_pzc) if report.phi_pzc is not None else None
        ),
    }
    return TaskResult(
        success=True,
        task="calibration_predict",
        outputs={},
        summary=summary,
    )


register(TaskDef(
    name="calibration_predict",
    category="calibration",
    description=(
        "Predict electrode potential from surface charge density using "
        "a saved calibration JSON.  Accepts scalar or array sigma; "
        "returns lists for both in every case.  Optional reference "
        "conversion to SHE / RHE / PZC."
    ),
    handler=_handle_calibration_predict,
    target_fn="md_analysis.workflows.calibration:run_calibration_predict",
    cli_codes=("233",),
    contract=_CALIBRATION_PREDICT_CONTRACT,
))

# 9. slowgrowth_quick (CLI 301) — contract-backed.
#     Boundary tightening: plot_style is now validated in the wrapper
#     (legacy slowgrowth_analysis silently skipped plots on unknown
#     values).  Returns structured summary metrics.


_SLOWGROWTH_QUICK_CONTRACT = TaskContract(
    inputs={
        "restart_path": FieldSpec(
            description="CP2K COLVAR restart file",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "log_path": FieldSpec(
            description="LagrangeMultLog file path",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "initial_step": FieldSpec(
            description=(
                "0-based initial array index of the slow-growth segment. "
                "If > final_step, the segment is reversed."
            ),
            json_schema={"type": "integer", "minimum": 0},
            type="int",
            required=False, default=0,
        ),
        "final_step": FieldSpec(
            description=(
                "Final (exclusive) array index; None → full length"
            ),
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "output_dir": FieldSpec(
            description=(
                "Output directory for CSV and PNG artifacts; defaults to "
                "the current working directory when None"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="dir",
            required=False, default=None,
        ),
        "plot_style": FieldSpec(
            description="Which plot(s) to produce",
            json_schema={
                "type": "string",
                "enum": ["quick", "publication", "both"],
            },
            type="str",
            choices=("quick", "publication", "both"),
            required=False, default="both",
        ),
        "colvar_id": FieldSpec(
            description=(
                "CP2K collective-variable ID; None uses the primary CV"
            ),
            json_schema={"type": ["integer", "null"]},
            type="int | None",
            required=False, default=None,
        ),
    },
    outputs_artifacts={
        "csv": FieldSpec(
            description="slowgrowth_data.csv file",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "quick_png": FieldSpec(
            description=(
                "slowgrowth_quick.png — present only for plot_style "
                "'quick' or 'both'"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "publication_png": FieldSpec(
            description=(
                "slowgrowth_publication.png — present only for plot_style "
                "'publication' or 'both'"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
    },
    outputs_metrics={
        "n_steps": FieldSpec(
            description="Length of the selected slow-growth segment",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "target_start_au": FieldSpec(
            description="CV target value at the segment start",
            json_schema={"type": "number"},
            type="float", unit="a.u.", category="metric",
        ),
        "target_end_au": FieldSpec(
            description="CV target value at the segment end",
            json_schema={"type": "number"},
            type="float", unit="a.u.", category="metric",
        ),
        "delta_F_eV": FieldSpec(
            description=(
                "Total free-energy change, "
                "free_energy_ev[-1] - free_energy_ev[0]"
            ),
            json_schema={"type": "number"},
            type="float", unit="eV", category="metric",
        ),
        "delta_F_barrier_eV": FieldSpec(
            description=(
                "Barrier height, max(free_energy_ev) - free_energy_ev[0]"
            ),
            json_schema={"type": "number"},
            type="float", unit="eV", category="metric",
        ),
        "barrier_step": FieldSpec(
            description=(
                "Original MD step at the barrier point (absolute pre-"
                "reversal step number, NOT the reset seg.steps index)"
            ),
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "is_reversed": FieldSpec(
            description="Whether initial_step > final_step triggered reversal",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "restart_path is a valid CP2K COLVAR restart file",
        "log_path is a valid LagrangeMultLog file",
        "plot_style ∈ {'quick', 'publication', 'both'}",
        "Selected segment is non-empty (initial_step != resolved final_step)",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Writes slowgrowth_data.csv into output_dir",
        "Writes 0, 1, or 2 PNG figures depending on plot_style",
        "Uses matplotlib Agg backend",
        "Does NOT mutate input restart/log files",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn=(
                "md_analysis.utils.formats.cp2k.colvar.ColvarParseError"
            ),
            triggered_by="Malformed restart or LagrangeMultLog content",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by="Missing restart_path or log_path",
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.KeyError",
            triggered_by="Invalid colvar_id for this restart",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid plot_style, empty selected segment, or malformed "
                "numeric input"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Empty or out-of-range selected segment",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not "
                "matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _handle_slowgrowth_quick(params: dict[str, Any]) -> TaskResult:
    """Handler: validate plot_style early, then delegate to the wrapper."""
    from ..enhanced_sampling.slowgrowth.SlowGrowthPlot import (
        slowgrowth_analysis_with_report,
    )

    report = slowgrowth_analysis_with_report(**params)

    outputs: dict[str, str] = {k: str(v) for k, v in report.artifacts.items()}
    summary: dict[str, Any] = {
        "n_steps": int(report.n_steps),
        "target_start_au": float(report.target_start_au),
        "target_end_au": float(report.target_end_au),
        "delta_F_eV": float(report.delta_F_eV),
        "delta_F_barrier_eV": float(report.delta_F_barrier_eV),
        "barrier_step": int(report.barrier_step),
        "is_reversed": bool(report.is_reversed),
    }
    return TaskResult(
        success=True,
        task="slowgrowth_quick",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="slowgrowth_quick",
    category="enhanced_sampling",
    description="Quick slow-growth free energy integration plot",
    handler=_handle_slowgrowth_quick,
    target_fn=(
        "md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot"
        ":slowgrowth_analysis_with_report"
    ),
    cli_codes=("301",),
    contract=_SLOWGROWTH_QUICK_CONTRACT,
))


# 10. config_show (CLI 900) — contract-backed, read-only.
#     Small boundary fix: honour ``config_path`` (previously ignored even
#     though the legacy signature-derived schema exposed it).


_CONFIG_SHOW_CONTRACT = TaskContract(
    inputs={
        "config_path": FieldSpec(
            description=(
                "Optional explicit path to a JSON config file.  Defaults "
                "to ``~/.config/md_analysis/config.json`` when omitted."
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
    },
    outputs_artifacts={},
    outputs_metrics={
        "config": FieldSpec(
            description=(
                "Loaded configuration dictionary; empty dict when no "
                "config file exists at the resolved path."
            ),
            json_schema={"type": "object"},
            type="dict", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "If config_path is provided: path is readable JSON or does not exist",
    ),
    side_effects=(
        "Reads one config JSON file if present",
        "Returns an empty dict if the resolved config file does not exist",
        "Does NOT create, modify, or delete any config file",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="md_analysis.config.ConfigError",
            triggered_by="Malformed JSON or unreadable config file",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by="Malformed numeric input at the agent boundary",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input type at the agent boundary",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="md_analysis.exceptions.MDAnalysisError",
            triggered_by=(
                "Fallback for any other MDAnalysisError subclass not "
                "matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _handle_config_show(params: dict[str, Any]) -> TaskResult:
    """Read-only: return the loaded config dict via TaskResult.summary."""
    from ..config import load_config

    config = load_config(params.get("config_path"))
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
    contract=_CONFIG_SHOW_CONTRACT,
))
