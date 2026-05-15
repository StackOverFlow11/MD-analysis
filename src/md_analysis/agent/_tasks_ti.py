"""Constrained-TI agent tasks.

Home for ``ti_full_analysis`` today; future home for ``ti_bader_status``
and ``ti_constant_potential_correction``.
"""

from __future__ import annotations

from typing import Any

from ._contracts import ExceptionMapping, FieldSpec, TaskContract
from ._core import TaskDef, TaskResult, register

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
        "parser": FieldSpec(
            description=(
                "Engine parser name. 'auto' sniffs registered parsers "
                "against the first matching subdirectory. Currently "
                "registered: 'cp2k'."
            ),
            json_schema={"type": "string"},
            type="str",
            required=False, default="auto",
        ),
        "dir_filter": FieldSpec(
            description=(
                "Optional glob pattern to restrict candidate subdirectories "
                "(e.g. 'ti_target_*'). null means: any subdirectory the "
                "parser recognises by content."
            ),
            json_schema={"type": ["string", "null"]},
            type="str | None",
            required=False, default=None,
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
        "strict": FieldSpec(
            description=(
                "If True (default, agent-facing), a matched constraint-"
                "point directory missing required files raises "
                "FileNotFoundError. If False, such directories are "
                "skipped with a warning (the CLI menu path uses False)."
            ),
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=True,
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
        "At least one subdirectory (or one matching dir_filter) is recognised "
        "by the resolved parser (CP2K: contains *.restart and *.LagrangeMultLog)",
        "dt is consistent across all loaded constraint points",
        "parser is 'auto' or names a registered parser",
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
            exception_fqn=(
                "md_analysis.engines.protocols.ParserInferenceError"
            ),
            triggered_by=(
                "parser='auto' but no registered parser recognises any "
                "subdirectory under root_dir"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "root_dir missing, or a candidate subdir fails to parse "
                "required files (strict mode)"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid point_slice / inconsistent dt / fewer than "
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
    """Thin handler: dispatch via workflows facade; route artifacts + metrics.

    Phase 6.5: routes through workflows.enhanced_sampling.run_ti_full_analysis
    instead of run_ti_full_from_root directly. ``result.extra`` is the
    same TIFullAnalysisReport the previous wrapper returned, so the
    TaskResult outputs / summary contract is byte-equal — agent callers
    cannot observe the routing change.
    """
    from ..workflows.enhanced_sampling import run_ti_full_analysis

    result = run_ti_full_analysis(**params)
    report = result.extra  # TIFullAnalysisReport

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
    target_fn="md_analysis.workflows.enhanced_sampling:run_ti_full_analysis",
    cli_codes=("312",),
    contract=_TI_FULL_CONTRACT,
))
