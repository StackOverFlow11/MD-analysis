"""Script / preparation agent tasks (``scripts/`` module counterparts).

Registration order: ``sp_gen_batch`` → ``ti_gen_batch`` → ``bader_gen_batch``,
matching the legacy ``_handlers.py`` order preserved by ``list_tasks()``.
"""

from __future__ import annotations

from typing import Any, Callable

from ._contracts import ExceptionMapping, FieldSpec, TaskContract
from ._core import TaskDef, TaskResult, register


# 12. sp_gen_batch (CLI 442) — Batch-generate CP2K SP work directories for
#     DeePMD training.  Contract-backed; delegates to
#     ``generate_sp_batch_with_report`` for preflight + structured metrics.


def _sp_gen_batch_summary(
    result: Any, params: dict[str, Any],
) -> dict[str, Any]:
    """Summary extractor for SpGenBatchReport."""
    if not hasattr(result, "workdirs"):
        return {}
    return {
        "n_frames": int(result.n_frames),
        "frame_indices": list(result.frame_indices),
        "steps": list(result.steps),
        "times_fs": list(result.times_fs),
    }


_SP_GEN_BATCH_CONTRACT = TaskContract(
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
                "Parent directory under which sp_t{time}_i{step}/ subdirs "
                "are created (existing same-name subdirs are overwritten)"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "inp_template_path": FieldSpec(
            description=(
                "DP SP inp template path; falls back to "
                "KEY_DP_SP_INP_TEMPLATE_PATH user config if omitted"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
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
                "Time mode: start time in fs (inclusive).  Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "time_end_fs": FieldSpec(
            description=(
                "Time mode: end time in fs (inclusive).  Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "time_step_fs": FieldSpec(
            description=(
                "Time mode: time step in fs (greedy matching).  Required when mode='time'"
            ),
            json_schema={"type": ["number", "null"], "exclusiveMinimum": 0},
            type="float | None", unit="fs",
            required=False, default=None,
        ),
        "script_path": FieldSpec(
            description=(
                "Submission script to copy as script.sh (falls back to "
                "KEY_CP2K_SCRIPT_PATH user config if not provided)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
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
            description="Created sp_t{time}_i{step}/ directories",
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
            description=(
                "MD step values from atoms.info['i'] per frame (fallback: "
                "frame index)"
            ),
            json_schema={"type": "array", "items": {"type": "integer"}},
            type="list[int]", shape="(N,)", category="metric",
        ),
        "times_fs": FieldSpec(
            description=(
                "MD time values from atoms.info['time'] per frame "
                "(fallback: 0.0)"
            ),
            json_schema={"type": "array", "items": {"type": "number"}},
            type="list[float]", unit="fs", shape="(N,)",
            category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "xyz_path exists and is a CP2K XYZ trajectory file",
        "inp_template_path exists (or KEY_DP_SP_INP_TEMPLATE_PATH resolves)",
        "cell_abc has length 3 with positive components",
        "mode ∈ {'index', 'time'}",
        "If mode='time': time_start_fs / time_end_fs / time_step_fs all provided",
        "If script_path is provided (or KEY_CP2K_SCRIPT_PATH resolves): it exists",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Creates one sp_t{time}_i{step}/ subdir per selected frame",
        "Writes init.xyz and sp.inp per workdir",
        "Copies script.sh if script_path (or KEY_CP2K_SCRIPT_PATH) is set",
        "Overwrites init.xyz / sp.inp / script.sh inside same-name "
        "sp_t{time}_i{step}/ directories that already exist",
        "Reads user config when template or script paths are omitted",
        "Does NOT submit CP2K jobs",
    ),
    exceptions=(
        # Ordered: specific subclasses first.
        ExceptionMapping(
            exception_fqn="md_analysis.scripts.SpGen.SpGenError",
            triggered_by=(
                "No inp template provided and KEY_DP_SP_INP_TEMPLATE_PATH "
                "not configured"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn=(
                "md_analysis.scripts._frame_selector.FrameSelectionError"
            ),
            triggered_by=(
                "Invalid frame/time selection parameters; time-mode missing "
                "any of time_start_fs / time_end_fs / time_step_fs; "
                "atoms.info missing 'time' metadata in time mode"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "xyz_path, resolved inp_template_path, or resolved "
                "script_path missing on disk"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "cell_abc wrong length, non-positive cell component, or "
                "malformed numeric input"
            ),
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
                "Fallback for any other MDAnalysisError subclass not "
                "matched above"
            ),
            error_type="analysis",
        ),
    ),
)


def _make_sp_gen_batch_handler() -> Callable[[dict[str, Any]], TaskResult]:
    """Handler for sp_gen_batch: route artifacts + metrics.

    Phase 6.agent-cleanup: dispatches through
    workflows.scripts.run_sp_batch (mirrors the 6.6 ti_gen_batch
    reroute). WorkflowResult has no .workdirs attribute, so outputs
    are built explicitly from result.artifacts (NOT _normalize_outputs,
    which would yield an empty dict) and summary is read from
    result.extra (the SpGenBatchReport). outputs (workdir_*) / summary
    contract stays byte-equal.
    """

    def handler(params: dict[str, Any]) -> TaskResult:
        from ..workflows.scripts import run_sp_batch

        result = run_sp_batch(**params)
        return TaskResult(
            success=True,
            task="sp_gen_batch",
            outputs={k: str(v) for k, v in result.artifacts.items()},
            summary=_sp_gen_batch_summary(result.extra, params),
        )

    return handler


register(TaskDef(
    name="sp_gen_batch",
    category="scripts",
    description=(
        "Batch-generate CP2K SP work directories for DeePMD training. "
        "Supports two frame-selection modes: mode='index' uses "
        "frame_start/end/step (0-based indices); mode='time' uses "
        "time_start_fs/end_fs/step_fs (all three required together, "
        "greedy matching on atoms.info['time']).  Does NOT submit jobs."
    ),
    handler=_make_sp_gen_batch_handler(),
    target_fn="md_analysis.workflows.scripts:run_sp_batch",
    cli_codes=("442",),
    contract=_SP_GEN_BATCH_CONTRACT,
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
        "colvar_id": FieldSpec(
            description=(
                "Constraint colvar_id the CV target is read from. "
                "null (default) = primary / first colvar."
            ),
            json_schema={"type": ["integer", "null"]},
            type="int | None",
            required=False, default=None,
        ),
        "overwrite": FieldSpec(
            description=(
                "If True, skip the pre-write collision guard and "
                "regenerate each ti_target_<cv>/ directory in place, "
                "overwriting an existing cMD.inp / init.xyz / script.sh "
                "(per-file overwrite via mkdir(exist_ok=True) — NOT a "
                "directory wipe; other files in a colliding directory are "
                "left untouched). Default False keeps the collision guard "
                "(a planned directory that already exists raises a "
                "validation error). CLI 422 passes True."
            ),
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
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
        "Unless overwrite=True: output_dir does not already contain "
        "ti_target_<cv>/ dirs matching any of the planned (snapped) targets",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Creates one ti_target_<cv>/ subdir per snapped target, containing "
        "init.xyz + cMD.inp [+ script.sh if script_path]",
        "If overwrite=False: raises ValueError BEFORE any filesystem write "
        "when a planned target dirname collides with an existing "
        "ti_target_<cv>/ directory",
        "If overwrite=True: skips the collision guard and rewrites "
        "cMD.inp / init.xyz / script.sh in a colliding directory in place "
        "(per-file overwrite; no directory removal/wipe)",
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
                "output_dir and overwrite=False (no writes performed)"
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
    """Dedicated handler so we can plug in the summary extractor cleanly.

    Phase 6.6: dispatches through workflows.scripts.run_ti_batch. The
    returned WorkflowResult has NO ``.workdirs`` attribute, so
    ``_normalize_outputs`` would fall through to its empty-dict branch —
    outputs are built explicitly from ``result.artifacts`` instead
    (``_flatten_workdirs_to_artifacts`` produces ``workdir_0/1/...``
    keys, byte-equal to the pre-migration ``_normalize_outputs``
    shape). Summary is read from ``result.extra`` (the TIGenBatchReport).
    """

    def handler(params: dict[str, Any]) -> TaskResult:
        from ..workflows.scripts import run_ti_batch

        result = run_ti_batch(**params)
        return TaskResult(
            success=True,
            task="ti_gen_batch",
            outputs={k: str(v) for k, v in result.artifacts.items()},
            summary=_ti_gen_batch_summary(result.extra, params),
        )

    return handler


register(TaskDef(
    name="ti_gen_batch",
    category="scripts",
    description=(
        "Batch-generate CP2K constrained-MD work directories for TI "
        "sampling points. Two target modes: numeric (explicit CV list in "
        "atomic units) or time (time_range object). Each target snaps to "
        "the nearest SG trajectory frame. colvar_id selects the "
        "constraint (None = primary). Raises a validation error if any "
        "planned target directory already exists, unless overwrite=True."
    ),
    handler=_make_ti_gen_batch_handler(),
    target_fn="md_analysis.workflows.scripts:run_ti_batch",
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
    """Handler for bader_gen_batch: route artifacts + metrics.

    Phase 6.agent-cleanup: dispatches through
    workflows.scripts.run_bader_batch (mirrors the 6.6 ti_gen_batch
    reroute). outputs built explicitly from result.artifacts (NOT
    _normalize_outputs — WorkflowResult has no .workdirs), summary from
    result.extra (the BaderGenBatchReport). Contract byte-equal.
    """

    def handler(params: dict[str, Any]) -> TaskResult:
        from ..workflows.scripts import run_bader_batch

        result = run_bader_batch(**params)
        return TaskResult(
            success=True,
            task="bader_gen_batch",
            outputs={k: str(v) for k, v in result.artifacts.items()},
            summary=_bader_gen_batch_summary(result.extra, params),
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
    target_fn="md_analysis.workflows.scripts:run_bader_batch",
    cli_codes=("412",),
    contract=_BADER_GEN_BATCH_CONTRACT,
))
