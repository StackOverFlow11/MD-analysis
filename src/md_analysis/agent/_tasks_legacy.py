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

from ..electrochemical.potential.config import DEFAULT_THICKNESS_ANG
from ..utils.constants import DEFAULT_LAYER_TOL_A
from ._contracts import ExceptionMapping, FieldSpec, TaskContract
from ._core import TaskDef, TaskResult, register
from ._handler_utils import _make_handler

# 1. water_three_panel (CLI 105) — contract-backed (Batch 4).


_WATER_THREE_PANEL_CONTRACT = TaskContract(
    inputs={
        "xyz_path": FieldSpec(
            description="CP2K XYZ trajectory file",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "md_inp_path": FieldSpec(
            description=(
                "CP2K md.inp used to read cell dimensions when "
                "cell_abc is not provided"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "cell_abc": FieldSpec(
            description="Orthogonal cell lengths [a, b, c]",
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "number", "exclusiveMinimum": 0},
                "minItems": 3, "maxItems": 3,
            },
            type="list[float] | None",
            unit="Angstrom", shape="(3,)",
            required=False, default=None,
        ),
        "output_dir": FieldSpec(
            description=(
                "Output directory for water CSV / TXT / PNG artifacts"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "layer_tol_A": FieldSpec(
            description="Layer clustering tolerance",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=DEFAULT_LAYER_TOL_A,
        ),
        "frame_start": FieldSpec(
            description="Optional slice start (0-based, inclusive)",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_end": FieldSpec(
            description="Optional slice exclusive upper bound",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_step": FieldSpec(
            description="Optional slice step",
            json_schema={"type": ["integer", "null"], "minimum": 1},
            type="int | None",
            required=False, default=None,
        ),
        "verbose": FieldSpec(
            description="If True, show tqdm progress bar",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
    },
    outputs_artifacts={
        "density_csv": FieldSpec(
            description="Water mass density CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "orientation_csv": FieldSpec(
            description="Water orientation-weighted density CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_profile_csv": FieldSpec(
            description="Adsorbed water layer profile CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_range_txt": FieldSpec(
            description="Adsorbed layer range / peak metrics TXT",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_theta_csv": FieldSpec(
            description="Adsorbed water θ angular distribution CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "plot_png": FieldSpec(
            description="Three-panel summary PNG",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
    },
    outputs_metrics={
        "output_dir": FieldSpec(
            description="Echo of resolved output_dir",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "n_artifacts": FieldSpec(
            description="Number of artifact files produced",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "frame_start": FieldSpec(
            description="Echo of applied slice start (None = use all)",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
        "frame_end": FieldSpec(
            description="Echo of applied slice end (None = use all)",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
        "frame_step": FieldSpec(
            description="Echo of applied slice step (None = 1)",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "xyz_path exists and is readable",
        "Either cell_abc is provided or md_inp_path can yield cell dimensions",
        "If cell_abc is provided: length 3 with positive values",
        "Selected frame slice is non-empty",
    ),
    side_effects=(
        "Creates output_dir if missing",
        "Writes water CSV / TXT / PNG artifacts into output_dir",
        "May overwrite files with the same names",
        "Does NOT mutate input trajectory or input files",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing xyz_path; or md_inp_path missing when needed "
                "for cell derivation"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Invalid frame slice indexing",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid cell_abc shape / values; malformed trajectory; "
                "layer detection fails for input data; invalid frame slice"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.RuntimeError",
            triggered_by=(
                "Lower-level analysis succeeded but one or more expected "
                "artifacts were not produced"
            ),
            error_type="analysis",
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


def _handle_water_three_panel(params: dict[str, Any]) -> TaskResult:
    """Handler — route WaterThreePanelReport fields."""
    from ..main import run_water_analysis_with_report

    report = run_water_analysis_with_report(**params)
    outputs: dict[str, str] = {k: str(v) for k, v in report.artifacts.items()}
    summary: dict[str, Any] = {
        "output_dir": str(report.output_dir),
        "n_artifacts": int(report.n_artifacts),
        "frame_start": (
            int(report.frame_start) if report.frame_start is not None else None
        ),
        "frame_end": (
            int(report.frame_end) if report.frame_end is not None else None
        ),
        "frame_step": (
            int(report.frame_step) if report.frame_step is not None else None
        ),
    }
    return TaskResult(
        success=True,
        task="water_three_panel",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="water_three_panel",
    category="water",
    description=(
        "Water density / orientation / adsorbed-layer three-panel "
        "analysis.  Reads a CP2K XYZ trajectory plus cell (either via "
        "cell_abc or md_inp_path) and writes six artifacts (3 CSV + "
        "1 TXT + 1 CSV + 1 PNG).  Does NOT decide whether sampling is "
        "statistically sufficient."
    ),
    handler=_handle_water_three_panel,
    target_fn="md_analysis.main:run_water_analysis_with_report",
    cli_codes=("105",),
    contract=_WATER_THREE_PANEL_CONTRACT,
))


# 2. potential_full (CLI 216) — contract-backed (Batch 4).


_POTENTIAL_FULL_CONTRACT = TaskContract(
    inputs={
        "output_dir": FieldSpec(
            description=(
                "Output directory; per-sub-analysis artifacts are "
                "written into output_dir/<sub>/ (center / fermi / "
                "electrode / phi_z / thickness_sensitivity)"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "cube_pattern": FieldSpec(
            description=(
                "Glob pattern for continuous-mode cube files; resolved "
                "relative to the process current working directory"
            ),
            json_schema={"type": "string"},
            type="str",
            required=False, default="md-POTENTIAL-v_hartree-1_*.cube",
        ),
        "md_out_path": FieldSpec(
            description=(
                "CP2K md.out for Fermi energy; optional in distributed "
                "mode (sp.out per-subdir is used instead)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "xyz_path": FieldSpec(
            description=(
                "XYZ trajectory for interface-based slab centering; "
                "required when center_mode='interface'"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "thickness_ang": FieldSpec(
            description="Slab thickness for centering / averaging",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=DEFAULT_THICKNESS_ANG,
        ),
        "center_mode": FieldSpec(
            description="Slab centering strategy",
            json_schema={"type": "string", "enum": ["interface", "cell"]},
            type="str",
            choices=("interface", "cell"),
            required=False, default="interface",
        ),
        "metal_elements": FieldSpec(
            description=(
                "Override default metal element list for layer / "
                "interface detection"
            ),
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "string"},
            },
            type="list[str] | None",
            required=False, default=None,
        ),
        "layer_tol_ang": FieldSpec(
            description="Layer clustering tolerance",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=DEFAULT_LAYER_TOL_A,
        ),
        "fermi_unit": FieldSpec(
            description="Fermi energy output unit",
            json_schema={"type": "string", "enum": ["au", "ev"]},
            type="str",
            choices=("au", "ev"),
            required=False, default="au",
        ),
        "compute_u": FieldSpec(
            description=(
                "If True and Fermi data is available, compute electrode "
                "potential U vs SHE"
            ),
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=True,
        ),
        "compute_phi_z": FieldSpec(
            description=(
                "If True, produce the φ(z) planar-average overlay PNG"
            ),
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=True,
        ),
        "max_curves": FieldSpec(
            description="Max number of φ(z) curves to overlay (0 = all)",
            json_schema={"type": "integer", "minimum": 0},
            type="int",
            required=False, default=0,
        ),
        "thickness_end": FieldSpec(
            description="Upper bound for thickness sensitivity sweep",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=15.0,
        ),
        "frame_start": FieldSpec(
            description="Optional slice start (0-based, inclusive)",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_end": FieldSpec(
            description="Optional slice exclusive upper bound",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_step": FieldSpec(
            description="Optional slice step",
            json_schema={"type": ["integer", "null"], "minimum": 1},
            type="int | None",
            required=False, default=None,
        ),
        "verbose": FieldSpec(
            description="If True, show tqdm progress bar",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
        "input_mode": FieldSpec(
            description=(
                "'continuous' = read a MD cube pattern; 'distributed' "
                "= read per-frame SP subdirectories under sp_root_dir"
            ),
            json_schema={
                "type": "string",
                "enum": ["continuous", "distributed"],
            },
            type="str",
            choices=("continuous", "distributed"),
            required=False, default="continuous",
        ),
        "sp_root_dir": FieldSpec(
            description=(
                "Distributed-mode SP root directory.  Required by "
                "behaviour when input_mode='distributed'"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="dir",
            required=False, default=None,
        ),
        "sp_dir_pattern": FieldSpec(
            description="Glob for SP subdirectories under sp_root_dir",
            json_schema={"type": "string"},
            type="str",
            required=False, default="potential_t*_i*",
        ),
        "sp_cube_filename": FieldSpec(
            description="Per-SP-subdir cube filename",
            json_schema={"type": "string"},
            type="str",
            required=False,
            default="sp_potential-v_hartree-1_0.cube",
        ),
        "sp_out_filename": FieldSpec(
            description="Per-SP-subdir stdout filename (for Fermi energy)",
            json_schema={"type": "string"},
            type="str",
            required=False, default="sp.out",
        ),
    },
    outputs_artifacts={
        "electrode_csv": FieldSpec(
            description=(
                "Electrode potential U vs SHE CSV; produced only when "
                "compute_u=True AND Fermi energy is available"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "center_csv": FieldSpec(
            description=(
                "Slab-centred Hartree potential CSV; produced only "
                "when compute_u=False OR Fermi is unavailable (then "
                "the component sub-analyses run separately)"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "fermi_csv": FieldSpec(
            description=(
                "Fermi energy time-series CSV; produced only when "
                "compute_u=False AND md_out_path is provided"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "phi_z_png": FieldSpec(
            description=(
                "φ(z) planar-average overlay PNG; produced only when "
                "compute_phi_z=True"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "thickness_sensitivity_csv": FieldSpec(
            description=(
                "Thickness sensitivity sweep CSV; produced only when "
                "Fermi energy is available"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
    },
    outputs_metrics={
        "output_dir": FieldSpec(
            description="Echo of resolved output_dir",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "input_mode": FieldSpec(
            description="Echo of resolved input_mode",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "has_fermi": FieldSpec(
            description=(
                "True iff any Fermi-dependent artifact (electrode / "
                "fermi / thickness sensitivity) was produced"
            ),
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_electrode": FieldSpec(
            description="Whether electrode_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_center": FieldSpec(
            description="Whether center_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_fermi": FieldSpec(
            description="Whether fermi_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_phi_z": FieldSpec(
            description="Whether phi_z_png was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_thickness_sensitivity": FieldSpec(
            description="Whether thickness_sensitivity_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "n_artifacts": FieldSpec(
            description="Number of artifacts produced",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "input_mode ∈ {'continuous', 'distributed'}",
        "continuous mode: cube_pattern matches at least one cube file "
        "relative to the current working directory",
        "If center_mode='interface': xyz_path is provided",
        "If compute_u=True in continuous mode: md_out_path is provided",
        "If input_mode='distributed': sp_root_dir exists with valid SP "
        "subdirectories containing sp_cube_filename",
        "fermi_unit ∈ {'au', 'ev'}",
        "thickness_ang / layer_tol_ang / thickness_end are positive",
    ),
    side_effects=(
        "Creates output_dir and per-sub-analysis subdirectories",
        "Writes CSV / PNG artifacts for sub-analyses that run",
        "May overwrite files with the same names",
        "Does NOT mutate input cube / md.out / xyz / SP files",
        "Does NOT parse generated CSVs for scientific metrics",
        "Does NOT judge whether the sampling is sufficient or the "
        "system is equilibrated",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing cube files; missing md_out_path / xyz_path / "
                "sp_root_dir / SP cube files; no valid frames after "
                "slicing"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Invalid frame slice indexing",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid input_mode / center_mode / fermi_unit; "
                "malformed cube/out/xyz input; invalid frame slicing"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.RuntimeError",
            triggered_by=(
                "No Fermi energy records or no matched frames despite "
                "apparently valid input"
            ),
            error_type="analysis",
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


def _handle_potential_full(params: dict[str, Any]) -> TaskResult:
    """Handler — route PotentialFullReport fields.

    Only artifacts that actually exist on disk make it into ``outputs``;
    the boolean ``ran_*`` summary flags mirror that existence check
    (not the caller's requested flags).
    """
    from ..main import run_potential_analysis_with_report

    report = run_potential_analysis_with_report(**params)
    outputs: dict[str, str] = {
        k: str(v) for k, v in report.artifacts.items()
    }
    summary: dict[str, Any] = {
        "output_dir": str(report.output_dir),
        "input_mode": str(report.input_mode),
        "has_fermi": bool(report.has_fermi),
        "ran_electrode": bool(report.ran_electrode),
        "ran_center": bool(report.ran_center),
        "ran_fermi": bool(report.ran_fermi),
        "ran_phi_z": bool(report.ran_phi_z),
        "ran_thickness_sensitivity": bool(report.ran_thickness_sensitivity),
        "n_artifacts": int(report.n_artifacts),
    }
    return TaskResult(
        success=True,
        task="potential_full",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="potential_full",
    category="potential",
    description=(
        "Full potential analysis bundle: slab-centred Hartree "
        "potential + Fermi energy + electrode potential U vs SHE + "
        "φ(z) profile + thickness sensitivity sweep.  Sub-analyses "
        "run conditionally on inputs (Fermi availability, compute_* "
        "flags, input_mode).  Only actually produced artifacts appear "
        "in TaskResult.outputs; summary booleans mirror that "
        "existence check, not the requested flags.  Does NOT judge "
        "electrochemical reliability or sampling sufficiency."
    ),
    handler=_handle_potential_full,
    target_fn="md_analysis.main:run_potential_analysis_with_report",
    cli_codes=("216",),
    contract=_POTENTIAL_FULL_CONTRACT,
))

# 3. charge_surface (CLI 221-223) — contract-backed.
#    Analyses existing Bader frame directories; does NOT generate Bader
#    dirs or submit jobs; does NOT decide whether constant-potential
#    correction is scientifically needed.


_CHARGE_COMMON_FRAME_FIELDS = {
    "dir_pattern": FieldSpec(
        description="Glob pattern for Bader frame subdirectories",
        json_schema={"type": "string"},
        type="str",
        required=False, default="bader_t*_i*",
    ),
    "frame_start": FieldSpec(
        description="Optional slice start (0-based, inclusive)",
        json_schema={"type": ["integer", "null"], "minimum": 0},
        type="int | None",
        required=False, default=None,
    ),
    "frame_end": FieldSpec(
        description="Optional slice exclusive upper bound",
        json_schema={"type": ["integer", "null"], "minimum": 0},
        type="int | None",
        required=False, default=None,
    ),
    "frame_step": FieldSpec(
        description="Optional slice step",
        json_schema={"type": ["integer", "null"], "minimum": 1},
        type="int | None",
        required=False, default=None,
    ),
    "verbose": FieldSpec(
        description="If True, show tqdm progress bar",
        json_schema={"type": "boolean"},
        type="bool",
        required=False, default=False,
    ),
}


_CHARGE_SURFACE_CONTRACT = TaskContract(
    inputs={
        "output_dir": FieldSpec(
            description=(
                "Parent output directory; actual outputs go into "
                "output_dir / method"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "root_dir": FieldSpec(
            description=(
                "Parent directory containing Bader frame subdirectories"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
            required=False, default=".",
        ),
        "metal_symbols": FieldSpec(
            description=(
                "Override default metal symbols for layer detection"
            ),
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "string"},
            },
            type="list[str] | None",
            required=False, default=None,
        ),
        "normal": FieldSpec(
            description="Cell axis perpendicular to the surface",
            json_schema={"type": "string", "enum": ["a", "b", "c"]},
            type="str",
            choices=("a", "b", "c"),
            required=False, default="c",
        ),
        "method": FieldSpec(
            description="Charge partitioning method",
            json_schema={"type": "string", "enum": ["counterion", "layer"]},
            type="str",
            choices=("counterion", "layer"),
            required=False, default="counterion",
        ),
        "layer_tol_A": FieldSpec(
            description="Layer clustering tolerance",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=DEFAULT_LAYER_TOL_A,
        ),
        "n_surface_layers": FieldSpec(
            description=(
                "Number of metal layers (per interface) summed for "
                "method='layer'"
            ),
            json_schema={"type": "integer", "minimum": 1},
            type="int",
            required=False, default=1,
        ),
        **_CHARGE_COMMON_FRAME_FIELDS,
    },
    outputs_artifacts={
        "charge_csv": FieldSpec(
            description="surface_charge.csv",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "charge_png": FieldSpec(
            description="surface_charge.png",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
    },
    outputs_metrics={
        "n_frames": FieldSpec(
            description="Number of frames analyzed",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "method": FieldSpec(
            description="Charge partitioning method used",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "sigma_aligned_mean": FieldSpec(
            description="Mean surface charge density on +axis side",
            json_schema={"type": "number"},
            type="float", unit="uC/cm^2", category="metric",
        ),
        "sigma_aligned_std": FieldSpec(
            description="Stdev of surface charge density on +axis side",
            json_schema={"type": "number"},
            type="float", unit="uC/cm^2", category="metric",
        ),
        "sigma_opposed_mean": FieldSpec(
            description="Mean surface charge density on −axis side",
            json_schema={"type": "number"},
            type="float", unit="uC/cm^2", category="metric",
        ),
        "sigma_opposed_std": FieldSpec(
            description="Stdev of surface charge density on −axis side",
            json_schema={"type": "number"},
            type="float", unit="uC/cm^2", category="metric",
        ),
        "phi_cumavg_last": FieldSpec(
            description=(
                "Last cumulative-average potential (only populated when "
                "single-side mode is used and a default-location "
                "calibration JSON exists)"
            ),
            json_schema={"type": ["number", "null"]},
            type="float | None", unit="V", category="metric",
        ),
        "phi_reference": FieldSpec(
            description=(
                "Reference scale used for extrapolated potential "
                "(only when calibration appended)"
            ),
            json_schema={"type": ["string", "null"]},
            type="str | None", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "root_dir exists and contains at least one matching frame subdir",
        "Each matching subdir contains POSCAR, ACF.dat, POTCAR",
        "normal ∈ {'a', 'b', 'c'}",
        "method ∈ {'counterion', 'layer'}",
        "n_surface_layers >= 1 (and <= total metal layers when method='layer')",
    ),
    side_effects=(
        "Creates output_dir / method",
        "Writes surface_charge.csv and surface_charge.png",
        "Reads the default-location calibration JSON if available to "
        "append phi_* columns (optional reporting behaviour, NOT a "
        "constant-potential correction)",
        "Does NOT mutate Bader frame directories",
        "Does NOT submit jobs or generate Bader inputs",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing root_dir, no frame dirs after slicing, or "
                "missing POSCAR / ACF.dat / POTCAR inside a frame"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Out-of-bounds atom or layer indexing",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid normal / method / n_surface_layers; malformed "
                "Bader content; invalid frame slicing"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
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


def _handle_charge_surface(params: dict[str, Any]) -> TaskResult:
    """Handler — route ChargeSurfaceReport fields."""
    from ..main import run_charge_analysis_with_report

    report = run_charge_analysis_with_report(**params)

    outputs: dict[str, str] = {
        "charge_csv": str(report.charge_csv),
        "charge_png": str(report.charge_png),
    }
    summary: dict[str, Any] = {
        "n_frames": int(report.n_frames),
        "method": str(report.method),
        "sigma_aligned_mean": float(report.sigma_aligned_mean),
        "sigma_aligned_std": float(report.sigma_aligned_std),
        "sigma_opposed_mean": float(report.sigma_opposed_mean),
        "sigma_opposed_std": float(report.sigma_opposed_std),
        "phi_cumavg_last": (
            float(report.phi_cumavg_last)
            if report.phi_cumavg_last is not None else None
        ),
        "phi_reference": (
            str(report.phi_reference)
            if report.phi_reference is not None else None
        ),
    }
    return TaskResult(
        success=True,
        task="charge_surface",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="charge_surface",
    category="charge",
    description=(
        "Bader surface charge density time series.  Analyses existing "
        "Bader frame directories (POSCAR + ACF.dat + POTCAR per frame); "
        "does NOT generate Bader dirs, submit VASP jobs, or decide "
        "whether constant-potential correction is scientifically "
        "needed.  When a calibration JSON is available at the default "
        "location, extrapolated φ columns are appended to the CSV as "
        "an optional reporting behaviour."
    ),
    handler=_handle_charge_surface,
    target_fn="md_analysis.main:run_charge_analysis_with_report",
    cli_codes=("221", "222", "223"),
    contract=_CHARGE_SURFACE_CONTRACT,
))


# 4. charge_tracked (CLI 225) — contract-backed.


_CHARGE_TRACKED_CONTRACT = TaskContract(
    inputs={
        "output_dir": FieldSpec(
            description=(
                "Parent output directory; actual outputs go into "
                "output_dir / 'tracked'"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "root_dir": FieldSpec(
            description=(
                "Parent directory containing Bader frame subdirectories"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
            required=False, default=".",
        ),
        "atom_indices_xyz": FieldSpec(
            description=(
                "0-based atom indices in ORIGINAL XYZ ordering (not "
                "POSCAR order); at least one"
            ),
            json_schema={
                "type": "array",
                "items": {"type": "integer", "minimum": 0},
                "minItems": 1,
            },
            type="list[int]",
        ),
        **_CHARGE_COMMON_FRAME_FIELDS,
    },
    outputs_artifacts={
        "tracked_charge_csv": FieldSpec(
            description="tracked_atom_charges.csv",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "tracked_charge_png": FieldSpec(
            description="tracked_atom_charges.png",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
    },
    outputs_metrics={
        "n_frames": FieldSpec(
            description="Number of frames analyzed",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "atom_indices_xyz": FieldSpec(
            description="Echo of the tracked indices (XYZ order)",
            json_schema={"type": "array", "items": {"type": "integer"}},
            type="list[int]", shape="(N,)", category="metric",
        ),
        "n_atoms_tracked": FieldSpec(
            description="Number of atoms tracked",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "root_dir exists and contains at least one matching frame subdir",
        "Each matching subdir contains POSCAR, ACF.dat, POTCAR",
        "atom_indices_xyz is non-empty and all indices are >= 0",
        "All indices are < total atoms in each frame",
    ),
    side_effects=(
        "Creates output_dir / 'tracked'",
        "Writes tracked_atom_charges.csv and tracked_atom_charges.png",
        "Does NOT mutate Bader frame directories",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing root_dir, no frame dirs after slicing, or "
                "missing POSCAR / ACF.dat / POTCAR inside a frame"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Atom index out of bounds",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Empty atom list; negative indices; malformed index "
                "map; inconsistent frame atom count"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
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


def _handle_charge_tracked(params: dict[str, Any]) -> TaskResult:
    """Handler — route TrackedChargeReport fields."""
    from ..main import run_tracked_charge_analysis_with_report

    report = run_tracked_charge_analysis_with_report(**params)
    outputs: dict[str, str] = {
        "tracked_charge_csv": str(report.tracked_charge_csv),
        "tracked_charge_png": str(report.tracked_charge_png),
    }
    summary: dict[str, Any] = {
        "n_frames": int(report.n_frames),
        "atom_indices_xyz": [int(i) for i in report.atom_indices_xyz],
        "n_atoms_tracked": int(report.n_atoms_tracked),
    }
    return TaskResult(
        success=True,
        task="charge_tracked",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="charge_tracked",
    category="charge",
    description=(
        "Track Bader net charges for specified XYZ atom indices "
        "(0-based, ORIGINAL XYZ order, NOT POSCAR order).  Analyses "
        "existing Bader frame directories; does NOT generate Bader "
        "dirs or submit jobs."
    ),
    handler=_handle_charge_tracked,
    target_fn="md_analysis.main:run_tracked_charge_analysis_with_report",
    cli_codes=("225",),
    contract=_CHARGE_TRACKED_CONTRACT,
))


# 5. charge_counterion (CLI 226) — contract-backed.


_CHARGE_COUNTERION_CONTRACT = TaskContract(
    inputs={
        "output_dir": FieldSpec(
            description=(
                "Parent output directory; actual outputs go into "
                "output_dir / 'counterion_tracking'"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "root_dir": FieldSpec(
            description=(
                "Parent directory containing Bader frame subdirectories"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
            required=False, default=".",
        ),
        "metal_symbols": FieldSpec(
            description=(
                "Override default metal symbols for layer / water "
                "detection"
            ),
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "string"},
            },
            type="list[str] | None",
            required=False, default=None,
        ),
        "normal": FieldSpec(
            description="Cell axis perpendicular to the surface",
            json_schema={"type": "string", "enum": ["a", "b", "c"]},
            type="str",
            choices=("a", "b", "c"),
            required=False, default="c",
        ),
        "layer_tol_A": FieldSpec(
            description="Layer clustering tolerance",
            json_schema={"type": "number", "exclusiveMinimum": 0},
            type="float", unit="Angstrom",
            required=False, default=DEFAULT_LAYER_TOL_A,
        ),
        **_CHARGE_COMMON_FRAME_FIELDS,
    },
    outputs_artifacts={
        "counterion_charge_csv": FieldSpec(
            description="counterion_charges.csv (per-frame)",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "counterion_summary_csv": FieldSpec(
            description="counterion_summary.csv (per-atom summary)",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "counterion_charge_png": FieldSpec(
            description=(
                "counterion_charges.png — written only when at least "
                "one counterion was detected in any frame; omitted "
                "from TaskResult.outputs otherwise"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
    },
    outputs_metrics={
        "n_frames": FieldSpec(
            description="Number of frames analyzed",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "n_unique_counterions": FieldSpec(
            description=(
                "Number of unique XYZ-indexed atoms detected as "
                "counterions across all frames"
            ),
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "root_dir exists and contains at least one matching frame subdir",
        "Each matching subdir contains POSCAR, ACF.dat, POTCAR",
        "normal ∈ {'a', 'b', 'c'}",
    ),
    side_effects=(
        "Creates output_dir / 'counterion_tracking'",
        "Writes counterion_charges.csv and counterion_summary.csv",
        "Writes counterion_charges.png if at least one counterion was "
        "detected in any frame (skipped when detection is empty)",
        "Does NOT mutate Bader frame directories",
        "Does NOT decide whether detected ions are chemically relevant",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing root_dir, no frame dirs after slicing, or "
                "missing POSCAR / ACF.dat / POTCAR inside a frame"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Out-of-bounds indexing in remap step",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid normal; malformed Bader content; malformed "
                "POSCAR index map; invalid layer detection inputs"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
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


def _handle_charge_counterion(params: dict[str, Any]) -> TaskResult:
    """Handler — route CounterionChargeReport fields.

    ``counterion_charge_png`` is an optional artifact (only written
    when detection is non-empty); we include it in ``outputs`` **only**
    when the underlying report carries a concrete path, to uphold the
    agent-layer invariant that ``TaskResult.outputs`` lists real files
    on disk.
    """
    from ..main import run_counterion_charge_analysis_with_report

    report = run_counterion_charge_analysis_with_report(**params)
    outputs: dict[str, str] = {
        "counterion_charge_csv": str(report.counterion_charge_csv),
        "counterion_summary_csv": str(report.counterion_summary_csv),
    }
    if report.counterion_charge_png is not None:
        outputs["counterion_charge_png"] = str(report.counterion_charge_png)
    summary: dict[str, Any] = {
        "n_frames": int(report.n_frames),
        "n_unique_counterions": int(report.n_unique_counterions),
    }
    return TaskResult(
        success=True,
        task="charge_counterion",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="charge_counterion",
    category="charge",
    description=(
        "Per-frame counterion detection (non-water, non-metal atoms "
        "with non-zero Bader net charge) with charge time evolution "
        "and summary.  Analyses existing Bader frame directories; "
        "does NOT generate Bader dirs, submit jobs, or decide whether "
        "detected ions are chemically relevant."
    ),
    handler=_handle_charge_counterion,
    target_fn="md_analysis.main:run_counterion_charge_analysis_with_report",
    cli_codes=("226",),
    contract=_CHARGE_COUNTERION_CONTRACT,
))

# 6. run_all (composite) — contract-backed (Batch 5).
#    Intentionally the narrowest public surface: only the legacy
#    signature fields.  Advanced potential controls belong in
#    ``potential_full`` directly; ``run_all`` stays a convenience
#    composite of water_three_panel + potential_full.


_RUN_ALL_CONTRACT = TaskContract(
    inputs={
        "xyz_path": FieldSpec(
            description="CP2K XYZ trajectory file",
            json_schema={"type": "string"},
            type="Path", path_kind="file",
        ),
        "md_inp_path": FieldSpec(
            description=(
                "CP2K md.inp used to derive water cell dimensions "
                "when cell_abc is not provided"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "cell_abc": FieldSpec(
            description="Orthogonal cell lengths [a, b, c]",
            json_schema={
                "type": ["array", "null"],
                "items": {"type": "number", "exclusiveMinimum": 0},
                "minItems": 3, "maxItems": 3,
            },
            type="list[float] | None",
            unit="Angstrom", shape="(3,)",
            required=False, default=None,
        ),
        "output_dir": FieldSpec(
            description=(
                "Output root; artifacts are written under "
                "output_dir/water/ and "
                "output_dir/electrochemical/potential/"
            ),
            json_schema={"type": "string"},
            type="Path", path_kind="dir",
        ),
        "cube_pattern": FieldSpec(
            description=(
                "Glob for continuous-mode potential cube files; "
                "resolved relative to process current working directory"
            ),
            json_schema={"type": "string"},
            type="str",
            required=False, default="md-POTENTIAL-v_hartree-1_*.cube",
        ),
        "md_out_path": FieldSpec(
            description=(
                "CP2K md.out for Fermi energy; when omitted, potential "
                "skips Fermi-dependent sub-analyses"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file",
            required=False, default=None,
        ),
        "frame_start": FieldSpec(
            description="Optional slice start (0-based, inclusive)",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_end": FieldSpec(
            description="Optional slice exclusive upper bound",
            json_schema={"type": ["integer", "null"], "minimum": 0},
            type="int | None",
            required=False, default=None,
        ),
        "frame_step": FieldSpec(
            description="Optional slice step",
            json_schema={"type": ["integer", "null"], "minimum": 1},
            type="int | None",
            required=False, default=None,
        ),
        "verbose": FieldSpec(
            description="If True, show tqdm progress bars",
            json_schema={"type": "boolean"},
            type="bool",
            required=False, default=False,
        ),
    },
    outputs_artifacts={
        # Water artifacts (all always produced on success).
        "density_csv": FieldSpec(
            description="Water mass density CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "orientation_csv": FieldSpec(
            description="Water orientation-weighted density CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_profile_csv": FieldSpec(
            description="Adsorbed water layer profile CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_range_txt": FieldSpec(
            description="Adsorbed layer range / peak metrics TXT",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "adsorbed_theta_csv": FieldSpec(
            description="Adsorbed water θ angular distribution CSV",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "plot_png": FieldSpec(
            description="Water three-panel summary PNG",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        # Potential artifacts (conditional; only present when produced).
        "electrode_csv": FieldSpec(
            description=(
                "Electrode potential CSV; only when compute_u=True "
                "and Fermi energy is available"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "center_csv": FieldSpec(
            description=(
                "Slab-centred Hartree potential CSV; only when "
                "electrode is not computed"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "fermi_csv": FieldSpec(
            description=(
                "Fermi energy time-series CSV; only when "
                "compute_u=False and md_out_path is provided"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
        "phi_z_png": FieldSpec(
            description="φ(z) planar-average overlay PNG",
            json_schema={"type": "string"},
            type="Path", path_kind="file", category="artifact",
        ),
        "thickness_sensitivity_csv": FieldSpec(
            description=(
                "Thickness sensitivity sweep CSV; only when Fermi "
                "energy is available"
            ),
            json_schema={"type": ["string", "null"]},
            type="Path | None", path_kind="file", category="artifact",
            required=False,
        ),
    },
    outputs_metrics={
        "output_dir": FieldSpec(
            description="Echo of resolved output_dir root",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "water_output_dir": FieldSpec(
            description="Resolved water sub-directory",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "potential_output_dir": FieldSpec(
            description="Resolved potential sub-directory",
            json_schema={"type": "string"},
            type="str", category="metric",
        ),
        "n_artifacts": FieldSpec(
            description="Total artifact count across both leaves",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "water_n_artifacts": FieldSpec(
            description="Number of water artifacts (expected: 6)",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "potential_n_artifacts": FieldSpec(
            description="Number of potential artifacts produced",
            json_schema={"type": "integer"},
            type="int", category="metric",
        ),
        "ran_water": FieldSpec(
            description="Whether water analysis produced artifacts",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "ran_potential": FieldSpec(
            description=(
                "Whether potential analysis produced at least one "
                "artifact"
            ),
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_has_fermi": FieldSpec(
            description=(
                "Whether any Fermi-dependent potential artifact was "
                "produced"
            ),
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_ran_electrode": FieldSpec(
            description="Whether potential/electrode_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_ran_center": FieldSpec(
            description="Whether potential/center_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_ran_fermi": FieldSpec(
            description="Whether potential/fermi_csv was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_ran_phi_z": FieldSpec(
            description="Whether potential/phi_z_png was produced",
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "potential_ran_thickness_sensitivity": FieldSpec(
            description=(
                "Whether potential/thickness_sensitivity_csv was produced"
            ),
            json_schema={"type": "boolean"},
            type="bool", category="metric",
        ),
        "frame_start": FieldSpec(
            description="Echo of applied slice start",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
        "frame_end": FieldSpec(
            description="Echo of applied slice end",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
        "frame_step": FieldSpec(
            description="Echo of applied slice step",
            json_schema={"type": ["integer", "null"]},
            type="int | None", category="metric",
        ),
    },
    outputs_raw_model={},
    preconditions=(
        "xyz_path exists and is readable",
        "Either cell_abc is provided or md_inp_path yields cell dimensions",
        "If cell_abc is provided: length 3 with positive values",
        "cube_pattern matches files relative to the current working directory",
        "If Fermi-dependent potential outputs are desired: md_out_path exists",
        "Selected frame slice is valid and non-empty for both leaves",
    ),
    side_effects=(
        "Creates output_dir, output_dir/water/, and "
        "output_dir/electrochemical/potential/ (+ sub-analysis subdirs)",
        "Writes water CSV / TXT / PNG artifacts",
        "Writes potential CSV / PNG artifacts for sub-analyses that run",
        "May overwrite files with the same names",
        "Does NOT mutate input trajectory / md.inp / cube / md.out files",
        "Does NOT run charge / Bader / calibration / TI / "
        "constant-potential workflows",
        "Does NOT submit jobs",
    ),
    exceptions=(
        ExceptionMapping(
            exception_fqn="builtins.FileNotFoundError",
            triggered_by=(
                "Missing xyz_path; missing md_inp_path when needed for "
                "cell derivation; missing cube files; missing md_out_path "
                "when needed for Fermi-dependent outputs"
            ),
            error_type="file_not_found",
        ),
        ExceptionMapping(
            exception_fqn="builtins.IndexError",
            triggered_by="Invalid frame slice indexing",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.ValueError",
            triggered_by=(
                "Invalid cell_abc; invalid frame slice; malformed "
                "trajectory / cube / output input"
            ),
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.TypeError",
            triggered_by="Wrong input shape/type",
            error_type="validation",
        ),
        ExceptionMapping(
            exception_fqn="builtins.RuntimeError",
            triggered_by=(
                "A leaf wrapper succeeded but expected artifacts are "
                "missing; or water + potential produced an artifact-key "
                "collision (the composite refuses to silently overwrite)"
            ),
            error_type="analysis",
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


def _handle_run_all(params: dict[str, Any]) -> TaskResult:
    """Handler — route RunAllReport fields.

    Only artifacts that actually exist on disk make it into ``outputs``
    (guaranteed by the Batch 4 leaf wrappers).  Summary booleans mirror
    actual produced artifacts, not requested flags.
    """
    from ..main import run_all_with_report

    report = run_all_with_report(**params)
    outputs: dict[str, str] = {
        k: str(v) for k, v in report.artifacts.items()
    }
    summary: dict[str, Any] = {
        "output_dir": str(report.output_dir),
        "water_output_dir": str(report.water_output_dir),
        "potential_output_dir": str(report.potential_output_dir),
        "n_artifacts": int(report.n_artifacts),
        "water_n_artifacts": int(report.water_n_artifacts),
        "potential_n_artifacts": int(report.potential_n_artifacts),
        "ran_water": bool(report.ran_water),
        "ran_potential": bool(report.ran_potential),
        "potential_has_fermi": bool(report.potential_has_fermi),
        "potential_ran_electrode": bool(report.potential_ran_electrode),
        "potential_ran_center": bool(report.potential_ran_center),
        "potential_ran_fermi": bool(report.potential_ran_fermi),
        "potential_ran_phi_z": bool(report.potential_ran_phi_z),
        "potential_ran_thickness_sensitivity": bool(
            report.potential_ran_thickness_sensitivity
        ),
        "frame_start": (
            int(report.frame_start) if report.frame_start is not None else None
        ),
        "frame_end": (
            int(report.frame_end) if report.frame_end is not None else None
        ),
        "frame_step": (
            int(report.frame_step) if report.frame_step is not None else None
        ),
    }
    return TaskResult(
        success=True,
        task="run_all",
        outputs=outputs,
        summary=summary,
    )


register(TaskDef(
    name="run_all",
    category="composite",
    description=(
        "Convenience composite: water_three_panel + potential_full "
        "under the canonical output layout "
        "(output_dir/water/ and output_dir/electrochemical/potential/)."
        "  This is NOT a prompt-level workflow coordinator and does "
        "NOT run charge / Bader / calibration / TI / constant-potential "
        "workflows.  Advanced potential controls belong to potential_full "
        "directly."
    ),
    handler=_handle_run_all,
    target_fn="md_analysis.main:run_all_with_report",
    contract=_RUN_ALL_CONTRACT,
))

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
    """Dedicated handler — route artifacts + metrics from the report."""
    from ..electrochemical.calibration.CalibrationWorkflow import (
        calibrate_with_report,
    )

    report = calibrate_with_report(**params)

    outputs: dict[str, str] = {"calibration_json": str(report.calibration_json)}
    if report.calibration_csv is not None:
        outputs["calibration_csv"] = str(report.calibration_csv)
    if report.calibration_png is not None:
        outputs["calibration_png"] = str(report.calibration_png)

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
    target_fn=(
        "md_analysis.electrochemical.calibration.CalibrationWorkflow"
        ":calibrate_with_report"
    ),
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
    """Dedicated handler — route JSON-serializable metrics from the report."""
    from ..electrochemical.calibration.CalibrationWorkflow import (
        predict_potential_with_report,
    )

    report = predict_potential_with_report(**params)

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
    target_fn=(
        "md_analysis.electrochemical.calibration.CalibrationWorkflow"
        ":predict_potential_with_report"
    ),
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
                "md_analysis.utils.RestartParser.ColvarParser.ColvarParseError"
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
