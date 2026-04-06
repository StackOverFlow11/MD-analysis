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
