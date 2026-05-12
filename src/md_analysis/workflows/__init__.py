"""Programmatic workflow facades for md_analysis.

This package partitions the public ``run_*`` programmatic entry points
by domain (water, potential, charge, calibration, enhanced sampling,
scripts, composite). Each ``run_*`` function returns a
:class:`~md_analysis.workflows.models.WorkflowResult` describing the
artifacts produced on disk and lightweight metadata about the run.

``md_analysis.main`` re-exports a stable subset of these entry points
for backwards compatibility; new code should prefer the explicit
``md_analysis.workflows.<domain>`` import paths.
"""

from __future__ import annotations

from .calibration import run_calibration_fit, run_calibration_predict
from .charge import run_counterion_charge, run_surface_charge, run_tracked_charge
from .composite import run_interface_analysis
from .enhanced_sampling import (
    run_slowgrowth_publication_plot,
    run_slowgrowth_quick_plot,
    run_ti_constant_potential_correction,
    run_ti_full_analysis,
    run_ti_single_diagnostics,
)
from .models import MissingArtifactError, WorkflowResult, require_artifacts_exist
from .potential import run_potential_full
from .scripts import (
    run_bader_batch,
    run_bader_single,
    run_potential_batch,
    run_potential_single,
    run_sp_batch,
    run_sp_single,
    run_ti_batch,
    run_ti_single,
)
from .water import run_water_three_panel

__all__ = [
    # Models
    "MissingArtifactError",
    "WorkflowResult",
    "require_artifacts_exist",
    # Leaf workflows (Phase 2)
    "run_water_three_panel",
    "run_potential_full",
    "run_surface_charge",
    "run_tracked_charge",
    "run_counterion_charge",
    # Calibration workflows (Phase 4.1)
    "run_calibration_fit",
    "run_calibration_predict",
    # Enhanced-sampling workflows (Phase 4.2)
    "run_slowgrowth_quick_plot",
    "run_slowgrowth_publication_plot",
    "run_ti_single_diagnostics",
    "run_ti_full_analysis",
    "run_ti_constant_potential_correction",
    # Scripts / work-directory generators (Phase 4.3)
    "run_bader_single",
    "run_bader_batch",
    "run_ti_single",
    "run_ti_batch",
    "run_potential_single",
    "run_potential_batch",
    "run_sp_single",
    "run_sp_batch",
    # Composite workflows (Phase 5)
    "run_interface_analysis",
]
