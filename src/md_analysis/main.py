"""Programmatic entry points for ``md_analysis`` analysis workflows.

This module is a **thin re-export facade** over
:mod:`md_analysis.workflows`. Every public name lives in
``workflows.<domain>`` and is re-exported here purely for ergonomic
top-level access (``from md_analysis.main import run_water_three_panel``
matches the historical ``from md_analysis.main import ...`` shape).

The legacy ``run_*_analysis`` and ``run_all`` names that previously
lived here were removed in Phase 7a of the entrance refactor. New
code should depend on either ``md_analysis.workflows`` or this
facade, not on the removed legacy names.

Every public ``run_*`` function returns a
:class:`md_analysis.workflows.WorkflowResult` with file artifacts on
``.artifacts`` and lightweight scalars on ``.metadata``; complex
report objects (TI / calibration / etc.) live on ``.extra``.
"""

from __future__ import annotations

from .workflows import (
    WorkflowResult,
    run_bader_batch,
    run_bader_single,
    run_calibration_fit,
    run_calibration_predict,
    run_counterion_charge,
    run_interface_analysis,
    run_potential_batch,
    run_potential_full,
    run_potential_single,
    run_slowgrowth_publication_plot,
    run_slowgrowth_quick_plot,
    run_sp_batch,
    run_sp_single,
    run_surface_charge,
    run_ti_batch,
    run_ti_constant_potential_correction,
    run_ti_full_analysis,
    run_ti_single,
    run_ti_single_diagnostics,
    run_tracked_charge,
    run_water_three_panel,
)

__all__ = [
    "WorkflowResult",
    # Water (1)
    "run_water_three_panel",
    # Potential (1)
    "run_potential_full",
    # Charge (3)
    "run_surface_charge",
    "run_tracked_charge",
    "run_counterion_charge",
    # Calibration (2)
    "run_calibration_fit",
    "run_calibration_predict",
    # Enhanced sampling (5)
    "run_slowgrowth_quick_plot",
    "run_slowgrowth_publication_plot",
    "run_ti_single_diagnostics",
    "run_ti_full_analysis",
    "run_ti_constant_potential_correction",
    # Scripts / work-directory generators (8)
    "run_bader_single",
    "run_bader_batch",
    "run_ti_single",
    "run_ti_batch",
    "run_potential_single",
    "run_potential_batch",
    "run_sp_single",
    "run_sp_batch",
    # Composite (1)
    "run_interface_analysis",
]
