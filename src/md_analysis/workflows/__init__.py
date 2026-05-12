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

from .models import MissingArtifactError, WorkflowResult, require_artifacts_exist

__all__ = [
    "MissingArtifactError",
    "WorkflowResult",
    "require_artifacts_exist",
]
