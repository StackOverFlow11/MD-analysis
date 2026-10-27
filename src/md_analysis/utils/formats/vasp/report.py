"""VASP REPORT file parser — placeholder.

The VASP ``REPORT`` file emitted by SHAKE/RATTLE constraint MD carries
the per-step Lagrange-multiplier series plus the constraint metadata
(target value, growth rate, timestep). It is the VASP counterpart of
the CP2K ``.restart`` + ``.LagrangeMultLog`` pair.

Phase 8 status
--------------
This module is a deliberate placeholder. Every parsing entry point
raises :class:`NotImplementedError` so that a future contributor can
fill in the conversion without first reverse-engineering call sites.
Until then:

- ``engines.vasp.VASPParser`` is NOT registered in the default parser
  registry (see ``engines/__init__.py``), so ``infer_parser`` cannot
  silently dispatch to a stub.
- Direct construction of ``VASPParser`` and these helpers is allowed,
  but every call raises :class:`NotImplementedError` with a clear
  pointer to the refactor plan.

Do NOT add real parsing here without rolling out the corresponding
``engines/vasp.py`` and registry update; the goal of the placeholder is
to make the extension point obvious, not to provide a half-working
default.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Type-only forward reference; the canonical home of these neutral
    # dataclasses is ``md_analysis.engines.models``. Importing them at
    # runtime would invert the architectural direction
    # (``engines`` -> ``utils/formats``); the TYPE_CHECKING guard keeps
    # the annotations meaningful for static analysis without creating a
    # runtime ``utils/formats`` -> ``engines`` edge.
    from ....engines.models import ConstraintMetadata, LambdaSeries

_NOT_IMPLEMENTED_MSG = (
    "VASP REPORT parsing is not implemented yet. See "
    "context4agent/requirements/overall_reconstruction_plan.md for the "
    "planned VASP rollout."
)


def parse_vasp_report_metadata(report_path: str | Path) -> "ConstraintMetadata":  # noqa: ARG001
    """Parse constraint metadata (target / growth / timestep / ...) from
    a VASP ``REPORT`` file.

    Returns the same engine-neutral :class:`ConstraintMetadata` shape
    that the CP2K parser does.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


def parse_vasp_report_lambda_series(report_path: str | Path) -> "LambdaSeries":  # noqa: ARG001
    """Parse the per-step Lagrange-multiplier (constraint force) series
    from a VASP ``REPORT`` file.

    Returns the same engine-neutral :class:`LambdaSeries` shape that the
    CP2K LagrangeMultLog parser does.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


__all__ = [
    "parse_vasp_report_metadata",
    "parse_vasp_report_lambda_series",
]
