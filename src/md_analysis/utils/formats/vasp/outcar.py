"""VASP OUTCAR file parser — placeholder.

The VASP ``OUTCAR`` file is the canonical run log, containing per-step
Fermi level, total energy, force statistics, and SCF history. It is the
VASP counterpart of the CP2K ``md.out`` / ``sp.out`` files parsed by
``utils.formats.cp2k.stdout``.

Phase 8 status
--------------
Deliberate placeholder. Every parsing entry point raises
:class:`NotImplementedError`. ``engines.vasp.VASPParser`` is NOT in the
default parser registry, so ``infer_parser`` cannot silently dispatch
to a stub.

See ``vasp_report.py`` for the broader Phase 8 policy on VASP
placeholders.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Type-only forward reference; the canonical home of ``FermiRecord``
    # is ``md_analysis.engines.models``. Importing at runtime would
    # invert the architectural direction (``engines`` ->
    # ``utils/formats``); the TYPE_CHECKING guard keeps the annotations
    # useful for static analysis without creating a runtime
    # ``utils/formats`` -> ``engines`` edge.
    from ....engines.models import FermiRecord

_NOT_IMPLEMENTED_MSG = (
    "VASP OUTCAR parsing is not implemented yet. See "
    "context4agent/requirements/overall_reconstruction_plan.md for the "
    "planned VASP rollout."
)


def parse_outcar_fermi(outcar_path: str | Path) -> "list[FermiRecord]":  # noqa: ARG001
    """Parse the per-step Fermi level (Hartree) series from a VASP
    ``OUTCAR`` file.

    Returns the same engine-neutral ``list[FermiRecord]`` shape that
    :func:`md_analysis.engines.cp2k.read_fermi_series` does.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


__all__ = [
    "parse_outcar_fermi",
]
