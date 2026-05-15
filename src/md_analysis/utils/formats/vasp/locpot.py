"""VASP LOCPOT file parser — placeholder.

The VASP ``LOCPOT`` file is the locally averaged electrostatic potential
written on a regular grid. It is the VASP counterpart of the CP2K
Hartree-potential ``.cube`` file parsed by ``utils.formats.common.cube``.

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

_NOT_IMPLEMENTED_MSG = (
    "VASP LOCPOT parsing is not implemented yet. See "
    "context4agent/requirements/overall_reconstruction_plan.md for the "
    "planned VASP rollout."
)


def read_locpot(locpot_path: str | Path) -> tuple:  # noqa: ARG001
    """Read a VASP ``LOCPOT`` file (planar-averaged electrostatic
    potential on a regular grid).

    Real implementation should return the same data quantities as the
    CP2K cube reader (header + 3-D values array), in whatever
    engine-neutral shape Phase 8 settles on; until then the call is a
    hard error.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


def read_locpot_plane_avg(locpot_path: str | Path, *, axis: str = "c") -> tuple:  # noqa: ARG001
    """Compute the plane-averaged ``φ(z)`` profile from a VASP ``LOCPOT``
    file along the requested cell axis.

    Engine-neutral counterpart of
    :func:`md_analysis.utils.formats.common.cube.plane_avg_phi_z_ev`.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


__all__ = [
    "read_locpot",
    "read_locpot_plane_avg",
]
