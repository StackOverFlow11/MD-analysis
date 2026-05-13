"""Engine-neutral data models.

Frozen dataclasses returned by engine adapter modules
(``engines.cp2k``, future ``engines.vasp``). Upper analysis layers
(``electrochemical``, ``enhanced_sampling``, ``water``) depend on these
neutral types rather than engine-specific parser output.

Status (Phase 7b1)
------------------
- ``PotentialFrame``      — moved here from ``electrochemical.potential
                             ._frame_source``; field set unchanged
- ``ConstraintMetadata``  — canonical name, renamed from
                             ``ColvarRestart`` in Phase 5b; field set
                             unchanged
- ``LambdaSeries``        — canonical name, renamed from
                             ``LagrangeMultLog`` in Phase 5b; field set
                             unchanged
- ``FermiRecord``         — new typed record for one ``(step, time_fs,
                             fermi_raw)`` triple, mirroring the legacy
                             dict shape so existing callers (e.g.
                             ``CenterPotential``) keep working unchanged

The legacy ``ColvarRestart`` / ``LagrangeMultLog`` names are kept as
module-level aliases in ``utils.formats.cp2k.colvar`` so existing tests
keep working during the transition. They will be removed in a later
cleanup phase.

Additional neutral types (``ConstraintPoint``, ``ConstraintSet``,
``ConstraintRun``) listed in the refactor plan are deferred until
concrete consumers exist.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..utils.formats.common.cube import CubeHeader
from ..utils.formats.cp2k.colvar import ConstraintMetadata, LambdaSeries

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]


@dataclass(frozen=True)
class PotentialFrame:
    """One frame of potential analysis data (immutable, engine-neutral).

    Produced by ``engines.cp2k`` when discovering cube/SP frames; consumed
    by ``electrochemical.potential`` analysis routines.

    Field set is intentionally unchanged from the Phase 4 definition
    that previously lived in ``electrochemical.potential._frame_source``.
    """

    step: int
    time_fs: float | None
    cube_path: Path
    header: CubeHeader
    values: np.ndarray
    fermi_raw: float | None  # Hartree; None if unavailable
    atoms: Atoms | None  # for interface detection (with cell set)


@dataclass(frozen=True)
class FermiRecord:
    """One ``(step, time_fs, fermi_raw)`` row from a CP2K stdout log.

    Engine-neutral typed counterpart to the legacy dict-of-records shape
    returned by :func:`md_analysis.utils.formats.cp2k.stdout.parse_md_out_fermi`.
    The fields mirror that dict exactly so legacy callers can keep
    using their dict-style access without change while new callers
    consume the typed model via the
    :func:`md_analysis.engines.cp2k.read_fermi_series` facade.

    ``fermi_raw`` is in Hartree.
    """

    step: int
    time_fs: float | None
    fermi_raw: float

    @classmethod
    def from_legacy_dict(cls, d: dict) -> "FermiRecord":
        """Bridge a legacy ``parse_md_out_fermi`` dict into a typed record.

        Strictly a one-way converter from the legacy dict shape; the
        parser function itself still returns ``list[dict]`` and is NOT
        modified.
        """
        return cls(
            step=int(d["step"]),
            time_fs=d["time_fs"],
            fermi_raw=float(d["fermi_raw"]),
        )


__all__ = [
    "PotentialFrame",
    "ConstraintMetadata",
    "LambdaSeries",
    "FermiRecord",
]
