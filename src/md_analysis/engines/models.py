"""Engine-neutral data models.

Frozen dataclasses returned by engine adapter modules
(``engines.cp2k``, future ``engines.vasp``). Upper analysis layers
(``electrochemical``, ``enhanced_sampling``, ``water``) depend on these
neutral types rather than engine-specific parser output.

Status (Phase 5b)
-----------------
- ``PotentialFrame``      — moved here from ``electrochemical.potential
                             ._frame_source``; field set unchanged
- ``ConstraintMetadata``  — canonical name, renamed from
                             ``ColvarRestart`` in Phase 5b; field set
                             unchanged
- ``LambdaSeries``        — canonical name, renamed from
                             ``LagrangeMultLog`` in Phase 5b; field set
                             unchanged

The legacy ``ColvarRestart`` / ``LagrangeMultLog`` names are kept as
module-level aliases in ``utils.formats.cp2k_colvar`` so the
``enhanced_sampling/_parsers.py`` shim and existing tests keep working
during the transition. They will be removed in a later cleanup phase.

Additional neutral types (``ConstraintPoint``, ``ConstraintSet``,
``ConstraintRun``, ``FermiRecord``) listed in the refactor plan are
deferred until concrete consumers exist (Phase 7 onward).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..utils.formats.cube import CubeHeader
from ..utils.formats.cp2k_colvar import ConstraintMetadata, LambdaSeries

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


__all__ = [
    "PotentialFrame",
    "ConstraintMetadata",
    "LambdaSeries",
]
