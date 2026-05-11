"""Engine-neutral data models.

Frozen dataclasses returned by engine adapter modules
(``engines.cp2k``, future ``engines.vasp``). Upper analysis layers
(``electrochemical``, ``enhanced_sampling``, ``water``) should depend on
these neutral types rather than engine-specific parser output.

Phase 5a status
---------------
This module currently re-uses two existing dataclasses without rename:

- ``PotentialFrame``     — moved here from ``electrochemical.potential
                            ._frame_source``; field set unchanged
- ``ColvarRestart``      — re-exported from ``utils.formats.cp2k_colvar``
                            for transitional convenience; will be
                            renamed to ``ConstraintMetadata`` in
                            Phase 5b
- ``LagrangeMultLog``    — re-exported from ``utils.formats.cp2k_colvar``;
                            will be renamed to ``LambdaSeries`` in
                            Phase 5b

Additional neutral types (``ConstraintPoint``, ``ConstraintSet``,
``ConstraintRun``, ``FermiRecord``) listed in the refactor plan are
deferred until concrete consumers exist (Phase 5b / Phase 7).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..utils.formats.cube import CubeHeader
# Phase 5b will rename these locally:
#   ColvarRestart    -> ConstraintMetadata
#   LagrangeMultLog  -> LambdaSeries
from ..utils.formats.cp2k_colvar import ColvarRestart, LagrangeMultLog

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
    "ColvarRestart",
    "LagrangeMultLog",
]
