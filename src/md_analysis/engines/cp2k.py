"""CP2K engine adapter.

Implements :class:`ConstraintMDParser` for CP2K constraint-MD point
directories (``.restart`` + ``.LagrangeMultLog`` file pair) and provides
a small module-level facade for the most frequently used CP2K output
readers.

Phase 7b1 status
----------------
- ``CP2KParser``            — Protocol implementation (Phase 5a).
- ``read_constraint_metadata`` / ``read_lambda_series`` —
  Path-friendly facade over ``CP2KParser`` that accepts ``str | Path``
  and constructs the parser internally. Use this when you only need the
  result; use ``CP2KParser`` directly when you need the Protocol
  surface or want to reuse the same parser across many directories.
- ``read_fermi_series``     — Reads ``md.out`` Fermi entries via
  ``utils.formats.cp2k_stdout.parse_md_out_fermi`` and converts each
  legacy dict into a typed ``FermiRecord``. The underlying parser
  function still returns ``list[dict]`` for the existing
  dict-based callers (see ``electrochemical.potential.CenterPotential``).

``read_cube_frames`` is intentionally NOT added yet: that helper would
require moving the multi-frame discovery logic currently in
``electrochemical.potential._frame_source`` down into ``engines``, which
is a layer-direction change deferred to Phase 7b2.
"""

from __future__ import annotations

from pathlib import Path

from ..utils.formats.cp2k_colvar import (
    parse_colvar_restart,
    parse_lagrange_mult_log,
)
from ..utils.formats.cp2k_stdout import parse_md_out_fermi
from .models import ConstraintMetadata, FermiRecord, LambdaSeries


class CP2KParser:
    """Parser for CP2K constraint-MD point directories.

    Recognises directories containing ``*.restart`` and
    ``*.LagrangeMultLog`` files (the standard CP2K output pair).
    """

    name = "cp2k"

    def is_constraint_directory(self, directory: Path) -> bool:
        if not directory.is_dir():
            return False
        try:
            self._find_restart(directory)
            self._find_log(directory)
        except FileNotFoundError:
            return False
        return True

    def parse_metadata(self, directory: Path) -> ConstraintMetadata:
        return parse_colvar_restart(self._find_restart(directory))

    def parse_lambda_series(self, directory: Path) -> LambdaSeries:
        return parse_lagrange_mult_log(self._find_log(directory))

    # ------------------------------------------------------------------
    # File discovery (private)
    # ------------------------------------------------------------------

    @staticmethod
    def _find_restart(directory: Path) -> Path:
        """Find the primary .restart file (skips .bak and .RESTART.wfn)."""
        candidates = sorted(directory.glob("*.restart"))
        candidates = [
            p for p in candidates
            if ".bak" not in p.name and "RESTART.wfn" not in p.name
        ]
        if not candidates:
            raise FileNotFoundError(f"No .restart file in {directory}")
        # Prefer the one with highest suffix number (e.g. cMD-1_1500.restart)
        return candidates[-1]

    @staticmethod
    def _find_log(directory: Path) -> Path:
        """Find the .LagrangeMultLog file."""
        candidates = list(directory.glob("*.LagrangeMultLog"))
        if not candidates:
            raise FileNotFoundError(f"No .LagrangeMultLog file in {directory}")
        return candidates[0]


# ---------------------------------------------------------------------------
# Module-level facade (Phase 7b1)
# ---------------------------------------------------------------------------


def read_constraint_metadata(
    directory: str | Path,
) -> ConstraintMetadata:
    """Read constraint-MD metadata from a CP2K point directory.

    Thin path-friendly wrapper around :meth:`CP2KParser.parse_metadata`;
    accepts ``str`` or ``Path`` and uses a fresh ``CP2KParser`` instance.
    """
    return CP2KParser().parse_metadata(Path(directory))


def read_lambda_series(
    directory: str | Path,
) -> LambdaSeries:
    """Read the Lagrange-multiplier (λ(t)) series from a CP2K point directory.

    Thin path-friendly wrapper around
    :meth:`CP2KParser.parse_lambda_series`; accepts ``str`` or ``Path``
    and uses a fresh ``CP2KParser`` instance.
    """
    return CP2KParser().parse_lambda_series(Path(directory))


def read_fermi_series(
    md_out_path: str | Path,
) -> list[FermiRecord]:
    """Read the Fermi-energy series from a CP2K ``md.out`` file.

    Returns a list of engine-neutral :class:`FermiRecord` rows
    (``step``, ``time_fs``, ``fermi_raw`` in Hartree). The underlying
    parser (:func:`utils.formats.cp2k_stdout.parse_md_out_fermi`) still
    returns the legacy ``list[dict]`` shape and is unchanged.
    """
    legacy = parse_md_out_fermi(Path(md_out_path))
    return [FermiRecord.from_legacy_dict(d) for d in legacy]


__all__ = [
    "CP2KParser",
    "read_constraint_metadata",
    "read_lambda_series",
    "read_fermi_series",
]
