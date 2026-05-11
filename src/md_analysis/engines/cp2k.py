"""CP2K engine adapter.

Implements :class:`ConstraintMDParser` for CP2K constraint-MD point
directories (``.restart`` + ``.LagrangeMultLog`` file pair).

Phase 5a status
---------------
Only the constraint-MD parser surface is hosted here today. Phase 7 will
extend this module with CP2K cube / md.out / xyz facade APIs (currently
those live in ``utils.formats.{cube,cp2k_stdout,cp2k_xyz}`` once Phase 7
finishes the split).
"""

from __future__ import annotations

from pathlib import Path

from ..utils.formats.cp2k_colvar import (
    parse_colvar_restart,
    parse_lagrange_mult_log,
)
from .models import ColvarRestart, LagrangeMultLog


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

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        return parse_colvar_restart(self._find_restart(directory))

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
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


__all__ = ["CP2KParser"]
