"""Private helpers for discovering and sorting Bader frame directories.

Delegates to the shared :mod:`md_analysis.utils.io._frame_discovery` helpers.
Kept as a thin compatibility layer so that call-sites in ``SurfaceCharge``,
``AtomCharges`` and ``BaderData`` continue to import from this module.
"""

from __future__ import annotations

from pathlib import Path

from ....utils.io._frame_discovery import (
    discover_frame_dirs as _discover_frame_dirs,
    extract_step_time_from_dirname,
)


def _extract_step_and_time(dirname: str) -> tuple[int, int]:
    """Extract ``(step, time_fs)`` from directory name ``bader_t{time}_i{step}``.

    Returns ``(0, 0)`` if the name does not match the expected pattern.
    """
    parsed = extract_step_time_from_dirname(dirname)
    return parsed if parsed is not None else (0, 0)


def _extract_t_value(dirname: str) -> int:
    """Extract the ``_t`` integer (time in fs) from a frame directory name.

    Preserved for ``SurfaceCharge`` CSV output — the ``"step"`` column there
    historically stores the ``_t`` value (time in fs), not the ``_i`` step.
    Returns ``0`` when the name does not match.
    """
    parsed = extract_step_time_from_dirname(dirname)
    return parsed[1] if parsed is not None else 0


def _sorted_frame_dirs(root: Path, dir_pattern: str) -> list[Path]:
    """Discover and numerically sort frame subdirectories by step (``_i`` value)."""
    return _discover_frame_dirs(root, dir_pattern, required=True)
