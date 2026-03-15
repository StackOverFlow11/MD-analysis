"""Private helpers for discovering and sorting Bader frame directories."""

from __future__ import annotations

import re
from pathlib import Path

_T_VALUE_RE = re.compile(r"_t(\d+)")
_I_VALUE_RE = re.compile(r"_i(\d+)")


def _extract_t_value(dirname: str) -> int:
    """Extract numeric t value from a directory name like ``bader_t50_i0``."""
    m = _T_VALUE_RE.search(dirname)
    return int(m.group(1)) if m else 0


def _extract_step_and_time(dirname: str) -> tuple[int, int]:
    """Extract ``(step, time_fs)`` from directory name ``bader_t{time}_i{step}``.

    Returns ``(0, 0)`` for fields not found.
    """
    m_i = _I_VALUE_RE.search(dirname)
    step = int(m_i.group(1)) if m_i else 0
    time_fs = _extract_t_value(dirname)
    return step, time_fs


def _sorted_frame_dirs(root: Path, dir_pattern: str) -> list[Path]:
    """Discover and numerically sort frame subdirectories by t value."""
    frame_dirs = sorted(root.glob(dir_pattern), key=lambda p: _extract_t_value(p.name))
    if not frame_dirs:
        raise FileNotFoundError(
            f"No subdirectories matching '{dir_pattern}' in {root}"
        )
    return frame_dirs
