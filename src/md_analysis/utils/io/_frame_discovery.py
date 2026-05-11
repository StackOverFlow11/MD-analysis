"""Private helpers for discovering and sorting frame work directories.

Used by:
- ``electrochemical.charge.Bader.*`` (Bader frame dirs ``bader_t{time}_i{step}``)
- ``electrochemical.potential._frame_source`` (SP frame dirs ``potential_t{time}_i{step}``)

Directory-name convention
-------------------------
Frames are organised in per-frame subdirectories whose names embed the CP2K
MD step (``i``) and time in fs (``t``) written by ``BaderGen``/``PotentialGen``:

    bader_t{time}_i{step}       (charge / Bader)
    potential_t{time}_i{step}   (potential / SP)

Both the ``_t(\\d+)`` and ``_i(\\d+)`` segments must be present. Sort order is
by ``(time_fs, step)`` — both are strictly monotonic within a single MD
trajectory, so the tuple ordering agrees with either single key in practice
while remaining well-defined when one of the two is held constant.
"""

from __future__ import annotations

import re
from pathlib import Path

# Matches ``_t{int}_i{int}`` anywhere in a directory name.
FRAME_DIR_STEP_TIME_RE = re.compile(r"_t(\d+)_i(\d+)")


def extract_step_time_from_dirname(dirname: str) -> tuple[int, int] | None:
    """Extract ``(step, time_fs)`` from a frame directory name.

    Returns
    -------
    ``(step, time_fs)`` tuple, or ``None`` if the name does not contain a
    ``_t<int>_i<int>`` segment.
    """
    m = FRAME_DIR_STEP_TIME_RE.search(dirname)
    if m is None:
        return None
    time_fs = int(m.group(1))
    step = int(m.group(2))
    return step, time_fs


def discover_frame_dirs(
    root: Path,
    glob_pattern: str,
    *,
    required: bool = True,
) -> list[Path]:
    """Discover frame sub-directories and sort them numerically by step.

    Directories whose names do not match ``_t<int>_i<int>`` are silently
    dropped. The remaining directories are sorted by the step value (``i``).

    Parameters
    ----------
    root
        Parent directory to search under.
    glob_pattern
        Glob pattern applied to ``root`` (e.g. ``"bader_t*_i*"``).
    required
        If True (default), raise ``FileNotFoundError`` when no matching
        sub-directories are discovered. If False, an empty list is returned.
    """
    dirs: list[tuple[Path, int, int]] = []
    for d in root.glob(glob_pattern):
        if not d.is_dir():
            continue
        parsed = extract_step_time_from_dirname(d.name)
        if parsed is None:
            continue
        step, time_fs = parsed
        dirs.append((d, time_fs, step))

    if not dirs:
        if required:
            raise FileNotFoundError(
                f"No subdirectories matching '{glob_pattern}' with "
                f"_t<int>_i<int> in {root}"
            )
        return []

    # Sort by (time_fs, step) — both monotonic for real trajectories.
    dirs.sort(key=lambda x: (x[1], x[2]))
    return [p for p, _, _ in dirs]
