"""File discovery and batch parsing for constrained-TI data.

Pure I/O layer — does NOT compute integration weights or trim equilibration.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

from .models import TIPointDefinition

logger = logging.getLogger(__name__)

# Patterns for directory-name -> xi extraction
_TI_TARGET_RE = re.compile(r"^ti_target_([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)$")
_XI_RE = re.compile(r"^xi_([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)$")


def _find_restart(directory: Path) -> Path:
    """Find the primary .restart file in a constraint-point directory."""
    # Prefer the base restart (no bak suffix), latest numbered if multiple
    candidates = sorted(directory.glob("*.restart"))
    # Filter out .bak files and checkpoint restarts
    candidates = [
        p
        for p in candidates
        if ".bak" not in p.name and "RESTART.wfn" not in p.name
    ]
    if not candidates:
        raise FileNotFoundError(f"No .restart file found in {directory}")
    # Prefer the one with highest suffix number (e.g., cMD-1_1500.restart)
    return candidates[-1]


def _find_log(directory: Path) -> Path:
    """Find the .LagrangeMultLog file in a constraint-point directory."""
    candidates = list(directory.glob("*.LagrangeMultLog"))
    if not candidates:
        raise FileNotFoundError(f"No .LagrangeMultLog file found in {directory}")
    return candidates[0]


def discover_ti_points(
    root_dir: Path,
    *,
    pattern: str = "auto",
) -> list[TIPointDefinition]:
    """Discover constraint-point directories and their files.

    Parameters
    ----------
    root_dir : Path
        Root directory containing constraint-point subdirectories.
    pattern : str
        Discovery pattern: "ti_target", "xi", or "auto".

    Returns
    -------
    list[TIPointDefinition]
        Sorted by xi value.

    Raises
    ------
    FileNotFoundError
        If no matching directories are found.
    """
    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory does not exist: {root}")

    points: list[TIPointDefinition] = []

    if pattern in ("ti_target", "auto"):
        for d in sorted(root.iterdir()):
            if not d.is_dir():
                continue
            m = _TI_TARGET_RE.match(d.name)
            if m:
                xi = float(m.group(1))
                try:
                    restart = _find_restart(d)
                    log = _find_log(d)
                    points.append(TIPointDefinition(xi=xi, restart_path=restart, log_path=log))
                except FileNotFoundError as e:
                    logger.warning("Skipping %s: %s", d.name, e)

    if not points and pattern in ("xi", "auto"):
        for d in sorted(root.iterdir()):
            if not d.is_dir():
                continue
            m = _XI_RE.match(d.name)
            if m:
                xi = float(m.group(1))
                try:
                    restart = _find_restart(d)
                    log = _find_log(d)
                    points.append(TIPointDefinition(xi=xi, restart_path=restart, log_path=log))
                except FileNotFoundError as e:
                    logger.warning("Skipping %s: %s", d.name, e)

    if not points:
        raise FileNotFoundError(
            f"No constraint-point directories found in {root} "
            f"(tried pattern={pattern!r})."
        )

    # Sort by xi
    points.sort(key=lambda p: p.xi)
    return points


def load_ti_series(
    point_defs: list[TIPointDefinition],
) -> list[tuple[float, np.ndarray, float]]:
    """Parse restart + LagrangeMultLog for each point.

    Returns full (untrimmed) series. Equilibration trimming is
    the workflow's responsibility.

    Parameters
    ----------
    point_defs : list[TIPointDefinition]

    Returns
    -------
    list[tuple[float, np.ndarray, float]]
        (xi, lambda_series, dt_fs) for each point.
    """
    from ...utils.RestartParser.ColvarParser import ColvarMDInfo

    results = []
    for pdef in point_defs:
        md_info = ColvarMDInfo.from_paths(
            str(pdef.restart_path), str(pdef.log_path)
        )
        lambda_series = md_info.lagrange.collective_shake
        dt_fs = float(md_info.restart.timestep_fs)
        results.append((pdef.xi, lambda_series, dt_fs))

    return results
