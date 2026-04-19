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
    reverse: bool = False,
    strict: bool = False,
) -> list[TIPointDefinition]:
    """Discover constraint-point directories and their files.

    Parameters
    ----------
    root_dir : Path
        Root directory containing constraint-point subdirectories.
    pattern : str
        Discovery pattern: "ti_target", "xi", or "auto".
    reverse : bool
        If True, sort by xi descending (initial state = max ξ).
    strict : bool
        If ``False`` (default, human-facing behaviour), matched directories
        missing a ``.restart`` or ``.LagrangeMultLog`` file are logged and
        skipped so the remaining points can still be analysed.

        If ``True`` (agent-facing behaviour), any matched directory missing
        required files raises :class:`FileNotFoundError` immediately — used
        by ``run_ti_full_from_root()`` so the agent layer can surface the
        problem as ``file_not_found`` instead of silently analysing a
        subset.

    Returns
    -------
    list[TIPointDefinition]
        Sorted by xi value (ascending by default, descending if reverse).

    Raises
    ------
    FileNotFoundError
        If no matching directories are found, or (when ``strict=True``) if
        a matched directory is missing a required file.
    """
    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory does not exist: {root}")

    points: list[TIPointDefinition] = []

    def _collect(regex: re.Pattern) -> None:
        for d in sorted(root.iterdir()):
            if not d.is_dir():
                continue
            m = regex.match(d.name)
            if not m:
                continue
            xi = float(m.group(1))
            try:
                restart = _find_restart(d)
                log = _find_log(d)
                points.append(
                    TIPointDefinition(xi=xi, restart_path=restart, log_path=log)
                )
            except FileNotFoundError as e:
                if strict:
                    raise FileNotFoundError(
                        f"Matched TI directory {d.name!r} is missing a "
                        f"required file: {e}"
                    ) from e
                logger.warning("Skipping %s: %s", d.name, e)

    if pattern in ("ti_target", "auto"):
        _collect(_TI_TARGET_RE)

    if not points and pattern in ("xi", "auto"):
        _collect(_XI_RE)

    if not points:
        raise FileNotFoundError(
            f"No constraint-point directories found in {root} "
            f"(tried pattern={pattern!r})."
        )

    # Sort by xi
    points.sort(key=lambda p: p.xi, reverse=reverse)
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
