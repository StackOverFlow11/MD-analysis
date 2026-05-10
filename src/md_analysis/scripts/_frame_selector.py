"""Unified trajectory frame selection for batch Gen scripts.

Provides a common ``FrameSelection`` dataclass and iterators used by
BaderGen / PotentialGen / SpGen to slice MD trajectories by either
frame index or simulation time (fs).

TIGen is NOT covered here because its "time mode" maps discrete target
times to nearest CV values (target-based), not range-based slicing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

from ase import Atoms
from ase.io import iread

from ..exceptions import MDAnalysisError

logger = logging.getLogger(__name__)


class FrameSelectionError(MDAnalysisError):
    """Raised when frame selection parameters or trajectory metadata are invalid."""


FrameMode = Literal["index", "time"]


@dataclass(frozen=True)
class FrameSelection:
    """Unified trajectory frame selection specification.

    Two mutually exclusive modes controlled by ``mode``:

    **Index mode** (``mode="index"``, default):
        Select frames at 0-based indices ``frame_start``, ``frame_start+step``,
        ... up to (and excluding) ``frame_end``. Matches the historic batch
        script behavior.

    **Time mode** (``mode="time"``):
        Requires all three of ``time_start_fs``, ``time_end_fs``, and
        ``time_step_fs``. Frames are selected greedily: at each target time
        ``t_k = t_start + k * t_step``, the first trajectory frame with
        ``time >= t_k`` is yielded (provided ``time <= t_end``), then the
        target advances by ``t_step``. Requires ``atoms.info["time"]``
        metadata in the XYZ comment lines (CP2K default).

    The ``mode`` field is an explicit discriminator (not inferred from
    which parameters are ``None``) so that JSON Schema auto-generation
    for agent dispatch cleanly exposes an enum choice.
    """

    mode: FrameMode = "index"

    # Index mode
    frame_start: int = 0
    frame_end: int | None = None
    frame_step: int = 1

    # Time mode — all three must be set together
    time_start_fs: float | None = None
    time_end_fs: float | None = None
    time_step_fs: float | None = None

    def __post_init__(self) -> None:
        if self.mode not in ("index", "time"):
            raise FrameSelectionError(
                f"Invalid mode {self.mode!r}; expected 'index' or 'time'"
            )

        if self.frame_step < 1:
            raise FrameSelectionError(
                f"frame_step must be >= 1, got {self.frame_step}"
            )

        if self.mode == "time":
            time_vals = (self.time_start_fs, self.time_end_fs, self.time_step_fs)
            if any(v is None for v in time_vals):
                raise FrameSelectionError(
                    "Time mode requires all of time_start_fs, time_end_fs, "
                    "time_step_fs"
                )
            if self.time_step_fs <= 0:
                raise FrameSelectionError(
                    f"time_step_fs must be > 0, got {self.time_step_fs}"
                )
            if self.time_start_fs > self.time_end_fs:
                raise FrameSelectionError(
                    f"time_start_fs ({self.time_start_fs}) must be <= "
                    f"time_end_fs ({self.time_end_fs})"
                )


def iter_selected_frames(
    xyz_path: str | Path,
    selection: FrameSelection,
) -> Iterator[tuple[int, Atoms]]:
    """Iterate selected ``(frame_idx, atoms)`` pairs from a CP2K XYZ trajectory.

    Parameters
    ----------
    xyz_path : str or Path
        CP2K XYZ trajectory file.
    selection : FrameSelection

    Yields
    ------
    tuple[int, Atoms]
        0-based frame index and the corresponding ase ``Atoms`` object.

    Raises
    ------
    FrameSelectionError
        If ``mode="time"`` and any frame lacks ``atoms.info["time"]``.
    """
    xyz_path = Path(xyz_path)

    if selection.mode == "index":
        yield from _iter_by_index(xyz_path, selection)
    else:
        yield from _iter_by_time(xyz_path, selection)


def _iter_by_index(
    xyz_path: Path, selection: FrameSelection
) -> Iterator[tuple[int, Atoms]]:
    next_yield = selection.frame_start
    for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
        if selection.frame_end is not None and idx >= selection.frame_end:
            break
        if idx == next_yield:
            yield idx, atoms
            next_yield += selection.frame_step


def _iter_by_time(
    xyz_path: Path, selection: FrameSelection
) -> Iterator[tuple[int, Atoms]]:
    t_start = float(selection.time_start_fs)   # type: ignore[arg-type]
    t_end = float(selection.time_end_fs)       # type: ignore[arg-type]
    t_step = float(selection.time_step_fs)     # type: ignore[arg-type]

    next_target = t_start

    for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
        time_fs = _require_time(atoms, idx)

        # Past the end of the requested interval → done
        if time_fs > t_end:
            break
        # Before the interval → skip
        if time_fs < t_start:
            continue

        # Greedy: yield first frame with time >= next_target
        if time_fs >= next_target:
            yield idx, atoms
            next_target = time_fs + t_step


def resolve_single_frame(
    xyz_path: str | Path,
    *,
    mode: FrameMode = "index",
    frame: int = 0,
    time_fs: float | None = None,
    time_tol_fs: float = 1e-6,
) -> tuple[int, Atoms, list[str]]:
    """Locate a single frame for ``generate_*_workdir`` single-frame commands.

    Parameters
    ----------
    xyz_path : str or Path
        CP2K XYZ trajectory file.
    mode : {"index", "time"}
        Selection mode.
    frame : int
        0-based frame index (used when ``mode="index"``).
    time_fs : float or None
        Target simulation time in fs (required when ``mode="time"``).
    time_tol_fs : float
        Tolerance for "exact" time match. If the nearest frame's time
        differs from ``time_fs`` by more than this value, a warning is
        appended to the returned warnings list.

    Returns
    -------
    tuple[int, Atoms, list[str]]
        ``(frame_index, atoms, warnings)``. The caller is responsible for
        surfacing warnings to the user (CLI: ``logger.warning``;
        agent: merge into ``TaskResult.warnings``).

    Raises
    ------
    FrameSelectionError
        If inputs are invalid or the trajectory lacks required metadata.
    FileNotFoundError
        If the trajectory file does not exist (from ase.io).
    """
    xyz_path = Path(xyz_path)

    if mode == "index":
        warnings: list[str] = []
        for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
            if idx == frame:
                return idx, atoms, warnings
        raise FrameSelectionError(
            f"Frame index {frame} not found in {xyz_path} "
            "(trajectory shorter than requested)"
        )

    if mode == "time":
        if time_fs is None:
            raise FrameSelectionError(
                "Time mode requires time_fs to be set"
            )
        return _find_nearest_by_time(xyz_path, float(time_fs), time_tol_fs)

    raise FrameSelectionError(
        f"Invalid mode {mode!r}; expected 'index' or 'time'"
    )


def _find_nearest_by_time(
    xyz_path: Path,
    target_fs: float,
    tol_fs: float,
) -> tuple[int, Atoms, list[str]]:
    """Find the frame whose ``atoms.info["time"]`` is closest to *target_fs*.

    Assumes time is monotonically non-decreasing; breaks early once we have
    passed the target and the gap starts growing.
    """
    best_idx: int | None = None
    best_atoms: Atoms | None = None
    best_time: float | None = None
    best_diff = float("inf")

    for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
        time_fs = _require_time(atoms, idx)
        diff = abs(time_fs - target_fs)

        if diff < best_diff:
            best_idx = idx
            best_atoms = atoms
            best_time = time_fs
            best_diff = diff
        elif time_fs > target_fs:
            # Past the target and gap is growing → done (monotonic assumption)
            break

    if best_idx is None or best_atoms is None or best_time is None:
        raise FrameSelectionError(
            f"No frames found in {xyz_path}"
        )

    warnings: list[str] = []
    if best_diff > tol_fs:
        warnings.append(
            f"No frame at exactly t={target_fs} fs; using nearest frame "
            f"index={best_idx} with t={best_time:.4f} fs "
            f"(Δt={best_time - target_fs:+.4f} fs)"
        )

    return best_idx, best_atoms, warnings


def _require_time(atoms: Atoms, idx: int) -> float:
    """Extract ``atoms.info['time']`` or raise ``FrameSelectionError``."""
    time_fs = atoms.info.get("time")
    if time_fs is None:
        raise FrameSelectionError(
            f"Frame {idx} is missing atoms.info['time']; "
            "time-based selection requires CP2K XYZ metadata "
            "(i = ..., time = ... , ...)"
        )
    return float(time_fs)
