"""Discovery and batch parsing for constrained-TI data.

Engine-agnostic: parsing is delegated to a :class:`ConstraintMDParser`.
This module only handles directory enumeration, sorting, and orchestration.
"""

from __future__ import annotations

import fnmatch
import logging
from pathlib import Path
from typing import Callable

import numpy as np

from .._parsers import (
    _REGISTRY,
    ConstraintMDParser,
    infer_parser,
    resolve_parser,
)
from .models import TIPointDefinition

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# dir_filter resolution
# ---------------------------------------------------------------------------

def _make_filter(
    dir_filter: str | Callable[[Path], bool] | None,
    parser: ConstraintMDParser,
) -> Callable[[Path], bool]:
    """Resolve a *dir_filter* argument into a callable predicate.

    Resolution rules:
      * ``None``     → use ``parser.is_constraint_directory`` (auto mode)
      * ``str``      → glob pattern matched against ``directory.name``
      * ``Callable`` → used as-is

    Note: ``str`` mode does *not* check that the directory contains the
    parser's expected files — the user is opting into name-based
    filtering and accepts that mismatched directories will fail later
    in :func:`load_ti_series`.  Use ``None`` if you want both name
    flexibility and automatic content validation.
    """
    if dir_filter is None:
        return parser.is_constraint_directory
    if isinstance(dir_filter, str):
        pattern = dir_filter

        def _glob(d: Path) -> bool:
            return d.is_dir() and fnmatch.fnmatch(d.name, pattern)

        return _glob
    return dir_filter


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def discover_ti_points(
    root_dir: Path,
    *,
    parser: ConstraintMDParser | str = "auto",
    dir_filter: str | Callable[[Path], bool] | None = None,
    reverse: bool = False,
    strict: bool = False,
) -> list[TIPointDefinition]:
    """Discover constraint-point directories under *root_dir*.

    Parameters
    ----------
    root_dir : Path
        Parent directory containing constraint-point subdirectories.
    parser : ConstraintMDParser | str, default ``"auto"``
        Parser to read each point.  ``"auto"`` triggers
        :func:`infer_parser` on the first matching candidate (or the
        first subdirectory if ``dir_filter`` is also auto).  Pass an
        instance (or the name of a registered parser) to skip sniffing.
    dir_filter : str | Callable | None, default ``None``
        How to decide which subdirectories are constraint points:

        * ``None`` — use ``parser.is_constraint_directory`` (recognises
          any directory containing the engine's expected files,
          regardless of name).
        * ``str``  — glob pattern matched against the directory name
          (e.g. ``"ti_target_*"``, ``"run_*"``).
        * Callable — custom predicate ``(Path) -> bool``.
    reverse : bool, default ``False``
        Sort by ξ descending (initial state = max ξ).
    strict : bool, default ``False``
        If ``False``, directories that match *dir_filter* but fail to
        parse (corrupt restart, missing log, etc.) are logged and
        skipped.  If ``True`` (agent-facing behaviour), the failure is
        re-raised so callers can surface it as ``file_not_found``.

    Returns
    -------
    list[TIPointDefinition]
        Sorted by ξ (ascending by default; descending if ``reverse``).
        Tiebreaker: directory name.

    Raises
    ------
    FileNotFoundError
        If *root_dir* doesn't exist, or no candidate directories yield
        valid points.  When ``strict=True``, also raised on per-point
        parse failures.
    ParserInferenceError
        If ``parser="auto"`` is requested but no registered parser
        recognises any candidate directory.
    """
    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory does not exist: {root}")

    candidates = sorted(d for d in root.iterdir() if d.is_dir())

    # ── Resolve parser ───────────────────────────────────────────────
    # Sniffing is for "auto" mode WITHOUT an explicit dir_filter — i.e.
    # the user wants us to figure out which engine produced these dirs.
    # If the user gave a name-based dir_filter, they have already
    # asserted the dirs are constraint points; we trust the assertion
    # and use the first registered parser as a sane default.  Pass an
    # explicit parser= to override.
    if parser == "auto":
        if dir_filter is None:
            if not candidates:
                raise FileNotFoundError(
                    f"No subdirectories under {root} to sniff a parser from."
                )
            resolved_parser: ConstraintMDParser | None = None
            for c in candidates:
                try:
                    resolved_parser = infer_parser(c)
                    break
                except Exception:
                    continue
            if resolved_parser is None:
                # Surface the original error against the first candidate.
                infer_parser(candidates[0])  # raises ParserInferenceError
            parser_obj: ConstraintMDParser = resolved_parser  # type: ignore[assignment]
        else:
            if not _REGISTRY:
                raise FileNotFoundError(
                    "No parser registered; cannot resolve parser='auto'."
                )
            # First-registered parser as default. CP2K is registered at
            # module import time, so this gives the historical CP2K
            # behaviour for free.
            parser_obj = next(iter(_REGISTRY.values()))()
    else:
        parser_obj = resolve_parser(parser)

    # ── Resolve dir_filter ───────────────────────────────────────────
    predicate = _make_filter(dir_filter, parser_obj)

    # ── Walk + parse metadata ────────────────────────────────────────
    points: list[TIPointDefinition] = []
    for d in candidates:
        if not predicate(d):
            continue
        try:
            metadata = parser_obj.parse_metadata(d)
        except Exception as e:
            if strict:
                raise FileNotFoundError(
                    f"Failed to parse metadata for {d.name!r}: {e}"
                ) from e
            logger.warning("Skipping %s: %s", d.name, e)
            continue
        points.append(
            TIPointDefinition(
                directory=d,
                parser=parser_obj,
                metadata=metadata,
            )
        )

    if not points:
        raise FileNotFoundError(
            f"No constraint-point directories found under {root} "
            f"(parser={parser_obj.name}, dir_filter={dir_filter!r})."
        )

    points.sort(key=lambda p: (p.xi, p.directory.name), reverse=reverse)
    return points


def load_ti_series(
    point_defs: list[TIPointDefinition],
) -> list[tuple[float, np.ndarray, float]]:
    """Parse the Lagrange-multiplier series for each point.

    Metadata (timestep, target, etc.) is already cached on each
    ``TIPointDefinition``; only the heavy λ(t) array is read here.
    Returns full (untrimmed) series — equilibration trimming is the
    workflow's responsibility.

    Returns
    -------
    list[tuple[float, np.ndarray, float]]
        ``(xi, lambda_series, dt_fs)`` for each point.
    """
    results = []
    for pdef in point_defs:
        log = pdef.parser.parse_lambda_series(pdef.directory)
        lambda_series = log.collective_shake
        dt_fs = float(pdef.metadata.timestep_fs)
        results.append((pdef.xi, lambda_series, dt_fs))
    return results
