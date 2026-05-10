"""Shared handler-factory helpers for agent task modules.

Convention: all analysis modules MUST be imported inside function bodies
(lazy import), never at module top level.  The ``_make_handler`` factory
enforces this structurally via ``importlib.import_module`` in the closure.
"""

from __future__ import annotations

import logging
from importlib import import_module
from pathlib import Path
from typing import Any, Callable

from ._core import TaskResult

logger = logging.getLogger(__name__)

# Type alias for summary extraction functions
SummaryExtractor = Callable[[Any, dict[str, Any]], dict[str, Any]]


def _make_handler(
    task_name: str,
    target_fn_path: str,
    summary_extractor: SummaryExtractor | None = None,
) -> Callable[[dict[str, Any]], TaskResult]:
    """Generate a pass-through handler that calls the target function.

    Exception handling is done by the ``dispatch()`` layer — handlers
    produced here do NOT wrap calls in try-except.
    """

    def handler(params: dict[str, Any]) -> TaskResult:
        module_path, fn_name = target_fn_path.rsplit(":", 1)
        fn = getattr(import_module(module_path), fn_name)

        result = fn(**params)

        outputs = _normalize_outputs(result)
        summary = summary_extractor(result, params) if summary_extractor else {}

        return TaskResult(
            success=True,
            task=task_name,
            outputs=outputs,
            summary=summary,
        )

    return handler


def _normalize_outputs(result: Any) -> dict[str, str]:
    """Normalize diverse return types to ``{name: path_string}``."""
    if isinstance(result, dict):
        return {k: str(v) for k, v in result.items()}
    if isinstance(result, Path):
        return {"output": str(result)}
    if isinstance(result, list):
        # list[Path] from batch scripts → enumerate as workdir_0, workdir_1, ...
        return {f"workdir_{i}": str(p) for i, p in enumerate(result)}
    # Dataclass reports with a `workdirs` tuple field (e.g. TIGenBatchReport)
    wds = getattr(result, "workdirs", None)
    if wds is not None and isinstance(wds, (tuple, list)):
        return {f"workdir_{i}": str(p) for i, p in enumerate(wds)}
    # SurfaceChargeResult and similar dataclasses with csv_path
    if hasattr(result, "csv_path"):
        out: dict[str, str] = {"csv": str(result.csv_path)}
        png = result.csv_path.parent / (result.csv_path.stem + ".png")
        if png.exists():
            out["png"] = str(png)
        return out
    logger.debug(
        "Unrecognized return type %s, outputs will be empty",
        type(result).__name__,
    )
    return {}
