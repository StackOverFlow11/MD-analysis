"""Unified dispatch, parameter coercion, and JSON Schema generation."""

from __future__ import annotations

import collections.abc
import inspect
import logging
import pathlib
import types
from importlib import import_module
from typing import Any, get_args, get_origin

from ..exceptions import MDAnalysisError
from ._core import (
    ERROR_ANALYSIS,
    ERROR_FILE_NOT_FOUND,
    ERROR_INTERNAL,
    ERROR_VALIDATION,
    TaskResult,
    get_task,
)

logger = logging.getLogger(__name__)

# ── Type mapping for JSON Schema ─────────────────────────────────────

_TYPE_MAP: dict[type, str] = {
    str: "string",
    int: "integer",
    float: "number",
    bool: "boolean",
}


# ── Public API ───────────────────────────────────────────────────────


def dispatch(task: str, params: dict[str, Any] | None = None) -> TaskResult:
    """Execute a named analysis task.

    Parameters
    ----------
    task : str
        Registered task name (e.g. ``"water_three_panel"``).
    params : dict, optional
        Task parameters.  Path-typed values accept plain strings.

    Returns
    -------
    TaskResult
    """
    # Resolve task (catch unknown task names gracefully)
    try:
        task_def = get_task(task)
    except ValueError as exc:
        return TaskResult(
            success=False,
            task=task,
            outputs={},
            summary={},
            error_type=ERROR_VALIDATION,
            errors=[str(exc)],
        )

    # Resolve target function and its signature once
    fn, sig = _resolve_fn(task_def.target_fn)
    coerced = _coerce_params(sig, params or {})

    # Layered exception handling
    try:
        return task_def.handler(coerced)
    except (ValueError, TypeError) as exc:
        logger.warning("Validation error in task %s: %s", task, exc)
        return TaskResult(
            success=False, task=task, outputs={}, summary={},
            error_type=ERROR_VALIDATION, errors=[str(exc)],
        )
    except FileNotFoundError as exc:
        logger.warning("File not found in task %s: %s", task, exc)
        return TaskResult(
            success=False, task=task, outputs={}, summary={},
            error_type=ERROR_FILE_NOT_FOUND, errors=[str(exc)],
        )
    except PermissionError as exc:
        logger.warning("Permission denied in task %s: %s", task, exc)
        return TaskResult(
            success=False, task=task, outputs={}, summary={},
            error_type=ERROR_FILE_NOT_FOUND, errors=[str(exc)],
        )
    except MDAnalysisError as exc:
        logger.error("Analysis error in task %s: %s", task, exc)
        return TaskResult(
            success=False, task=task, outputs={}, summary={},
            error_type=ERROR_ANALYSIS, errors=[str(exc)],
        )
    except Exception as exc:
        logger.error("Unexpected error in task %s", task, exc_info=True)
        return TaskResult(
            success=False, task=task, outputs={}, summary={},
            error_type=ERROR_INTERNAL,
            errors=[f"{type(exc).__name__}: {exc}"],
        )


def get_task_schema(task: str) -> dict[str, Any]:
    """Generate JSON Schema from the target function's signature.

    Uses ``typing.get_type_hints()`` to resolve stringified annotations
    (caused by ``from __future__ import annotations``).
    """
    task_def = get_task(task)
    fn, sig = _resolve_fn(task_def.target_fn)

    # get_type_hints resolves string annotations → real types
    try:
        hints = _get_type_hints_safe(fn)
    except Exception:
        hints = {}

    properties: dict[str, Any] = {}
    required: list[str] = []

    for pname, param in sig.parameters.items():
        if pname in ("self", "cls"):
            continue
        if param.kind == param.VAR_KEYWORD:
            continue
        if param.kind == param.VAR_POSITIONAL:
            continue

        annotation = hints.get(pname, param.annotation)
        prop = _annotation_to_schema(annotation)

        if pname in task_def.param_descriptions:
            prop["description"] = task_def.param_descriptions[pname]

        if pname in task_def.param_choices:
            prop["enum"] = task_def.param_choices[pname]

        if param.default is not inspect.Parameter.empty:
            if param.default is not None:
                prop["default"] = _serialize_default(param.default)
        elif not _is_optional(annotation):
            required.append(pname)

        properties[pname] = prop

    return {
        "name": task,
        "description": task_def.description,
        "parameters": {
            "type": "object",
            "properties": properties,
            "required": required,
        },
    }


# ── Internal helpers ─────────────────────────────────────────────────


def _resolve_fn(target_fn: str):
    """Import and return (callable, signature) for a ``module:name`` path."""
    module_path, fn_name = target_fn.rsplit(":", 1)
    mod = import_module(module_path)
    fn = getattr(mod, fn_name)
    sig = inspect.signature(fn)
    return fn, sig


def _get_type_hints_safe(fn):
    """``typing.get_type_hints`` with fallback for edge cases."""
    import typing
    return typing.get_type_hints(fn)


def _coerce_params(sig: inspect.Signature, params: dict[str, Any]) -> dict[str, Any]:
    """Convert JSON-compatible values to Python types expected by the target."""
    coerced: dict[str, Any] = {}

    try:
        hints = _get_type_hints_safe_from_sig(sig)
    except Exception:
        hints = {}

    for key, value in params.items():
        if key in sig.parameters:
            annotation = hints.get(key, sig.parameters[key].annotation)
            coerced[key] = _coerce_value(key, value, annotation)
        else:
            # Preserve unknown params — they may pass through **kwargs
            coerced[key] = value

    return coerced


def _get_type_hints_safe_from_sig(sig: inspect.Signature) -> dict[str, Any]:
    """Best-effort type resolution from signature parameters."""
    hints: dict[str, Any] = {}
    for pname, param in sig.parameters.items():
        ann = param.annotation
        if ann is not inspect.Parameter.empty:
            hints[pname] = ann
    return hints


def _coerce_value(key: str, value: Any, annotation: Any) -> Any:
    """Coerce a single parameter value based on its type annotation."""
    if value is None:
        return None

    # str → Path
    if _is_path_type(annotation):
        return pathlib.Path(value).resolve()

    origin = get_origin(annotation)

    # Handle Optional[X] — unwrap to X
    if _is_optional(annotation):
        inner = _unwrap_optional(annotation)
        return _coerce_value(key, value, inner)

    # list → tuple (for cell_abc: tuple[float, float, float])
    if origin is tuple and isinstance(value, list):
        return tuple(value)

    # list → set (for metal_elements: set[str])
    if origin is set and isinstance(value, list):
        return set(value)

    # list stays as list for Iterable[X] / Sequence[X] — no conversion needed
    # (Python list satisfies both Iterable and Sequence)

    return value


def _annotation_to_schema(annotation: Any) -> dict[str, Any]:
    """Convert a Python type annotation to a JSON Schema property."""
    # inspect.Parameter.empty or unresolved string
    if annotation is inspect.Parameter.empty or isinstance(annotation, str):
        return {"type": "string"}

    origin = get_origin(annotation)

    # Union types: X | Y (types.UnionType) or typing.Union[X, Y]
    if origin is types.UnionType or _is_typing_union(origin):
        args = get_args(annotation)
        non_none = [a for a in args if a is not type(None)]
        if len(non_none) == 1:
            return _annotation_to_schema(non_none[0])
        # Multi-type union (e.g. str | Path) — prefer str for JSON
        if str in non_none:
            return {"type": "string"}
        return {"type": "string"}

    # list[X]
    if origin is list:
        item_args = get_args(annotation)
        items = _annotation_to_schema(item_args[0]) if item_args else {}
        return {"type": "array", "items": items}

    # tuple[X, ...] — treat as array
    if origin is tuple:
        item_args = get_args(annotation)
        if item_args:
            items = _annotation_to_schema(item_args[0])
        else:
            items = {"type": "number"}
        return {"type": "array", "items": items}

    # set[X]
    if origin is set:
        item_args = get_args(annotation)
        items = _annotation_to_schema(item_args[0]) if item_args else {"type": "string"}
        return {"type": "array", "items": items, "uniqueItems": True}

    # Iterable[X] / Sequence[X]
    if origin in (collections.abc.Iterable, collections.abc.Sequence):
        item_args = get_args(annotation)
        items = _annotation_to_schema(item_args[0]) if item_args else {}
        return {"type": "array", "items": items}

    # Path
    if _is_path_type(annotation):
        return {"type": "string", "description": "file or directory path"}

    # Primitive types
    if annotation in _TYPE_MAP:
        return {"type": _TYPE_MAP[annotation]}

    # Fallback
    return {"type": "string"}


def _is_path_type(annotation: Any) -> bool:
    """Check if annotation is pathlib.Path (handles various import forms)."""
    try:
        return isinstance(annotation, type) and issubclass(annotation, pathlib.Path)
    except TypeError:
        return False


def _is_optional(annotation: Any) -> bool:
    """Check if annotation is Optional[X] (i.e. X | None)."""
    origin = get_origin(annotation)
    if origin is types.UnionType or _is_typing_union(origin):
        return type(None) in get_args(annotation)
    return False


def _is_typing_union(origin: Any) -> bool:
    """Check if origin is typing.Union."""
    import typing
    return origin is getattr(typing, "Union", None)


def _unwrap_optional(annotation: Any) -> Any:
    """Extract T from Optional[T]."""
    args = get_args(annotation)
    non_none = [a for a in args if a is not type(None)]
    return non_none[0] if len(non_none) == 1 else annotation


def _serialize_default(value: Any) -> Any:
    """Convert a default value to JSON-compatible form."""
    if isinstance(value, pathlib.Path):
        return str(value)
    if isinstance(value, (set, frozenset)):
        return sorted(value)
    if isinstance(value, tuple):
        return list(value)
    return value
