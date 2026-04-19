"""Core data structures: TaskResult, TaskHandler Protocol, TaskDef, registry."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

if TYPE_CHECKING:
    from ._contracts import TaskContract

logger = logging.getLogger(__name__)

# ── Error type constants ─────────────────────────────────────────────

ERROR_VALIDATION = "validation"
ERROR_FILE_NOT_FOUND = "file_not_found"
ERROR_ANALYSIS = "analysis"
ERROR_INTERNAL = "internal"


# ── TaskResult ───────────────────────────────────────────────────────

@dataclass(frozen=True)
class TaskResult:
    """Unified return type for all agent tasks."""

    success: bool
    task: str
    outputs: dict[str, str]
    summary: dict[str, Any]
    error_type: str | None = None
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable dict."""
        return {
            "success": self.success,
            "task": self.task,
            "outputs": self.outputs,
            "summary": self.summary,
            "error_type": self.error_type,
            "errors": self.errors,
            "warnings": self.warnings,
        }


# ── TaskHandler Protocol ─────────────────────────────────────────────

@runtime_checkable
class TaskHandler(Protocol):
    """Callable signature that all task handlers must satisfy."""

    def __call__(self, params: dict[str, Any]) -> TaskResult: ...


# ── TaskDef ──────────────────────────────────────────────────────────

@dataclass(frozen=True)
class TaskDef:
    """Descriptor for a dispatchable analysis task.

    Two schema/coercion sources, in priority order:

    1. **``contract``** (preferred, Tools-layer ``TaskContract``) — when
       present, ``get_task_schema()`` and ``dispatch()`` derive schema and
       parameter coercion from the contract, ignoring ``target_fn`` signature.
    2. **``target_fn`` signature + ``param_descriptions`` / ``param_choices``**
       (legacy fallback) — used when ``contract`` is ``None``.  Omitting
       ``param_descriptions`` / ``param_choices`` only reduces schema
       readability, never causes runtime errors.

    ``reference_fn`` is an optional dotted path ("module.path:function_name")
    used by static tests to verify that ``contract.inputs`` aligns with a
    real function signature.  Composite handlers whose inner ``target_fn``
    has a different signature should set ``reference_fn`` to the wrapper
    that matches the contract.
    """

    name: str
    category: str
    description: str
    handler: TaskHandler
    target_fn: str  # dotted path "module.path:function_name"
    cli_codes: tuple[str, ...] = ()
    param_descriptions: dict[str, str] = field(default_factory=dict)
    param_choices: dict[str, list[str]] = field(default_factory=dict)
    contract: "TaskContract | None" = None
    reference_fn: str | None = None


# ── Registry ─────────────────────────────────────────────────────────

_TASK_REGISTRY: dict[str, TaskDef] = {}


def register(task_def: TaskDef) -> None:
    """Register a task definition (idempotent)."""
    _TASK_REGISTRY[task_def.name] = task_def


def get_task(name: str) -> TaskDef:
    """Look up a registered task by name.

    Raises ``ValueError`` for unknown tasks.
    """
    if name not in _TASK_REGISTRY:
        available = sorted(_TASK_REGISTRY.keys())
        raise ValueError(f"Unknown task: {name!r}. Available: {available}")
    return _TASK_REGISTRY[name]


def list_tasks() -> list[dict[str, Any]]:
    """Return a summary list of all registered tasks."""
    return [
        {
            "name": t.name,
            "category": t.category,
            "description": t.description,
            "cli_codes": t.cli_codes,
        }
        for t in _TASK_REGISTRY.values()
    ]


def _reset_registry() -> None:
    """Clear the registry (for testing only)."""
    _TASK_REGISTRY.clear()
