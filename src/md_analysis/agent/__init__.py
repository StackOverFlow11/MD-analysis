"""Agent-friendly programmatic entry point for md_analysis.

Public API
----------
- :func:`dispatch`         — execute a named analysis task
- :func:`list_tasks`       — enumerate all registered tasks
- :func:`get_task_schema`  — JSON Schema for a task's parameters
- :class:`TaskResult`      — unified return type
- :class:`TaskDef`         — task registration descriptor

Extension API (Evolving)
------------------------
Extension authors declaring a new contract-backed task can import the
Tools-layer contract types directly from the (private) ``_contracts``
module::

    from md_analysis.agent._contracts import (
        FieldSpec, ExceptionMapping, TaskContract,
    )

These symbols are intentionally **not** re-exported from ``__all__``
while the contract shape is still evolving toward Resources / Prompts /
MCP-server layers.  See ``interface_exposure.md`` §2.5 for details.
"""

from ._core import TaskDef, TaskHandler, TaskResult, list_tasks, register
from ._dispatch import dispatch, get_task_schema

__all__ = [
    "dispatch",
    "list_tasks",
    "get_task_schema",
    "TaskResult",
    "TaskDef",
    "TaskHandler",
    "register",
]

# Trigger task registration (side-effect import)
from . import _handlers as _handlers  # noqa: F401, E402
