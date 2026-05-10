"""Task registration aggregator.

Importing this module registers all built-in agent tasks via module-level
side effects in the task modules below.  The import order here dictates
the ``list_tasks()`` enumeration order.

Reload semantics
----------------
Test code occasionally does::

    from md_analysis.agent._core import _reset_registry
    _reset_registry()
    importlib.reload(md_analysis.agent._handlers)

to re-populate the registry from scratch.  After the split in
``10_split_agent_handlers_plan.md``, a bare ``reload(_handlers)`` would
only re-execute *this* aggregator; the ``from . import _tasks_*`` lines
then find the task modules already cached in ``sys.modules`` and would
**not** re-run their ``register(...)`` side effects.

To make reload re-populate the registry **without** double-executing
each task module on a fresh first import, we distinguish the two cases
via ``sys.modules``: on reload the task modules are already cached and
we explicitly ``importlib.reload()`` them; on first import they are
absent, the subsequent ``from . import`` lines load them exactly once.
"""

import importlib as _importlib
import sys as _sys

_PKG = __name__.rsplit(".", 1)[0]
_TASK_MODULE_NAMES = (
    f"{_PKG}._tasks_legacy",
    f"{_PKG}._tasks_ti",
    f"{_PKG}._tasks_scripts",
)

# Reload path: only triggers when the task modules are already cached,
# i.e. this aggregator is being reloaded.  First-import path skips this
# loop entirely, so each task module body runs exactly once via the
# ``from . import`` statements below.
for _name in _TASK_MODULE_NAMES:
    if _name in _sys.modules:
        _importlib.reload(_sys.modules[_name])

from . import _tasks_legacy  # noqa: E402, F401
from . import _tasks_ti  # noqa: E402, F401
from . import _tasks_scripts  # noqa: E402, F401

del _importlib, _sys, _name, _PKG, _TASK_MODULE_NAMES
