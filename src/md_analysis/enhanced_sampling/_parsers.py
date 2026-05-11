"""Transitional shim — Phase 5 / Phase 6 migration target.

Phase 5a moved the canonical implementation to ``md_analysis.engines``:

- Protocol + registry + sniffer → ``md_analysis.engines.protocols``
- CP2K parser implementation     → ``md_analysis.engines.cp2k``
- Engine-neutral dataclasses     → ``md_analysis.engines.models``

This module re-exports the same names for backwards-compatible imports
so callers under ``enhanced_sampling`` keep working unchanged during the
transition.  Phase 6 will rewire those callers to import from
``md_analysis.engines`` directly and delete this shim.

Do not add new code here — extend ``md_analysis.engines.*`` instead.
"""

from __future__ import annotations

# Importing the package triggers the default CP2K parser registration
# in ``engines/__init__.py`` (see that module's top-level code).
from ..engines import protocols as _protocols  # noqa: F401  (side-effect import)

from ..engines.cp2k import CP2KParser
from ..engines.protocols import (
    ConstraintMDParser,
    ParserInferenceError,
    _REGISTRY,  # noqa: F401  (private — exposed transitionally; will be removed in Phase 6)
    get_parser,
    infer_parser,
    register_parser,
    resolve_parser,
)


__all__ = [
    "ConstraintMDParser",
    "ParserInferenceError",
    "CP2KParser",
    "register_parser",
    "get_parser",
    "infer_parser",
    "resolve_parser",
]
