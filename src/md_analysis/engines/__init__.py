"""Engine facade layer for `md_analysis`.

Provides engine-neutral data models and per-engine adapter modules:

- ``engines.models``    — engine-neutral frozen dataclasses (constraint
  metadata, λ(t) series, potential frames, ...) returned by engine
  parsers
- ``engines.protocols`` — ``ConstraintMDParser`` Protocol + parser
  registry / discovery (``register_parser``, ``infer_parser``, ...)
- ``engines.cp2k``      — CP2K adapter; reads CP2K restart /
  LagrangeMultLog / cube / md.out output into the neutral models
- ``engines.vasp``      — VASP adapter placeholder (not auto-registered,
  explicit ``NotImplementedError`` until populated)

Dependency direction:

    engines/ ─► utils/formats, utils/constants, utils/io, exceptions

``engines`` MUST NOT import from ``electrochemical``, ``water``,
``enhanced_sampling``, ``scripts``, or ``cli`` (those are upper layers
that consume engine output).

The public API surface (Protocol + registry + CP2K parser + neutral
dataclasses) is re-exported here so callers can write::

    from md_analysis.engines import get_parser, infer_parser, PotentialFrame

without reaching into the ``engines.protocols`` / ``engines.models``
submodules.
"""

from __future__ import annotations

from .cp2k import CP2KParser
from .models import PotentialFrame
from .protocols import (
    ConstraintMDParser,
    ParserInferenceError,
    get_parser,
    infer_parser,
    register_parser,
    resolve_parser,
)

# Default engine registration (import-time side-effect): every public
# entry point (``infer_parser``, ``get_parser("cp2k")`` …) sees the CP2K
# adapter as soon as the package is imported.
register_parser("cp2k", CP2KParser)


__all__ = [
    "ConstraintMDParser",
    "ParserInferenceError",
    "CP2KParser",
    "PotentialFrame",
    "register_parser",
    "get_parser",
    "infer_parser",
    "resolve_parser",
]
