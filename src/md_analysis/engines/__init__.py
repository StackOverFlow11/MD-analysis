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
"""

from __future__ import annotations

# Triggering imports register the default CP2K parser at package
# import time so callers can use ``infer_parser`` / ``get_parser("cp2k")``
# without an explicit registration step.
from .protocols import register_parser as _register_parser
from .cp2k import CP2KParser as _CP2KParser

_register_parser("cp2k", _CP2KParser)


__all__: list[str] = []
