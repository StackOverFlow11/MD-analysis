"""Low-level utilities package for `md_analysis`.

This package intentionally does NOT re-export symbols. Consumers should
import directly from the relevant submodule:

- physical constants        -> ``md_analysis.utils.constants``
- structure / layer parsers -> ``md_analysis.utils.structure.{layer,water,cluster}``
- file format parsers       -> ``md_analysis.utils.formats.{cube,bader,cp2k_cell,cp2k_colvar}``

Rationale: a centralised re-export hub must be kept in sync with every
submodule change and offers no additional information over a direct
import path. Direct submodule imports make the call sites self-documenting.
"""

from __future__ import annotations

__all__: list[str] = []
