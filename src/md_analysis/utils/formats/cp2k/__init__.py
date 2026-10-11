"""CP2K single-file parsers.

Modules:

- ``cell`` — parse orthorhombic ``ABC`` cell parameters from a CP2K
  ``.restart`` or ``md.inp`` file.
- ``colvar`` — parse COLVAR restart metadata and ``LagrangeMultLog``
  time series for slow-growth / constrained-MD.
- ``inp`` — parse and generate CP2K input files via
  ``cp2k-input-tools`` (cit tree); the only module in the package that
  imports cit and ``pint`` (both lazily).  Provides the tree primitives
  shared by the ``cell`` / ``colvar`` readers and by the
  ``md_analysis.scripts`` generators ("parse -> tree edit ->
  generate").
- ``stdout`` — parse Fermi-level series from CP2K ``md.out`` /
  single-point ``sp.out``.
- ``xyz`` — read selected steps from a CP2K trajectory ``-pos-1.xyz``.

Each submodule returns frozen dataclasses, numpy arrays, or (``inp``)
plain dict trees; none of them discover directories or pick which file
to read (those concerns belong to ``md_analysis.utils.io``).

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.formats.cp2k.cell import parse_abc_from_restart
    from md_analysis.utils.formats.cp2k.colvar import parse_colvar_restart
"""

from __future__ import annotations

__all__: list[str] = []
