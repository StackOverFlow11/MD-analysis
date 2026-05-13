"""Bader charge analysis parsers.

Modules:

- ``acf`` — parse ACF.dat from the Bader code and attach charges to an
  ASE ``Atoms`` object (``load_bader_atoms``).
- ``potcar`` — extract ``(element, ZVAL)`` pairs from a VASP POTCAR.
- ``_errors`` — shared ``BaderParseError`` exception.

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.formats.bader.acf import load_bader_atoms
"""

from __future__ import annotations

__all__: list[str] = []
