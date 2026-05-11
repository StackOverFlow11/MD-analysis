"""File-format parsers for `md_analysis`.

Single-file parsers for specific formats (Gaussian cube, VASP ACF, CP2K
restart / md.inp / colvar log, etc.). Modules here:

- Parse one file or one format family
- Return numpy arrays / frozen dataclasses
- Do NOT discover directories, walk filesystems, or pick which file to read
  (those concerns belong to ``md_analysis.utils.io``)
- Do NOT carry engine-level workflow logic (those belong to
  ``md_analysis.engines``)

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.formats.cube import (
        read_cube_header_and_values,
        slab_average_potential_ev,
    )
    from md_analysis.utils.formats.bader import load_bader_atoms
    from md_analysis.utils.formats.cp2k_cell import parse_abc_from_restart
    from md_analysis.utils.formats.cp2k_colvar import parse_colvar_restart
"""

from __future__ import annotations

__all__: list[str] = []
