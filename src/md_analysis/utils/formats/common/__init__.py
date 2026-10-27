"""Engine-neutral single-file parsers.

Modules:

- ``cube`` — Gaussian cube file reader (header + 3-D volumetric data),
  plus slab-averaged Hartree potential utilities. Used by both CP2K
  (cube output) and VASP (LOCPOT-as-cube) consumers.

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.formats.common.cube import (
        read_cube_header_and_values,
        slab_average_potential_ev,
    )
"""

from __future__ import annotations

__all__: list[str] = []
