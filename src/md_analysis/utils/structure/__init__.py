"""Structure analysis utilities for `md_analysis`.

Geometric and chemical-semantics helpers operating on ASE Atoms or numpy
arrays (metal layer detection, water molecule topology, 1D periodic
clustering). Modules here:

- Operate on in-memory structural data
- Do NOT parse engine-specific files (those belong to
  ``md_analysis.utils.formats``)
- Do NOT discover directories or do filesystem walking (those belong to
  ``md_analysis.utils.io``)

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.structure.layer import detect_interface_layers
    from md_analysis.utils.structure.water import detect_water_molecule_indices
    from md_analysis.utils.structure.cluster import cluster_1d_periodic
"""

from __future__ import annotations

__all__: list[str] = []
