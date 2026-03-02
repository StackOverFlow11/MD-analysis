# src.potential — Implementation Guidelines

## Layer dependency

- `src.potential` depends on `src.utils` (CubeParser, ClusterUtils, config)
- `src.potential` does NOT depend on `src.water`

## Module layout

Flat structure (no sub-packages):
- `config.py` — default output filenames
- `CenterPotential.py` — center slab potential + Fermi + electrode potential
- `PhiZProfile.py` — φ(z) plane-averaged profile analysis

## Cube file conventions

- CP2K V_HARTREE_CUBE: z is the fastest-running index → reshape as `(nx, ny, nz)`
- Units in cube files: Bohr (positions) and Hartree (values)
- Conversion constants are in `src.utils.config`

## Interface detection for slab centering

Uses `ClusterUtils.cluster_1d_periodic` + `find_largest_gap_periodic` to
detect metal layers from xyz trajectory z-coordinates, then finds the
water region midpoint as the slab center.
