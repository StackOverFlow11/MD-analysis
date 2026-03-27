# utils.RestartParser — Implementation Guidelines

## Role

CP2K restart file and LagrangeMultLog parsing. Provides cell parameters and COLVAR constraint data for both slow-growth and constrained-TI workflows.

## File Organization

- `CellParser.py` — `parse_abc_from_restart()`, `parse_abc_from_md_inp()`: extract orthogonal cell lengths from `.restart` or `md.inp`
- `ColvarParser.py` — `parse_colvar_restart()`, `parse_lagrange_mult_log()`, `ColvarMDInfo.from_paths()`: parse COLVAR restart metadata and Lagrange multiplier time series

## Key Conventions

- All CV values in atomic units (CP2K default)
- Timestep in fs
- `TARGET_GROWTH` in per-a.u.-time (converted to per-step by consumers)
- `_safe_float()` converts CP2K overflow markers (`***`) to `np.nan`

## Dependencies

- numpy, re, pathlib (stdlib + numpy only, no ASE)
