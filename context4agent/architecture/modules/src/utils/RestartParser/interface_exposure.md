# utils.RestartParser — Interface Exposure

## Public API

### CellParser

| Function | Description |
|----------|-------------|
| `parse_abc_from_restart(path)` | Extract (a, b, c) cell lengths from CP2K `.restart` file |
| `parse_abc_from_md_inp(path)` | Extract (a, b, c) from `md.inp` ABC line |

### ColvarParser

| Symbol | Description |
|--------|-------------|
| `ColvarRestart` | Frozen dataclass: timestep_fs, time_start_fs, constraints list |
| `ConstraintInfo` | Frozen dataclass: target_au, target_growth_au, colvar metadata |
| `ColvarMDInfo` | Combined restart + log data; `from_paths(restart, log)` factory |
| `parse_colvar_restart(path)` | Parse `.restart` -> `ColvarRestart` |
| `parse_lagrange_mult_log(path)` | Parse `.LagrangeMultLog` -> `LagrangeMultLog` |

## Stability

Stable — used by slowgrowth, constrained_ti, and TIGen.
