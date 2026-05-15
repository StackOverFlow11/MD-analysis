# constrained_ti — Interface Exposure

## Public API (from __init__.py)

### Models & Exceptions

`ConstraintPointInput`, `ConstraintPointReport`, `TIReport`, `TIPointDefinition`,
`AutocorrResult`, `BlockAverageResult`, `RunningAverageResult`, `GewekeResult`,
`ConstantPotentialCorrection`, `ConstantPotentialResult`,
`ConvergenceError`, `InsufficientSamplingError`

### Workflow (workflow.py)

| Function | Description |
|----------|-------------|
| `analyze_standalone(series, ...)` | Single-point diagnostics (no TI context) |
| `analyze_ti(xi_values, lambda_list, dt, ...)` | Multi-point TI analysis -> TIReport |
| `standalone_diagnostics(restart, log, ...)` | Parse + analyze + plot + CSV |
| `write_convergence_csv(ti_report, ...)` | Per-point convergence CSV |
| `write_free_energy_csv(ti_report, ...)` | Free energy profile CSV |

### I/O (io.py)

Engine-agnostic — parsing delegated to `ConstraintMDParser` Protocol
exposed by `md_analysis.engines.protocols` (the legacy
`enhanced_sampling/_parsers.py` shim has been removed; CP2K
implementation lives in `md_analysis.engines.cp2k.CP2KParser` and is
auto-registered at engines package import).

| Function | Description |
|----------|-------------|
| `discover_ti_points(root, *, parser="auto", dir_filter=None, reverse=False, strict=False)` | Discover constraint-point dirs; `parser`="auto" sniffs registered parsers, `dir_filter` is None/glob/callable |
| `load_ti_series(point_defs)` | Parse Lagrange-multiplier series for each point (metadata is already cached on each `TIPointDefinition`) |

`TIPointDefinition` fields: `directory: Path`,
`parser: md_analysis.engines.protocols.ConstraintMDParser`,
`metadata: md_analysis.engines.models.ConstraintMetadata`. `xi` is a
property derived from `metadata.colvars.primary.target_au` — single
source of truth, never parsed from directory name.

### Plot (plot.py)

| Function | Description |
|----------|-------------|
| `plot_point_diagnostics(report, ...)` | 2x2 diagnostic PNG per point |
| `plot_free_energy_profile(ti_report, ...)` | Free energy profile PNG |

### Correction (correction.py)

| Function | Description |
|----------|-------------|
| `compute_constant_potential_correction(ti_report, point_defs, mapper, ...)` | Norskov correction |
| `write_corrected_free_energy_csv(result, ...)` | Corrected free energy CSV |
| `plot_corrected_free_energy_profile(result, ...)` | Corrected free energy PNG |

## Re-export

NOT re-exported from `md_analysis.__init__`. Import directly:
`from md_analysis.enhanced_sampling.constrained_ti import ...`
