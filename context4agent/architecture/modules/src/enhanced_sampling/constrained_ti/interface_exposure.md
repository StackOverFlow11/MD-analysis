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

| Function | Description |
|----------|-------------|
| `discover_ti_points(root_dir, *, pattern, reverse)` | Scan for ti_target_*/xi_* dirs |
| `load_ti_series(point_defs)` | Parse restart + log for each point |

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
