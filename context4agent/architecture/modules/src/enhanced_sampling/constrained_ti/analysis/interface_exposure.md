# constrained_ti.analysis — Interface Exposure

## Public API

| Function | Module | Returns |
|----------|--------|---------|
| `analyze_autocorrelation(series, ...)` | autocorrelation.py | `AutocorrResult` |
| `analyze_block_average(series, ...)` | block_average.py | `BlockAverageResult` |
| `analyze_running_average(series, ...)` | running_average.py | `RunningAverageResult` |
| `analyze_geweke(series, ...)` | geweke.py | `GewekeResult` |

## Internal (used by autocorrelation + geweke)

| Function | Module | Description |
|----------|--------|-------------|
| `compute_acf(series)` | _acf_core.py | FFT-accelerated normalized ACF |
| `compute_iat(acf, n)` | _acf_core.py | Integrated autocorrelation time |
| `compute_sem_corrected(series, tau)` | _acf_core.py | Corrected standard error |

## Stability

Internal to constrained_ti — accessed via workflow.py orchestrator.
