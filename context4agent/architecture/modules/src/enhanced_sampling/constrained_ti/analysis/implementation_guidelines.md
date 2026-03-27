# constrained_ti.analysis — Implementation Guidelines

## Role

Four independent analysis engines for TI convergence diagnostics. All are pure functions (no logging, no I/O, no side effects). Only depend on numpy.

## Engines

### _acf_core.py (shared)
- `compute_acf(series)` — FFT-accelerated normalized ACF, C(0)=1
- `compute_iat(acf, n)` — Sokal (1997) self-consistent IAT: M >= alpha*tau, hard cap M <= N*0.5
- `compute_sem_corrected(series, tau)` — SEM = sigma * sqrt(2*tau/N)

### autocorrelation.py (Step 2)
- `analyze_autocorrelation(series, *, sem_max, ...)` -> `AutocorrResult`
- Wraps _acf_core + pass/fail: N_eff >= 50, SEM <= SEM_max

### block_average.py (Step 3)
- `analyze_block_average(series, *, sem_max, min_blocks=4, n_consecutive=2)` -> `BlockAverageResult`
- F&P (1989): pow2 block sizes, delta_SEM = SEM(B) / sqrt(2*(n_b-1))
- Plateau: n_consecutive levels where SEM increase < delta_SEM

### running_average.py (Step 1)
- `analyze_running_average(series, *, sem, drift_factor=3.0)` -> `RunningAverageResult`
- Drift D = max-min of cumulative mean over [N/2, N]; pass if D < factor*SEM

### geweke.py (Step 4)
- `analyze_geweke(series, *, f_a=0.1, f_b=0.5, z_crit=1.96)` -> `GewekeResult`
- z-test on front vs rear segments using spectral variance from ACF
- Reliability flag: N_eff_a >= min_neff_subseries

## Constraints

- No engine imports another engine
- No logging, no file I/O
- Only numpy dependency
