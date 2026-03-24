# Constrained TI Convergence Analysis — Architecture Design

## 0. Context

### 0.1 What this module does

Implements the four-step convergence diagnostic pipeline for **equilibrium constrained thermodynamic integration** (Blue Moon ensemble). Given a set of constraint points along a collective variable ξ, each with a Lagrange multiplier time series λ(t), the module:

1. Checks cumulative-mean stability (running average drift)
2. Estimates effective sample count via autocorrelation analysis
3. Determines the true standard error via block averaging (Flyvbjerg–Petersen)
4. Tests stationarity via the Geweke diagnostic

Then integrates the per-point free-energy gradients (梯形 rule, non-uniform spacing) with rigorous error propagation to yield ΔA ± σ_A.

The mathematical specification lives in `convergence_criteria.md` (repository root). This document covers **software architecture only**.

### 0.2 Relationship to existing code

```
enhanced_sampling/
├── slowgrowth/           # non-equilibrium slow-growth TI (existing)
│   ├── SlowGrowth.py     #   data structures
│   ├── SlowGrowthPlot.py #   plotting + CSV + unified entry point
│   └── config.py         #   filename constants
│
└── constrained_ti/       # equilibrium constrained TI (NEW)
    └── ...               #   this document describes the layout below
```

**Slow growth** = single continuous trajectory, CV dragged at constant rate, cumulative integral → ΔA.
**Constrained TI** = K independent fixed-ξ simulations, ensemble average at each point, numerical quadrature → ΔA.

The two share the parser layer (`ColvarParser`) but diverge entirely in analysis logic.

### 0.3 Design goals

1. **Separation of concerns** — each statistical test is an independent pure function; the orchestrator is the only module that knows execution order.
2. **Consistency with existing codebase** — frozen dataclasses, `config.py` for defaults, `_`-prefixed internals, lazy imports in CLI, same test pyramid (unit with synthetic data, integration with real data).
3. **All magic numbers configurable** — every threshold from `convergence_criteria.md` is exposed as a keyword argument with a default from `config.py`.
4. **NumPy-only core** — no dependencies beyond numpy for the analysis engines. Matplotlib is imported lazily in the plotting module only.

---

## 1. File Structure

```
src/md_analysis/enhanced_sampling/constrained_ti/
├── __init__.py                # re-exports public API
├── config.py                  # default thresholds + output filename constants
├── models.py                  # frozen dataclasses + exception hierarchy
│
├── analysis/                  # four independent analysis engines
│   ├── __init__.py            #   re-exports engine entry points
│   ├── _acf_core.py           #   shared ACF computation (FFT-accelerated)
│   ├── _arctan_fit.py         #   arctan SEM(B) extrapolation (pure function, WLS fit)
│   ├── autocorrelation.py     #   Step 2: τ_corr, N_eff, SEM_auto
│   ├── block_average.py       #   Step 3: SEM(B) plateau detection + arctan extrapolation
│   ├── running_average.py     #   Step 1: cumulative-mean drift check
│   └── geweke.py              #   Step 4: stationarity test
│
├── integration.py             # trapezoid weights, error propagation, optimal allocation
├── workflow.py                # orchestrator + standalone_diagnostics + CSV exports
├── plot.py                    # pure plotting (diagnostic PNGs + free-energy profile)
└── io.py                      # multi-point file discovery + batch parsing (pure I/O)
```

**Module responsibility boundaries** (one-way dependency constraints):
- `workflow.py` may import `plot.py`; `plot.py` must **never** import `workflow.py`.
- `io.py` must **never** import `integration.py`; weight computation belongs in `workflow.py`.
- No analysis engine imports another analysis engine.

### Corresponding test layout

```
test/
├── unit/enhanced_sampling/constrained_ti/
│   ├── test_acf_core.py
│   ├── test_autocorrelation.py
│   ├── test_block_average.py
│   ├── test_running_average.py
│   ├── test_geweke.py
│   ├── test_integration.py
│   ├── test_workflow.py
│   └── test_models.py
└── integration/enhanced_sampling/
    └── test_constrained_ti.py       # uses data_example/ti/ real data
```

### CLI extension

```
cli/_enhanced_sampling.py      # add new command classes (303, 304, 305) alongside existing 301, 302
```

### `main.py` extension

```python
# new entry point, parallel to run_water_analysis / run_potential_analysis
def run_constrained_ti_analysis(...) -> dict[str, Path]: ...
```

---

## 2. Data Structures (`models.py`)

All dataclasses are `frozen=True`, matching the `Slowgrowth` / `SlowgrowthFull` convention.

### 2.1 Per-engine result dataclasses

```python
@dataclass(frozen=True)
class RunningAverageResult:
    """Output of Step 1: cumulative-mean drift analysis."""
    running_mean: np.ndarray     # shape (N,) — λ̄(n) curve for plotting
    drift_D: float               # max - min of λ̄(n) over [N/2, N]
    drift_limit: float           # drift_factor × SEM (the comparison target)
    passed: bool

@dataclass(frozen=True)
class AutocorrResult:
    """Output of Step 2: autocorrelation analysis."""
    acf: np.ndarray              # C(j), truncated at self-consistent cutoff M
    tau_corr: float              # integrated autocorrelation time (frames); see §4.0 for definition
    n_eff: float                 # effective independent sample count
    sem_auto: float              # σ_λ √(2τ/N)
    passed_neff: bool            # N_eff ≥ N_EFF_MIN
    passed_sem: bool | None      # SEM_auto ≤ SEM_max; None if sem_max not provided
    t_min_frames: int | None     # if N_eff < 50, minimum frames needed; else None (meaning "sufficient")

@dataclass(frozen=True)
class ArctanFitResult:
    """Result of arctan extrapolation on the SEM(B) curve (empirical model)."""
    sem_asymptote: float         # A · π/2 — extrapolated asymptotic SEM
    A: float                     # amplitude parameter
    B: float                     # rate parameter
    r2: float                    # goodness-of-fit R²
    reliable: bool               # R² ≥ threshold AND A>0, B>0 AND pcov non-singular
    tau_corr_implied: float      # back-calculated τ = N·SEM²/(2·σ²)
    fit_curve: np.ndarray | None # fitted values at input block_sizes (for plotting)

@dataclass(frozen=True)
class BlockAverageResult:
    """Output of Step 3: block-averaging (Flyvbjerg–Petersen)."""
    block_sizes: np.ndarray      # B values tested
    sem_curve: np.ndarray        # SEM(B) for each B
    sem_plateau: float           # s̄ from plateau detection
    sem_at_max_B: float          # SEM at largest B (for fallback transparency)
    plateau_rtol: float          # δ_s / s̄
    plateau_reached: bool        # plateau_rtol < PLATEAU_RTOL_MAX
    cross_valid_ok: bool         # |SEM_block − SEM_auto| / SEM_block < threshold
    passed: bool | None          # plateau_reached AND SEM_plateau ≤ SEM_max; None if sem_max not provided
    arctan: ArctanFitResult | None = None  # None = not attempted (scipy unavailable / too few points)

@dataclass(frozen=True)
class GewekeResult:
    """Output of Step 4: Geweke stationarity test."""
    z: float                     # test statistic
    mean_a: float                # front-segment mean
    mean_b: float                # rear-segment mean
    n_eff_a: float | None        # effective samples in front segment (None if unreliable)
    passed: bool                 # |z| < z_crit
    reliable: bool               # False if front/rear sub-series too short for spectral variance
```

### 2.2 Per-point composite structures

```python
@dataclass(frozen=True)
class ConstraintPointInput:
    """Everything needed to diagnose one constraint point.

    Supports two modes:
    - **TI context** (multi-point): weight and sem_max are computed from
      trapezoid quadrature and ε_tol; pass/fail is fully determined.
    - **Standalone** (single-point): weight and sem_max are None;
      all four diagnostics still run but pass/fail for the SEM threshold
      is skipped (reported as "N/A — no precision target").
    """
    xi: float                    # CV constraint value (or user label)
    lambda_series: np.ndarray    # shape (N,) — Lagrange multiplier time series (post-equilibration)
    dt: float                    # frame interval (fs)
    weight: float | None         # trapezoid quadrature weight w_k (None in standalone mode)
    sem_max: float | None        # precision target SEM_{max,k} (None in standalone mode)
    point_index: int | None      # 0-based index into K points; None in standalone mode

@dataclass(frozen=True)
class ConstraintPointReport:
    """Complete diagnostic output for one constraint point."""
    xi: float
    point_index: int
    lambda_mean: float           # λ̄ = dA/dξ estimate at this point
    sigma_lambda: float          # sample standard deviation of λ

    # Engine results (carried verbatim for downstream inspection)
    autocorr: AutocorrResult
    block_avg: BlockAverageResult
    running_avg: RunningAverageResult
    geweke: GewekeResult

    # Summary
    sem_final: float             # 3-tier: arctan (primary) → plateau → SEM_auto (fallback)
    sem_max: float | None        # precision target for reference (None in standalone)
    passed: bool | None          # all four steps passed; None if sem_max not set
    failure_reasons: tuple[str, ...]   # empty if passed
```

### 2.3 Multi-point report

```python
@dataclass(frozen=True)
class TIReport:
    """Full thermodynamic-integration convergence report."""
    point_reports: tuple[ConstraintPointReport, ...]
    xi_values: np.ndarray        # shape (K,)
    weights: np.ndarray          # shape (K,) — trapezoid weights
    forces: np.ndarray           # shape (K,) — λ̄ at each point
    force_errors: np.ndarray     # shape (K,) — SEM_final at each point

    delta_A: float               # integrated free-energy difference
    sigma_A: float               # propagated statistical error
    epsilon_tol: float           # user-specified tolerance

    all_passed: bool             # every point passed
    failing_indices: tuple[int, ...]  # 0-based indices into point_reports

    # Optional: optimal allocation suggestion
    suggested_time_ratios: np.ndarray | None  # shape (K,) or None
```

### 2.4 I/O data structures

```python
@dataclass(frozen=True)
class TIPointDefinition:
    """Discovered constraint point file paths (output of io.discover_ti_points)."""
    xi: float                    # CV constraint value parsed from directory name
    restart_path: Path           # path to *.restart file
    log_path: Path               # path to *.LagrangeMultLog file
```

### 2.5 Exception hierarchy

Defined in `models.py` alongside the data contracts (part of the module's public interface).

```python
class ConvergenceError(MDAnalysisError):
    """Base for constrained-TI convergence errors."""

class InsufficientSamplingError(ConvergenceError):
    """Raised when sampling is too short for any reliable estimate (N < N_MIN)."""
```

---

## 3. Configuration (`config.py`)

Follows the pattern of `slowgrowth/config.py` (filename constants) and `utils/config.py` (physical constants). Statistical thresholds specific to constrained TI are defined locally here (not promoted to `utils/config.py`).

```python
"""Default thresholds for constrained-TI convergence diagnostics.

Every threshold corresponds to a boxed equation in convergence_criteria.md.
All are exposed as keyword arguments in the analysis functions; the values
here are the recommended defaults.

Unit convention for epsilon_tol:
  User-facing API accepts eV (DEFAULT_EPSILON_TOL_EV).
  workflow.analyze_ti converts to Hartree internally before computing
  SEM targets, so that SEM_max has the same unit as lambda_series (a.u.).
  The conversion chain: user eV → EV_TO_HARTREE → Hartree (a.u.).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

# Hard minimum: raise InsufficientSamplingError if N < this.
DEFAULT_N_MIN: int = 10

# Tiered warning thresholds for short series.
# N < N_WARN_UNRELIABLE: strong warning — ACF/block results may be unreliable.
# N < N_WARN_SHORT: mild warning — results should be interpreted with caution.
DEFAULT_N_WARN_UNRELIABLE: int = 100
DEFAULT_N_WARN_SHORT: int = 500

# Maximum fraction of NaN values allowed before raising InsufficientSamplingError.
DEFAULT_NAN_FRACTION_MAX: float = 0.1

# ---------------------------------------------------------------------------
# Step 2: Autocorrelation
# ---------------------------------------------------------------------------

# Self-consistent truncation multiplier (§2.2: M ≥ α · τ_corr).
# Reference: Sokal (1997), also default in ALPS and emcee.
DEFAULT_ACF_ALPHA: int = 5

# Hard upper limit for ACF truncation window to prevent noise-dominated estimates.
# M_max = N // 2 is enforced regardless of the self-consistent criterion.
# This protects against unreliable IAT estimates when τ_corr > N/10.
DEFAULT_ACF_M_MAX_FRACTION: float = 0.5

# Minimum effective independent samples (§2.5: N_eff ≥ 50).
DEFAULT_NEFF_MIN: int = 50

# ---------------------------------------------------------------------------
# Step 3: Block averaging
# ---------------------------------------------------------------------------

# Maximum relative spread for plateau detection (§3.3: δ_s / s̄ < 0.2).
DEFAULT_PLATEAU_RTOL: float = 0.2

# Minimum number of blocks at largest B (§3.3: n_b ≥ 4).
DEFAULT_MIN_BLOCKS: int = 4

# Number of largest-B points used in plateau check (§3.3: "3 to 4 points").
DEFAULT_PLATEAU_WINDOW: int = 4

# Cross-validation tolerance between SEM_block and SEM_auto (§3.6: 30%).
# Common sources of disagreement: (1) ACF truncation error, (2) plateau
# identification ambiguity, (3) finite-sample fluctuation. 30% avoids
# excessive false-positive warnings for typical 5k–50k frame trajectories.
DEFAULT_CROSS_VALID_RTOL: float = 0.3

# ---------------------------------------------------------------------------
# Step 1: Running average drift
# ---------------------------------------------------------------------------

# Drift factor (§1.3: D < factor × SEM).
DEFAULT_DRIFT_FACTOR: float = 3.0

# ---------------------------------------------------------------------------
# Step 4: Geweke stationarity
# ---------------------------------------------------------------------------

# Front/rear segment fractions (§4.2: f_A = 0.1, f_B = 0.5).
DEFAULT_GEWEKE_FA: float = 0.1
DEFAULT_GEWEKE_FB: float = 0.5

# Critical z-value, two-sided 95% (§4.3: |z| < 1.96).
DEFAULT_GEWEKE_ZCRIT: float = 1.96

# Minimum N_eff in front segment for reliable spectral variance estimation.
# If N_eff_A < this, GewekeResult.reliable is set to False.
DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES: int = 10

# ---------------------------------------------------------------------------
# Global precision
# ---------------------------------------------------------------------------

# Default free-energy tolerance (0.05 eV ≈ 1 kcal/mol chemical accuracy).
DEFAULT_EPSILON_TOL_EV: float = 0.05

# Unit conversion: eV → Hartree (for internal consistency with CP2K a.u.).
EV_TO_HARTREE: float = 1.0 / 27.211386245988

# ---------------------------------------------------------------------------
# Output filenames
# ---------------------------------------------------------------------------

DEFAULT_REPORT_CSV_NAME: str = "ti_convergence_report.csv"
DEFAULT_DIAGNOSTICS_PNG_PREFIX: str = "ti_diag"
DEFAULT_FE_PROFILE_PNG_NAME: str = "ti_free_energy.png"
DEFAULT_FE_CSV_NAME: str = "ti_free_energy.csv"
DEFAULT_SUMMARY_TXT_NAME: str = "ti_summary.txt"

# Standalone single-point output filenames
DEFAULT_STANDALONE_PNG_NAME: str = "ti_single_point_diag.png"
DEFAULT_STANDALONE_CSV_NAME: str = "ti_single_point.csv"
```

---

## 4. Analysis Engines (`analysis/`)

### 4.0 Shared core: `_acf_core.py`

Both `autocorrelation.py` (full series) and `geweke.py` (sub-series variance estimation) need autocorrelation. This module provides the shared implementation.

**τ_corr definition** (Sokal 1997 convention): `τ_int = 1/2 + Σ_{j=1}^{M} C(j)`, where C(j) is the normalized ACF. This is the integrated autocorrelation time. With this definition, `N_eff = N / (2τ_int)` is exact. The `compute_iat` function returns τ_int per this convention (including the 1/2 term).

```
_acf_core.py
├── compute_acf(series) → np.ndarray
│       FFT-accelerated normalized ACF. O(N log N) instead of O(N²).
│       Returns C(j) for j = 0, ..., N-1 with C(0) = 1.
│       Raises ValueError if series has zero variance (constant input).
│
├── compute_iat(acf, alpha, n) → float
│       Self-consistent truncation IAT (§2.2).
│       Iteratively sums C(j) until j ≥ α · τ_corr.
│       Hard upper limit: M ≤ N // 2 (prevents noise-dominated estimates).
│       Returns τ_int (in frames), including the 1/2 base term.
│
└── compute_sem_corrected(series, tau_corr) → float
        σ_λ √(2τ/N). Used by both autocorrelation and Geweke engines.
```

**Design decision**: FFT-based ACF is essential. For a 10,000-frame series, naive O(N²) ≈ 10⁸ ops; FFT ≈ 10⁵ ops. Since AIMD trajectories are typically 5k–50k frames, the difference is noticeable in batch analysis across K points.

**Caution on short series**: When τ_corr > N/10, the IAT estimate itself becomes unreliable (ACF noise at large lags). The M ≤ N//2 hard cap mitigates this but cannot fully eliminate the bias. Downstream consumers (workflow) should flag this condition.

### 4.1 `autocorrelation.py` — Step 2

```
analyze_autocorrelation(
    series:   np.ndarray,          # λ time series, shape (N,)
    *,
    sem_max:  float | None = None, # precision target (None = skip SEM pass/fail)
    alpha:    int   = DEFAULT_ACF_ALPHA,
    neff_min: int   = DEFAULT_NEFF_MIN,
) → AutocorrResult
```

Logic:
1. Call `_acf_core.compute_acf` → C(j)
2. Call `_acf_core.compute_iat` → τ_corr (with M ≤ N//2 hard cap)
3. N_eff = N / (2τ_corr)
4. SEM_auto = σ_λ √(2τ_corr / N)
5. passed_sem = (SEM_auto ≤ sem_max) if sem_max is not None else None
6. If N_eff < neff_min → compute T_min = 2 · neff_min · τ_corr (= 100·τ when neff_min=50)
7. Return `AutocorrResult`

**Note on T_min derivation**: T_min = N_min · Δt, where N_min = 2 · neff_min · τ_corr. The factor 2 comes from N_eff = N/(2τ), requiring N ≥ 2 · neff_min · τ_corr to achieve N_eff ≥ neff_min. This assumes τ_corr remains approximately constant as the simulation is extended.

### 4.2 `block_average.py` — Step 3

```
analyze_block_average(
    series:          np.ndarray,
    *,
    sem_max:         float | None = None,   # precision target (None = skip SEM pass/fail)
    sem_auto:        float | None = None,   # for cross-validation
    min_blocks:      int   = DEFAULT_MIN_BLOCKS,
    plateau_window:  int   = DEFAULT_PLATEAU_WINDOW,
    plateau_rtol:    float = DEFAULT_PLATEAU_RTOL,
    cross_valid_rtol: float = DEFAULT_CROSS_VALID_RTOL,
) → BlockAverageResult
```

Logic:
1. Generate block sizes B = 1, 2, 4, 8, ..., up to N // min_blocks
2. For each B: compute n_b block means, then SEM(B) per §3.2
3. Record `sem_at_max_B` = SEM at the largest B tested (for fallback transparency)
4. Plateau check: take last `plateau_window` entries where n_b ≥ min_blocks, compute δ_s / s̄
5. If `sem_auto` provided: cross-validate |SEM_block − SEM_auto| / SEM_block
6. passed = (plateau_reached AND SEM_plateau ≤ sem_max) if sem_max is not None else None
7. Return `BlockAverageResult`

**Fallback transparency**: When `sem_final` falls back to `SEM_auto` (plateau not reached), `sem_at_max_B` is reported alongside to help the user judge whether SEM(B) was still growing — a strong signal that more simulation time is needed.

### 4.3 `running_average.py` — Step 1

```
analyze_running_average(
    series: np.ndarray,
    *,
    sem:          float,                    # from Step 2 or 3
    drift_factor: float = DEFAULT_DRIFT_FACTOR,
) → RunningAverageResult
```

Logic:
1. Compute λ̄(n) for n = 1, ..., N (cumulative mean)
2. D = max − min over [N/2, N]
3. Compare D vs drift_factor × SEM
4. Return `RunningAverageResult`

**Note**: This module takes SEM as input (not computing it internally), because of the dependency inversion described in `convergence_criteria.md` §5 and Q4. The workflow orchestrator is responsible for passing the correct value.

### 4.4 `geweke.py` — Step 4

```
analyze_geweke(
    series: np.ndarray,
    *,
    f_a:    float = DEFAULT_GEWEKE_FA,
    f_b:    float = DEFAULT_GEWEKE_FB,
    z_crit: float = DEFAULT_GEWEKE_ZCRIT,
    alpha:  int   = DEFAULT_ACF_ALPHA,
    min_neff_subseries: int = DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES,
) → GewekeResult
```

Logic:
1. Slice front A = series[:N_A] and rear B = series[N-N_B:]
2. For each sub-series: call `_acf_core` to get τ_corr, then compute Ŝ² = 2τ · σ² / n
3. Compute N_eff_A = n_A / (2τ_A). If N_eff_A < min_neff_subseries, set `reliable = False`
4. z = (λ̄_A − λ̄_B) / √(Ŝ²_A + Ŝ²_B)
5. passed = |z| < z_crit (computed regardless of reliability)
6. Return `GewekeResult`

**Reliability flag**: For N=5000 with f_A=0.1, the front segment has only 500 frames. If τ_corr is large (e.g., τ=50), N_eff_A ≈ 5, making the spectral variance estimate unreliable. When `reliable=False`, the Geweke result should be interpreted with caution. The workflow includes this in its diagnostic output but does not override the pass/fail from other engines.

---

## 5. Integration (`integration.py`)

Handles the multi-point layer: quadrature weights, per-point precision targets, error propagation, and optional optimal resource allocation.

**Public** (independently useful):
```
compute_trapezoid_weights(xi: np.ndarray) → np.ndarray
    Non-uniform trapezoid weights (§1.4.1).
    w_1 = h_1/2, w_k = (h_{k-1} + h_k)/2, w_K = h_{K-1}/2.
    Validates: sum(w) ≈ xi[-1] − xi[0].
    Raises ValueError if K < 2 (trapezoid rule requires at least 2 points).
    Raises ValueError if xi is not strictly monotonic (non-sorted input rejected).
```

**Private** (called only by `workflow.analyze_ti`):
```
_compute_sem_targets(
    weights: np.ndarray,
    epsilon_tol_au: float,        # ε_tol in Hartree (same unit as λ)
) → np.ndarray
    Per-point SEM upper bound (§1.4.3 Scheme A).
    SEM_{max,k} = ε_tol_au / (|w_k| √K).
    Safety floor: SEM_{max,k} is capped at max(SEM_{max,k}, C · ε_tol_au)
    to prevent degenerate tolerance when w_k → 0 on non-uniform grids.

_integrate_free_energy(
    forces: np.ndarray,           # λ̄ at each point (dA/dξ estimates)
    weights: np.ndarray,
    sems:   np.ndarray,           # SEM_final at each point
) → tuple[float, float]
    Returns (ΔA, σ_A) where:
    ΔA = Σ w_k F_k
    σ_A = √(Σ w_k² SEM_k²)
    Correct for independent constraint points (no inter-point correlation).

_suggest_time_allocation(
    weights: np.ndarray,
    sigmas:  np.ndarray,          # σ_λ per point (from trial run)
    tau_corrs: np.ndarray,        # τ_corr per point
) → np.ndarray
    Optimal time ratios ∝ |w_k| σ_k √τ_k (§1.4.3 Scheme C).
    Returns normalized array that sums to 1.
    Derivation: Lagrange-multiplier minimization of σ_A² subject to ΣT_k = T_total.
```

---

## 6. Orchestrator (`workflow.py`)

The **only** module that knows the execution order of Steps 1–4 and the inter-step data dependencies. No analysis engine imports another engine; all data flows through the orchestrator.

Additionally, `workflow.py` hosts the unified entry points (`standalone_diagnostics`, CSV export functions) that combine analysis + plotting + I/O — following the pattern of `surface_charge_analysis` in `SurfaceCharge.py`. This keeps `plot.py` as a pure rendering module.

### 6.0 Pre-analysis validation (all paths)

Before calling any engine, the orchestrator performs:

1. **NaN handling**: Count NaN values in `lambda_series`. If NaN fraction > `DEFAULT_NAN_FRACTION_MAX` (10%), raise `InsufficientSamplingError`. Otherwise, truncate to the last contiguous non-NaN segment and emit a warning with the number of discarded frames.
2. **Length check**: If remaining N < `DEFAULT_N_MIN` (10), raise `InsufficientSamplingError`.
3. **Short-series warnings**: If N < `DEFAULT_N_WARN_UNRELIABLE` (100), append "extremely short series — results unreliable" to `failure_reasons`. If N < `DEFAULT_N_WARN_SHORT` (500), log a cautionary warning.
4. **Equilibration trimming**: Discard first `equilibration` frames. This is the **sole location** where trimming occurs — `io.py` always returns full (untrimmed) series.

### 6.0.1 Engine override dispatch

`**engine_overrides` are distributed to engines by parameter name matching. Shared parameters (e.g., `alpha`) are forwarded to all engines that accept them:

| Override key | Forwarded to |
|---|---|
| `alpha` | `analyze_autocorrelation`, `analyze_geweke` |
| `neff_min` | `analyze_autocorrelation` |
| `drift_factor` | `analyze_running_average` |
| `f_a`, `f_b`, `z_crit` | `analyze_geweke` |
| `min_blocks`, `plateau_window`, `plateau_rtol`, `cross_valid_rtol` | `analyze_block_average` |

### 6.0.2 Unit conversion for epsilon_tol

`analyze_ti` accepts `epsilon_tol_ev` (user-facing, eV) and converts internally:
```python
epsilon_tol_au = epsilon_tol_ev * EV_TO_HARTREE
```
All downstream SEM targets and comparisons use Hartree (a.u.), matching the unit of `lambda_series` from CP2K.

### 6.1 Single-point analysis (TI context)

```
analyze_single_point(inp: ConstraintPointInput, **overrides) → ConstraintPointReport
```

Execution order (matches `convergence_criteria.md` §5 but respects the data dependency):

```
┌─────────────────────────────────────────────────────────┐
│ 1. autocorrelation.analyze_autocorrelation              │
│    → AutocorrResult (τ_corr, N_eff, SEM_auto)           │
│                                                         │
│ 2. block_average.analyze_block_average                  │
│    → BlockAverageResult (SEM_block, plateau check)      │
│    cross-validates against SEM_auto from step 1         │
│                                                         │
│ 3. running_average.analyze_running_average               │
│    → RunningAverageResult (drift D)                     │
│    receives SEM = SEM_block from step 2 as input        │
│                                                         │
│ 4. geweke.analyze_geweke                                │
│    → GewekeResult (z statistic)                         │
│                                                         │
│ 5. Aggregate: all four passed? Collect failure reasons.  │
└─────────────────────────────────────────────────────────┘
```

**Note on SEM selection (3-tier hierarchy)**: `sem_final` is selected as:
1. `SEM_arctan` (arctan asymptote) when `block_avg.arctan` is not None and reliable — the primary estimator.
2. `SEM_block` (plateau mean) when arctan is unavailable/unreliable but plateau is reached.
3. `SEM_auto` (ACF-based) as the final fallback when neither arctan nor plateau is usable.
Each fallback appends a diagnostic message to `failure_reasons`.

### 6.2 Standalone single-point analysis (no TI context)

**Motivation**: During production runs, individual constraint points finish at different times. The user needs to diagnose a single point's sampling quality immediately — without defining the full set of K points, without knowing trapezoid weights, and without specifying a global ε_tol.

```
analyze_standalone(
    lambda_series: np.ndarray,
    *,
    dt:          float = 1.0,              # frame interval (fs)
    xi:          float = 0.0,              # CV label (for reporting/plot titles)
    sem_target:  float | None = None,      # optional user-specified precision target
    equilibration: int = 0,                # frames to discard from start
    **engine_overrides,
) → ConstraintPointReport
```

This is a **thin wrapper** around `analyze_single_point` that constructs a `ConstraintPointInput` with `weight=None` and `sem_max=sem_target` (which may be None). The internal logic is identical — all four engines run — but pass/fail semantics differ:

| Condition | With `sem_target` | Without `sem_target` |
|-----------|-------------------|---------------------|
| N_eff ≥ 50 | pass/fail as normal | pass/fail as normal |
| SEM ≤ SEM_max | pass/fail against `sem_target` | **skipped** — reported as `None` |
| Drift D < 3×SEM | pass/fail as normal | pass/fail as normal |
| Geweke \|z\| < 1.96 | pass/fail as normal | pass/fail as normal |

When `sem_target` is None, `ConstraintPointReport.passed` is also `None` (rather than True/False) to explicitly signal "diagnostics ran but precision adequacy is undefined". The four engine results are always fully populated.

**Use cases**:
- Quick sanity check while simulations are still running
- Exploring how many more frames are needed (inspect N_eff and τ_corr)
- Teaching/demonstration without setting up a full TI project

### 6.3 Full TI analysis

```
analyze_ti(
    xi_values:      np.ndarray,              # shape (K,)
    lambda_series:  list[np.ndarray],        # K arrays, each shape (N_k,)
    dt:             float,
    *,
    epsilon_tol_ev: float = DEFAULT_EPSILON_TOL_EV,  # user-facing unit
    equilibration:  int | list[int] = 0,     # frames to discard per point
    **engine_overrides,
) → TIReport
```

Logic:
1. Convert `epsilon_tol_ev` → `epsilon_tol_au` (Hartree)
2. Validate `xi_values`: must be strictly monotonic, K ≥ 2
3. Compute trapezoid weights and SEM targets via `integration._compute_sem_targets`
4. Apply equilibration trimming + NaN handling (§6.0) to each λ series
5. Construct `ConstraintPointInput` with 0-based `point_index`
6. Call `analyze_single_point` for each of the K points
7. Call `_integrate_free_energy` with the per-point λ̄ and SEM_final
8. Optionally call `_suggest_time_allocation` if any point failed
9. Assemble `TIReport`

---

## 7. I/O (`io.py`)

Pure I/O layer: bridges between filesystem conventions and typed data structures. **Does not** compute integration weights or SEM targets — that is the orchestrator's responsibility (§6.3).

### 7.1 Expected directory layout for constrained TI

**Primary pattern** — subdirectory mode (matches TIGen output, see `scripts/TIGen.py`):

```
ti_project/
├── ti_target_0.031369/         # generated by TIGen (format: ti_target_{cv:.6f})
│   ├── cMD-1.restart
│   └── cMD-constraint_force.dat-1.LagrangeMultLog
├── ti_target_-0.239239/
│   ├── cMD-1.restart
│   └── cMD-constraint_force.dat-1.LagrangeMultLog
├── ...
└── ti_target_1.114938/
    ├── cMD-1.restart
    └── cMD-constraint_force.dat-1.LagrangeMultLog
```

**Alternative pattern** — user-defined `xi_*` subdirectories:

```
ti_project/
├── xi_1.0/
│   ├── *.restart
│   └── *.LagrangeMultLog
└── xi_5.5/
    ├── *.restart
    └── *.LagrangeMultLog
```

**Reference data**: `data_example/ti/` contains 8 real constraint points in the `ti_target_*` format, suitable for integration tests.

### 7.2 Public functions

```
discover_ti_points(
    root_dir: Path,
    *,
    pattern: str = "auto",        # "ti_target" | "xi" | "auto"
) → list[TIPointDefinition]
    Returns list of TIPointDefinition sorted by xi.
    Auto-detection: tries ti_target_* first (TIGen convention), then xi_*.
    Extracts xi by parsing the numeric suffix from the directory name.
    Raises FileNotFoundError if no matching directories found.

load_ti_series(
    point_defs: list[TIPointDefinition],
) → list[tuple[float, np.ndarray, float]]
    Calls ColvarParser for each point.
    Returns list of (xi, lambda_series, dt) tuples with FULL (untrimmed) series.
    No equilibration trimming here — that is workflow's responsibility (§6.0).
```

**Reuse**: internally calls `ColvarParser.parse_colvar_restart` and `ColvarParser.parse_lagrange_mult_log` — the same functions used by slow growth. No new parser code needed.

**Design rationale for split**: The original `parse_ti_points` combined I/O + weight computation + trimming. This violated separation of concerns — `io.py` had to import `integration.py`. Now `io.py` is pure I/O (file discovery + parsing), and the orchestrator (`workflow.analyze_ti`) handles weight computation, trimming, and `ConstraintPointInput` construction.

---

## 8. Plotting (`plot.py`) and Output (`workflow.py`)

`plot.py` is a **pure rendering module** — it accepts dataclasses, produces PNGs. It never calls analysis functions or ColvarParser. Follows `SlowGrowthPlot.py` conventions: lazy `matplotlib` import, `Agg` backend, `dpi=180`, returns `Path`.

CSV exports and the unified `standalone_diagnostics` entry point live in `workflow.py` (the orchestrator), since they coordinate analysis + plotting + I/O.

### 8.1 Per-point diagnostics (one PNG per point) — `plot.py`

A 2×2 subplot figure containing:

| Panel | Content | Source |
|-------|---------|--------|
| Top-left | λ̄(n) cumulative mean + horizontal band at λ̄ ± SEM | `RunningAverageResult.running_mean` |
| Top-right | C(j) autocorrelation decay + vertical line at M = α·τ | `AutocorrResult.acf` |
| Bottom-left | SEM(B) vs B + horizontal line at SEM_max + plateau region shaded | `BlockAverageResult.sem_curve` |
| Bottom-right | Text summary table: τ_corr, N_eff, SEM_auto, SEM_block, z, pass/fail | All results |

```
plot_point_diagnostics(
    report: ConstraintPointReport,
    *,
    output_dir: Path | None = None,
) → Path
    Saves: {output_dir}/{DEFAULT_DIAGNOSTICS_PNG_PREFIX}_xi{xi:.4f}.png
```

**Standalone adaptation**: When `report.sem_max is None`, the bottom-left panel omits the SEM_max horizontal line and the bottom-right panel shows "SEM target: N/A" instead of a pass/fail verdict for precision. All other panels render identically.

**Fallback SEM display**: When block average plateau was not reached, the bottom-right panel shows both `SEM_auto` (used as `sem_final`) and `sem_at_max_B` to help the user judge whether SEM(B) was still growing.

### 8.2 Standalone single-point unified entry point — `workflow.py`

Convenience function that combines parsing + `analyze_standalone` + `plot_point_diagnostics` + CSV in a single call. This is the primary entry point for the standalone use case, both from the Python API and from CLI command 305.

```
standalone_diagnostics(
    restart_path:   str,
    log_path:       str,
    *,
    equilibration:  int = 0,
    sem_target:     float | None = None,      # optional precision target (same unit as λ)
    colvar_id:      int | None = None,
    output_dir:     Path | None = None,
) → dict[str, Path | ConstraintPointReport]
    Returns {"report": ConstraintPointReport, "diagnostics_png": Path, "csv": Path}.
```

Internally:
1. Parse restart + LagrangeMultLog via `ColvarParser` (reuse existing)
2. Call `analyze_standalone` → `ConstraintPointReport`
3. Call `plot.plot_point_diagnostics` → diagnostic PNG
4. Call `write_single_point_csv` → CSV with time series statistics
5. Print console summary (τ_corr, N_eff, SEM, z, pass/fail)

**Design rationale**: Placing this in `workflow.py` (not `plot.py`) ensures the dependency is `workflow → plot` (one-way), and CLI commands 303/304/305 all import from `workflow` consistently.

### 8.3 Free-energy profile — `plot.py`

Single figure with:
- dA/dξ (force) vs ξ with error bars (left axis)
- Integrated A(ξ) vs ξ with shaded σ_A band (right axis)
- Color-coded x-ticks: green = passed, red = failed

```
plot_free_energy_profile(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) → Path
    Saves: {output_dir}/{DEFAULT_FE_PROFILE_PNG_NAME}
```

### 8.4 CSV exports — `workflow.py`

```
write_convergence_csv(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) → Path
    One row per constraint point. Energy quantities in eV. Columns:
    xi, lambda_mean_eV, sigma_lambda_eV, tau_corr, n_eff,
    sem_auto_eV, sem_block_eV, delta_sem_block_eV, plateau_B, plateau_reached,
    sem_final_eV, sem_final_method, sem_max_eV, geweke_z, geweke_reliable,
    drift_D_eV, passed, failure_reasons

write_free_energy_csv(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) → Path
    Columns: xi, weight, dA_dxi_eV, sem_eV, A_integrated_eV, sigma_A_cumulative_eV

write_single_point_csv(
    report: ConstraintPointReport,
    *,
    output_dir: Path | None = None,
) → Path
    Single-row CSV with all diagnostic quantities for one point.
    Used by standalone_diagnostics (§8.2) and optionally by the
    per-point loop in the multi-point pipeline.
    Columns: xi, lambda_mean_eV, sigma_lambda_eV, tau_corr, n_eff,
    sem_auto_eV, sem_block_eV, delta_sem_block_eV, plateau_B, plateau_reached,
    sem_final_eV, sem_final_method, sem_max_eV, geweke_z, geweke_reliable,
    drift_D_eV, passed, failure_reasons
```

---

## 9. CLI Integration

### 9.1 New menu entries

In `cli/_enhanced_sampling.py`, add alongside existing 301/302:

| Code | Label | Function |
|------|-------|----------|
| 303 | Constrained TI — convergence check | Run 4-step diagnostics on all points, print summary, save report |
| 304 | Constrained TI — free-energy profile | Full pipeline: diagnose + integrate + plot |
| 305 | Constrained TI — single-point diagnostics | Diagnose one constraint point independently |

### 9.2 Command class sketch

```python
class TIConvergenceCmd(MenuCommand):
    """303: Run convergence diagnostics on constrained TI data."""

    def _collect_all_params(self) -> dict:
        # Prompt for root directory containing constraint-point subdirectories
        # Prompt for ε_tol (default: 1.0 kcal/mol)
        # Prompt for equilibration frames to discard
        # Optional: advanced params (α, N_eff_min, etc.)
        ...

    def execute(self, ctx: dict) -> None:
        analyze_ti = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "analyze_ti",
        )
        ...

class TIFreeEnergyCmd(MenuCommand):
    """304: Full constrained TI pipeline with plots."""
    # Extends 303 by also calling plot_point_diagnostics + plot_free_energy_profile
    ...

class TISinglePointCmd(MenuCommand):
    """305: Standalone single-point convergence diagnostics."""

    def _collect_all_params(self) -> dict:
        # Reuse _discover_restart_file / _discover_log_file from slow-growth commands
        # Prompt for equilibration frames (default: 0)
        # Prompt for SEM target (optional — user can press Enter to skip)
        # Prompt for output directory
        ...

    def execute(self, ctx: dict) -> None:
        standalone_diagnostics = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "standalone_diagnostics",
        )
        standalone_diagnostics(
            ctx[K.RESTART_PATH],
            ctx[K.LOG_PATH],
            equilibration=ctx.get(K.EQUILIBRATION, 0),
            sem_target=ctx.get(K.SEM_TARGET),
            colvar_id=ctx.get(K.COLVAR_ID),
            output_dir=Path(ctx[K.OUTDIR]),
        )
```

**Note**: `TISinglePointCmd` reuses the same file-discovery helpers (`_discover_restart_file`, `_discover_log_file`) and metadata display (`_print_sg_info`) already defined for the slow-growth commands. No duplication.

### 9.3 Menu tree change

```python
# In cli/__init__.py, build_menu_tree():
es_group.add(
    SGQuickPlotCmd("301", "Slow Growth — Quick Plot"),
    SGPublicationPlotCmd("302", "Slow Growth — Publication Plot"),
    "─" * 40,   # separator
    TIConvergenceCmd("303", "Constrained TI — Convergence Check"),
    TIFreeEnergyCmd("304", "Constrained TI — Free Energy Profile"),
    TISinglePointCmd("305", "Constrained TI — Single-Point Diagnostics"),
)
```

---

## 10. Dependency Graph

```
         ┌──────────────────────────────────────────┐
         │       cli (303/304/305)                   │
         └──────┬──────────────────┬────────────────┘
                │                  │
   ┌────────────▼──────┐   ┌──────▼──────────────────┐
   │ main.py           │   │ 305: standalone path     │
   │ run_constrained_ti│   │ (skips integration.py)   │
   └────────┬──────────┘   └──────┬──────────────────┘
            │                     │
            ▼                     ▼
  ┌───────────────────────────────────────────────────┐
  │          workflow.py (orchestrator + entry points) │
  │                                                   │
  │  analyze_ti()              ← 303/304              │
  │    calls analyze_single_point per K pts           │
  │    calls integration.py for ΔA, σ_A               │
  │                                                   │
  │  analyze_standalone()      ← 305                  │
  │    calls analyze_single_point once                │
  │    NO integration.py dependency                   │
  │                                                   │
  │  standalone_diagnostics()  ← 305 unified entry    │
  │    calls analyze_standalone + plot + CSV           │
  │                                                   │
  │  write_*_csv()             ← CSV export helpers   │
  └──┬─────┬─────┬─────┬──────────┬──────────────────┘
     │     │     │     │          │
     ▼     ▼     ▼     ▼          ▼
  ┌──────┐┌──────┐┌──────┐┌──────────────┐  ┌──────────────┐
  │run_  ││auto  ││geweke││block_average │  │  plot.py     │
  │avg   ││corr  ││      ││              │  │  (pure render)│
  └──────┘└──┬───┘└──┬───┘└──────────────┘  └──────────────┘
             │       │
        ┌────▼───────▼──┐
        │  _acf_core.py │
        └───────────────┘

  ┌──────────────────────────┐
  │   integration.py         │   ← used by workflow.analyze_ti only
  │   (weights, ΔA, σ_A)    │      NOT by analyze_standalone
  └──────────────────────────┘

  ┌──────────────────────────┐
  │   io.py (pure I/O)       │   ← file discovery + ColvarParser
  │   discover_ti_points     │      NO dependency on integration.py
  │   load_ti_series         │      NO equilibration trimming
  └──────────┬───────────────┘
             │ reuses
  ┌──────────▼───────────────┐
  │ utils/RestartParser/     │
  │   ColvarParser (existing)│
  └──────────────────────────┘

  ┌──────────────────────────┐
  │   models.py              │
  │   (shared by all above)  │
  └──────────────────────────┘
```

Key constraints:
- No analysis engine imports another analysis engine. Data flows exclusively through the orchestrator.
- `workflow.py` may import `plot.py`; `plot.py` must **never** import `workflow.py` (one-way dependency).
- `io.py` must **never** import `integration.py` — weight computation belongs in the orchestrator.
- `analyze_standalone` has **no dependency** on `integration.py` — it can run without knowing K, weights, or ε_tol.
- `plot_point_diagnostics` is shared between the multi-point and standalone paths; it adapts rendering based on whether `sem_max` is None.
- Equilibration trimming occurs **only** in `workflow.py` (§6.0) — `io.py` always returns full series.

---

## 11. Error Handling

### 11.1 Exception hierarchy

Defined in `models.py` (see §2.5):

```python
class ConvergenceError(MDAnalysisError):
    """Base for constrained-TI convergence errors."""

class InsufficientSamplingError(ConvergenceError):
    """Raised when sampling is too short for any reliable estimate."""
```

### 11.2 Error vs. warning policy

| Condition | Severity | Behavior |
|-----------|----------|----------|
| K < 2 in `analyze_ti` | **Error** | Raise `ValueError` (trapezoid rule needs ≥ 2 points) |
| xi not strictly monotonic | **Error** | Raise `ValueError` in `compute_trapezoid_weights` |
| N < 10 (after NaN removal + equilibration trim) | **Error** | Raise `InsufficientSamplingError` |
| NaN fraction > 10% in lambda_series | **Error** | Raise `InsufficientSamplingError` |
| equilibration ≥ N (nothing left after trim) | **Error** | Raise `InsufficientSamplingError` |
| Series has zero variance (constant λ) | **Error** | Raise `ValueError` in `compute_acf` |
| NaN values present (≤ 10%) | **Warning** | Truncate to last contiguous non-NaN segment, log count |
| N < 100 (after trim) | **Warning** | Append "extremely short series — results unreliable" to `failure_reasons` |
| N < 500 (after trim) | **Warning** | Log cautionary message |
| N_eff < 50 | **Warning** | Continue analysis, flag in report, compute T_min |
| τ_corr > N/10 | **Warning** | Flag "IAT estimate may be unreliable" in report |
| Plateau not reached | **Warning** | Fall back to SEM_auto, report both SEM_auto and sem_at_max_B |
| Geweke \|z\| ≥ 1.96 | **Warning** | Flag in report, suggest increasing equilibration |
| Geweke sub-series N_eff_A < 10 | **Warning** | Set `GewekeResult.reliable = False`, note in report |
| SEM > SEM_max | **Warning** | Flag in report, suggest extending simulation |
| ACF and block SEM disagree > 30% | **Warning** | Flag in report (possible truncation issue) |

The orchestrator never short-circuits on warnings. All four steps always run to completion so the report is maximally informative. Only truly unrecoverable conditions (degenerate input) raise exceptions.

---

## 12. Testing Strategy

### 12.1 Synthetic data generators (test fixtures)

```python
def make_ar1(n: int, phi: float, sigma: float, rng) -> np.ndarray:
    """AR(1) process with known τ_corr = (1+φ)/(1−φ) / 2."""

def make_drifting(n: int, slope: float, sigma: float, rng) -> np.ndarray:
    """Linear drift + white noise — should fail Geweke."""

def make_iid(n: int, mu: float, sigma: float, rng) -> np.ndarray:
    """i.i.d. Gaussian — τ_corr = 0.5, N_eff = N."""

def make_step_change(n: int, mu1: float, mu2: float, sigma: float, rng) -> np.ndarray:
    """Two-segment stationary series with different means — should fail Geweke."""

def make_with_nans(n: int, nan_positions: list[int], rng) -> np.ndarray:
    """Inject NaN at specified positions into otherwise clean AR(1)."""

def make_short_ar1(n: int = 15, phi: float = 0.8, sigma: float = 1.0, rng = ...) -> np.ndarray:
    """Short correlated series for boundary testing (N_warn_unreliable zone)."""
```

All generators accept `rng = np.random.default_rng(seed)` for reproducibility.

### 12.2 Unit test targets

| Module | Key assertions |
|--------|---------------|
| `_acf_core` | AR(1) N≥5000 → τ_corr within 5% of analytic; N=500 → within 20%; i.i.d. → abs(τ_corr - 0.5) < 0.1; constant series → raises ValueError; high-phi (0.95) → self-consistent iteration converges |
| `autocorrelation` | N_eff = N for i.i.d. (rtol=1e-6); SEM_auto = σ/√N for i.i.d. (rtol=1e-6); sem_max=None → passed_sem is None; T_min = 2·neff_min·τ |
| `block_average` | Plateau at correct SEM for AR(1); SEM(B=1) = naive SEM for i.i.d. (exact); sem_max=None → passed is None; sem_at_max_B always populated |
| `running_average` | Converged i.i.d. → D < 3×SEM; drifting → D ≥ 3×SEM; field renamed to drift_limit |
| `geweke` | Stationary → \|z\| < 1.96; drifting → \|z\| > 1.96; step_change → \|z\| > 1.96; short front segment → reliable=False |
| `integration` | K=2 exact formula; uniform-spacing weights correct; K=1 → ValueError; non-monotonic xi → ValueError; SEM_max safety floor enforced |
| `workflow` — `analyze_standalone` | sem_max=None → passed is None; sem_target set → passed is True/False; NaN input handled; equilibration ≥ N → error; failure_reasons non-empty when passed=False, empty when passed=True |
| `workflow` — `analyze_ti` | All-pass on clean AR(1); at least one fail on drifting; unequal-length series (N_1≠N_2) handled; epsilon_tol_ev correctly converted to Hartree |
| `models` | Frozen dataclass immutability; TIPointDefinition fields typed |

### 12.3 Integration tests

Use **`data_example/ti/`** (8 real constraint points in `ti_target_*` format). End-to-end:

1. `discover_ti_points(data_example/ti/)` → verify 8 `TIPointDefinition` entries, sorted by xi
2. `load_ti_series(point_defs)` → verify all lambda_series are non-empty numpy arrays
3. `analyze_ti(xi_values, lambda_series, dt)` → verify `TIReport` has K=8 points, no crashes
4. `standalone_diagnostics(restart_path, log_path)` → verify output files exist (PNG + CSV)

For early milestones (Phase 8a/10a), existing `data_example/sg/` LagrangeMultLog can serve as single-point test data (same file format, different physics).

---

## 13. Implementation Order

Recommended build sequence, each step testable before proceeding:

| Phase | Files | Deliverable |
|-------|-------|-------------|
| 0 | Commit `DESIGN.md` to version control | Design frozen as reference |
| 1 | `config.py`, `models.py` (incl. exceptions + `TIPointDefinition`) | Frozen dataclasses importable; all defaults defined; **field names frozen before Phase 2** |
| 2 | `analysis/_acf_core.py` + tests | FFT-ACF and IAT verified against AR(1) analytics; hard cap M ≤ N//2 tested |
| 3 | `analysis/autocorrelation.py` + tests | Step 2 engine (depends on Phase 2 interface) |
| 4 | `analysis/block_average.py` + tests | Step 3 engine (independent of Phase 2) |
| 5 | `analysis/running_average.py` + tests | Step 1 engine (independent of Phase 2) |
| 6 | `analysis/geweke.py` + tests | Step 4 engine (depends on Phase 2 interface) |
| 7 | `integration.py` + tests | Weights, targets, quadrature; K<2 and non-monotonic xi rejection |
| 8a | `workflow.py` — `analyze_standalone` + NaN handling + tests | Single-point orchestration end-to-end |
| 10a | `plot.py` — `plot_point_diagnostics` | Single-point diagnostic PNG (visual QA) |
| 8a+ | `workflow.py` — `standalone_diagnostics` + `write_single_point_csv` | Unified entry point for standalone path |
| 11a | `cli/_enhanced_sampling.py` — CLI 305 only | **★ Milestone B: CLI 305 end-to-end usable** |
| 8b | `workflow.py` — `analyze_single_point` + `analyze_ti` + tests | Multi-point orchestration wired |
| 9 | `io.py` + tests | File discovery (`ti_target_*` pattern) + `load_ti_series` |
| 10b | `plot.py` — `plot_free_energy_profile` | Multi-point profile PNG |
| 8b+ | `workflow.py` — `write_convergence_csv` + `write_free_energy_csv` | Multi-point CSV exports |
| 11b | `cli/_enhanced_sampling.py` — CLI 303/304 | Multi-point CLI commands |
| 12 | `main.py` | `run_constrained_ti_analysis` entry point |
| 13 | `context4agent/` documentation sync | Update architecture, interface_exposure, data_contract, glossary_units, short_term |
| 14 | Integration tests | End-to-end with `data_example/ti/` real data |

### 13.1 Parallelism constraints

**Phase 2 is a prerequisite for Phases 3 and 6** (they import `_acf_core`). Phases 4 and 5 are fully independent of Phase 2. So the actual parallelism is:

```
Phase 1 (must complete first — field names frozen)
    │
    ├── Phase 2 (_acf_core)
    │       │
    │       ├── Phase 3 (autocorrelation)  ─┐
    │       └── Phase 6 (geweke)           ─┤  these 4 can be parallel
    ├── Phase 4 (block_average)            ─┤  after Phase 2 completes
    └── Phase 5 (running_average)          ─┘
```

### 13.2 Milestones

**Milestone A** (Phase 8a): `analyze_standalone` API works with synthetic data — all four engines wired, NaN handling, short-series warnings.

**Milestone B** (Phase 11a): CLI 305 end-to-end — user can run `md-analysis` → 305, point it at a single constraint directory, get a diagnostic PNG + CSV. This is the **minimal viable feature**, usable without multi-point infrastructure.

**Milestone C** (Phase 11b): Full TI pipeline — CLI 303/304 with `data_example/ti/` as reference. Diagnose all K points, integrate ΔA ± σ_A, produce free-energy profile.

### 13.3 Documentation sync (Phase 13)

Per project CLAUDE.md rules, the following `context4agent/` files must be updated:

| Change | File to update |
|--------|---------------|
| New `constrained_ti/` sub-package | `architecture/README.md` + `architecture/modules/src/enhanced_sampling/` |
| New public API | `architecture/modules/src/enhanced_sampling/interface_exposure.md` |
| Implementation patterns | `architecture/modules/src/enhanced_sampling/implementation_guidelines.md` |
| New CSV columns, PNG filenames | `architecture/modules/data_contract.md` |
| kcal/mol → Hartree conversion | `architecture/modules/glossary_units.md` |
| New CLI commands 303/304/305 | `architecture/modules/src/cli/implementation_guidelines.md` |
| Project capabilities | `requirements/short_term.md` |
