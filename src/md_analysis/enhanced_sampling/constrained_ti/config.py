"""Default thresholds for constrained-TI convergence diagnostics.

Every threshold corresponds to a boxed equation in convergence_criteria.md.
All are exposed as keyword arguments in the analysis functions; the values
here are the recommended defaults.

Unit convention for epsilon_tol
-------------------------------
User-facing API accepts kcal/mol (DEFAULT_EPSILON_TOL_KCAL).
``workflow.analyze_ti`` converts to Hartree internally before computing
SEM targets, so that SEM_max has the same unit as lambda_series (a.u.).
The conversion chain: user kcal/mol → KCAL_TO_HARTREE → Hartree (a.u.).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

# Hard minimum: raise InsufficientSamplingError if N < this.
DEFAULT_N_MIN: int = 10

# Tiered warning thresholds for short series.
DEFAULT_N_WARN_UNRELIABLE: int = 100
DEFAULT_N_WARN_SHORT: int = 500

# Maximum fraction of NaN values allowed before raising error.
DEFAULT_NAN_FRACTION_MAX: float = 0.1

# ---------------------------------------------------------------------------
# Step 2: Autocorrelation
# ---------------------------------------------------------------------------

# Self-consistent truncation multiplier (M >= alpha * tau_corr).
# Reference: Sokal (1997), also default in ALPS and emcee.
DEFAULT_ACF_ALPHA: int = 5

# Hard upper limit fraction for ACF truncation window (M <= N * fraction).
DEFAULT_ACF_M_MAX_FRACTION: float = 0.5

# Minimum effective independent samples (N_eff >= 50).
DEFAULT_NEFF_MIN: int = 50

# ---------------------------------------------------------------------------
# Step 3: Block averaging
# ---------------------------------------------------------------------------

# Maximum relative spread for plateau detection (delta_s / s_bar < 0.2).
DEFAULT_PLATEAU_RTOL: float = 0.2

# Minimum number of blocks at largest B (n_b >= 4).
DEFAULT_MIN_BLOCKS: int = 4

# Number of largest-B points used in plateau check.
DEFAULT_PLATEAU_WINDOW: int = 4

# Cross-validation tolerance between SEM_block and SEM_auto (30%).
DEFAULT_CROSS_VALID_RTOL: float = 0.3

# ---------------------------------------------------------------------------
# Step 3b: Arctan extrapolation
# ---------------------------------------------------------------------------

# Minimum R² for arctan fit to be considered reliable.
# Set below 0.95 because dense sampling (B=1..20 continuous) produces many
# precise low-B points where the arctan model is weakest, lowering R²
# without degrading the asymptotic estimate.
DEFAULT_ARCTAN_R2_MIN: float = 0.75

# Minimum number of (non-NaN) block sizes needed to attempt arctan fit.
DEFAULT_ARCTAN_MIN_POINTS: int = 5

# ---------------------------------------------------------------------------
# Step 1: Running average drift
# ---------------------------------------------------------------------------

# Drift factor (D < factor * SEM).
DEFAULT_DRIFT_FACTOR: float = 3.0

# ---------------------------------------------------------------------------
# Step 4: Geweke stationarity
# ---------------------------------------------------------------------------

# Front/rear segment fractions (f_A = 0.1, f_B = 0.5).
DEFAULT_GEWEKE_FA: float = 0.1
DEFAULT_GEWEKE_FB: float = 0.5

# Critical z-value, two-sided 95%.
DEFAULT_GEWEKE_ZCRIT: float = 1.96

# Minimum N_eff in front segment for reliable spectral variance.
DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES: int = 10

# ---------------------------------------------------------------------------
# Global precision
# ---------------------------------------------------------------------------

# Default free-energy tolerance (1 kcal/mol ~ chemical accuracy).
DEFAULT_EPSILON_TOL_KCAL: float = 1.0

# Unit conversion: kcal/mol -> Hartree.
KCAL_TO_HARTREE: float = 1.0 / 627.509474

# Unit conversion: kcal/mol -> eV.
KCAL_TO_EV: float = 0.0433641

# ---------------------------------------------------------------------------
# Output filenames
# ---------------------------------------------------------------------------

DEFAULT_REPORT_CSV_NAME: str = "ti_convergence_report.csv"
DEFAULT_DIAGNOSTICS_PNG_PREFIX: str = "ti_diag"
DEFAULT_FE_PROFILE_PNG_NAME: str = "ti_free_energy.png"
DEFAULT_FE_CSV_NAME: str = "ti_free_energy.csv"
DEFAULT_SUMMARY_TXT_NAME: str = "ti_summary.txt"

# Standalone single-point output filenames.
DEFAULT_STANDALONE_PNG_NAME: str = "ti_single_point_diag.png"
DEFAULT_STANDALONE_CSV_NAME: str = "ti_single_point.csv"
