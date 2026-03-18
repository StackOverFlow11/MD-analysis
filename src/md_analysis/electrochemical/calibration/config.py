"""Constants and defaults for charge-potential calibration."""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Default calibration file location (independent of global config.json)
# ---------------------------------------------------------------------------

DEFAULT_CALIBRATION_DIR: Path = Path.home() / ".config" / "md_analysis"
DEFAULT_CALIBRATION_FILE: Path = DEFAULT_CALIBRATION_DIR / "calibration.json"

# ---------------------------------------------------------------------------
# Output filenames
# ---------------------------------------------------------------------------

DEFAULT_CALIBRATION_CSV_NAME: str = "calibration_data.csv"
DEFAULT_CALIBRATION_PNG_NAME: str = "calibration_fit.png"

# ---------------------------------------------------------------------------
# Nernst equation constants (CODATA 2018) — for RHE conversion
# ---------------------------------------------------------------------------

R_J_PER_MOL_K: float = 8.314462618        # Gas constant (J / mol / K)
F_C_PER_MOL: float = 96485.33212          # Faraday constant (C / mol)
LN10: float = 2.302585092994046           # ln(10)

# ---------------------------------------------------------------------------
# Fitting methods
# ---------------------------------------------------------------------------

FITTING_LINEAR: str = "linear"
FITTING_POLYNOMIAL: str = "polynomial"
FITTING_SPLINE: str = "spline"
DEFAULT_FITTING_METHOD: str = FITTING_LINEAR
DEFAULT_POLY_DEGREE: int = 2

# ---------------------------------------------------------------------------
# Potential reference
# ---------------------------------------------------------------------------

DEFAULT_REFERENCE: str = "SHE"
DEFAULT_TEMPERATURE_K: float = 298.15
DEFAULT_PH: float = 0.0
