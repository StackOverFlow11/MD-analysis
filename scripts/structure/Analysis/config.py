"""Default configuration for structure analysis utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

# Default output directory:
# for future CLI usage, default to the caller's current working directory.
DEFAULT_OUTPUT_DIR: Path = Path.cwd()

# Which metal interface to use as the profile start point along c.
# - "low_c": interface with smaller c fractional coordinate
# - "high_c": interface with larger c fractional coordinate
DEFAULT_START_INTERFACE: Literal["low_c", "high_c"] = "low_c"

# Default output CSV filename for water mass-density analysis.
DEFAULT_WATER_MASS_DENSITY_CSV_NAME: str = "water_mass_density_z_distribution_analysis.csv"

# Default output CSV filename for water orientation-weighted-density analysis.
DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME: str = (
    "water_orientation_weighted_density_z_distribution_analysis.csv"
)

# Default output filenames for adsorbed-water orientation analysis.
DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME: str = "adsorbed_water_orientation_profile.csv"
DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME: str = "adsorbed_water_layer_range.txt"
DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME: str = "adsorbed_water_theta_distribution_0_180.csv"

# Default figure filename for integrated three-panel water analysis plot.
DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME: str = "water_three_panel_analysis.png"
