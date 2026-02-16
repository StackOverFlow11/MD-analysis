"""Public interface for `scripts.structure.Analysis`."""

from __future__ import annotations

from .config import DEFAULT_OUTPUT_DIR
from .config import DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
from .config import DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME
from .config import DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME
from .config import DEFAULT_START_INTERFACE
from .config import DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from .config import DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME
from .config import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from .WaterAnalysis import ad_water_orientation_analysis
from .WaterAnalysis import compute_adsorbed_water_theta_distribution
from .WaterAnalysis import detect_adsorbed_layer_range_from_density_profile
from .WaterAnalysis import water_mass_density_z_distribution_analysis
from .WaterAnalysis import water_orientation_weighted_density_z_distribution_analysis
from .Water import plot_water_three_panel_analysis

__all__ = [
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME",
    "DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME",
    "DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME",
    "DEFAULT_START_INTERFACE",
    "DEFAULT_WATER_MASS_DENSITY_CSV_NAME",
    "DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME",
    "DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME",
    "detect_adsorbed_layer_range_from_density_profile",
    "ad_water_orientation_analysis",
    "compute_adsorbed_water_theta_distribution",
    "water_mass_density_z_distribution_analysis",
    "water_orientation_weighted_density_z_distribution_analysis",
    "plot_water_three_panel_analysis",
]
