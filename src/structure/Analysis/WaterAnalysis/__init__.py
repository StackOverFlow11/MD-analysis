"""Public interface for `src.structure.Analysis.WaterAnalysis`."""

from __future__ import annotations

# --- Internal re-exports from _common (used by Water.py, NOT part of public API) ---
from ._common import StartInterface
from ._common import _compute_density_orientation_ensemble

from .AdWaterOrientation import ad_water_orientation_analysis
from .AdWaterOrientation import compute_adsorbed_water_theta_distribution
from .AdWaterOrientation import detect_adsorbed_layer_range_from_density_profile
from .WaterDensity import water_mass_density_z_distribution_analysis
from .WaterOrientation import water_orientation_weighted_density_z_distribution_analysis

__all__ = [
    "detect_adsorbed_layer_range_from_density_profile",
    "ad_water_orientation_analysis",
    "compute_adsorbed_water_theta_distribution",
    "water_mass_density_z_distribution_analysis",
    "water_orientation_weighted_density_z_distribution_analysis",
]
