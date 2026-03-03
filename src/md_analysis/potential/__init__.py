"""Public interface for `src.potential`."""

from __future__ import annotations

from .CenterPotential import (
    center_slab_potential_analysis,
    fermi_energy_analysis,
    electrode_potential_analysis,
    thickness_sensitivity_analysis,
)
from .PhiZProfile import phi_z_planeavg_analysis
from .config import (
    DEFAULT_CENTER_POTENTIAL_CSV_NAME,
    DEFAULT_CENTER_POTENTIAL_PNG_NAME,
    DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME,
    DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME,
    DEFAULT_FERMI_ENERGY_CSV_NAME,
    DEFAULT_FERMI_ENERGY_PNG_NAME,
    DEFAULT_PHI_Z_PNG_NAME,
    DEFAULT_PHI_Z_STATS_CSV_NAME,
    DEFAULT_SLAB_CENTER_CSV_NAME,
    DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME,
    DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME,
    DEFAULT_THICKNESS_ANG,
)

__all__ = [
    "center_slab_potential_analysis",
    "fermi_energy_analysis",
    "electrode_potential_analysis",
    "thickness_sensitivity_analysis",
    "phi_z_planeavg_analysis",
    "DEFAULT_CENTER_POTENTIAL_CSV_NAME",
    "DEFAULT_CENTER_POTENTIAL_PNG_NAME",
    "DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME",
    "DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME",
    "DEFAULT_FERMI_ENERGY_CSV_NAME",
    "DEFAULT_FERMI_ENERGY_PNG_NAME",
    "DEFAULT_PHI_Z_PNG_NAME",
    "DEFAULT_PHI_Z_STATS_CSV_NAME",
    "DEFAULT_SLAB_CENTER_CSV_NAME",
    "DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME",
    "DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME",
    "DEFAULT_THICKNESS_ANG",
]
