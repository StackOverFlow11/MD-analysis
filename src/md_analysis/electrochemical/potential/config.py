"""Default output filenames for potential analysis."""

from __future__ import annotations

DEFAULT_CENTER_POTENTIAL_CSV_NAME: str = "center_potential.csv"
DEFAULT_CENTER_POTENTIAL_PNG_NAME: str = "center_potential.png"
DEFAULT_FERMI_ENERGY_CSV_NAME: str = "fermi_energy.csv"
DEFAULT_FERMI_ENERGY_PNG_NAME: str = "fermi_energy.png"
DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME: str = "electrode_potential_U_vs_SHE.csv"
DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME: str = "electrode_potential_U_vs_SHE.png"
DEFAULT_SLAB_CENTER_CSV_NAME: str = "slab_center_and_interfaces.csv"
DEFAULT_PHI_Z_STATS_CSV_NAME: str = "phi_z_planeavg_stats.csv"
DEFAULT_PHI_Z_PNG_NAME: str = "phi_z_planeavg_all_frames.png"
DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME: str = "thickness_sensitivity.csv"
DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME: str = "thickness_sensitivity.png"

DEFAULT_THICKNESS_ANG: float = 7.5

# Distributed single-point defaults
DEFAULT_SP_DIR_PATTERN: str = "potential_t*_i*"
DEFAULT_SP_CUBE_FILENAME: str = "sp_potential-v_hartree-1_0.cube"
DEFAULT_SP_OUT_FILENAME: str = "sp.out"
