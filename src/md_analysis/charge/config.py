"""Constants and defaults for charge analysis."""

from __future__ import annotations

# 1 e/Å² = 1.602176634e3 μC/cm²
E_PER_A2_TO_UC_PER_CM2: float = 1.602176634e3

DEFAULT_DIR_PATTERN = "bader_t*_i*"
DEFAULT_STRUCTURE_FILENAME = "POSCAR"
DEFAULT_ACF_FILENAME = "ACF.dat"
DEFAULT_POTCAR_FILENAME = "POTCAR"

DEFAULT_SURFACE_CHARGE_CSV_NAME = "surface_charge.csv"
DEFAULT_SURFACE_CHARGE_PNG_NAME = "surface_charge.png"
