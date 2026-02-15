"""Configuration for structure-layer analysis."""

from __future__ import annotations

# Transition-metal symbols used as the default metal set in structure analysis.
# To avoid group-3 convention ambiguity, this set includes both La/Ac and Lu/Lr.
TRANSITION_METAL_SYMBOLS: tuple[str, ...] = (
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "La",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Lu",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Ac",
    "Lr",
)

# Current project default: treat all transition metals as slab-metal candidates.
DEFAULT_METAL_SYMBOLS: tuple[str, ...] = TRANSITION_METAL_SYMBOLS

# Default z-axis bin width for profile/distribution calculations (Angstrom).
DEFAULT_Z_BIN_WIDTH_A: float = 0.25

# Default theta bin width for angular probability-density profiles (degree).
DEFAULT_THETA_BIN_DEG: float = 5.0

# Maximum O-H distance used to identify water connectivity (Angstrom).
DEFAULT_WATER_OH_CUTOFF_A: float = 1.25

# Water molar mass (g/mol), used for converting water counts to mass density.
WATER_MOLAR_MASS_G_PER_MOL: float = 18.01528
