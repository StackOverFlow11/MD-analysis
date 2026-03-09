"""Shared configuration constants for structure-layer and potential analysis."""

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
DEFAULT_Z_BIN_WIDTH_A: float = 0.1

# Default theta bin width for angular probability-density profiles (degree).
DEFAULT_THETA_BIN_DEG: float = 5.0

# Maximum O-H distance used to identify water connectivity (Angstrom).
DEFAULT_WATER_OH_CUTOFF_A: float = 1.25

# Water molar mass (g/mol), used for converting water counts to mass density.
WATER_MOLAR_MASS_G_PER_MOL: float = 18.01528

# --------------------------------------------------------------------------
# Unit conversion constants (potential analysis)
# --------------------------------------------------------------------------

# Hartree to electron-volt (CODATA 2018).
HA_TO_EV: float = 27.211386245988

# Bohr radius to Angstrom.
BOHR_TO_ANG: float = 0.529177210903

# --------------------------------------------------------------------------
# Computational SHE constants (eV) — from method_intro/*
# --------------------------------------------------------------------------

# Deprotonation free energy of H3O+(aq) in water (eV).
DP_A_H3O_W_EV: float = 15.35

# Chemical potential of H+ in the gas-phase standard state (eV).
MU_HPLUS_G0_EV: float = 15.81

# Zero-point energy correction (eV).
DELTA_E_ZP_EV: float = 0.35

# --------------------------------------------------------------------------
# Axis mapping and area-vector indices
# --------------------------------------------------------------------------

# Map cell-axis label to integer index.
AXIS_MAP: dict[str, int] = {"a": 0, "b": 1, "c": 2}

# Map surface-normal axis to the two cell-vector indices spanning the surface plane.
AREA_VECTOR_INDICES: dict[str, tuple[int, int]] = {"a": (1, 2), "b": (0, 2), "c": (0, 1)}

# --------------------------------------------------------------------------
# Interface label and charge method constants
# --------------------------------------------------------------------------

INTERFACE_NORMAL_ALIGNED: str = "normal_aligned"
INTERFACE_NORMAL_OPPOSED: str = "normal_opposed"
CHARGE_METHOD_COUNTERION: str = "counterion"
CHARGE_METHOD_LAYER: str = "layer"
