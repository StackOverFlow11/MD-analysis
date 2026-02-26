"""Public package interface for `src.structure.utils`."""

from __future__ import annotations

# --- Internal shared helpers (used by Analysis layer, NOT part of public API) ---
# Intentionally excluded from __all__ but explicitly collected here to centralise
# cross-layer internal dependencies.  If a private function is renamed or moved,
# only this file needs updating.
from .LayerParser import _circular_mean_fractional
from .WaterParser import _compute_bisector_cos_theta_vec
from .WaterParser import _oxygen_to_hydrogen_map
from .WaterParser import _theta_bin_count_from_ndeg

from .config import DEFAULT_METAL_SYMBOLS
from .config import DEFAULT_THETA_BIN_DEG
from .config import DEFAULT_WATER_OH_CUTOFF_A
from .config import DEFAULT_Z_BIN_WIDTH_A
from .config import TRANSITION_METAL_SYMBOLS
from .config import WATER_MOLAR_MASS_G_PER_MOL
from .LayerParser import Layer
from .LayerParser import SurfaceDetectionResult
from .LayerParser import SurfaceGeometryError
from .LayerParser import detect_interface_layers
from .LayerParser import format_detection_summary
from .WaterParser import WaterTopologyError
from .WaterParser import detect_water_molecule_indices
from .WaterParser import get_water_oxygen_indices_array

__all__ = [
    "TRANSITION_METAL_SYMBOLS",
    "DEFAULT_METAL_SYMBOLS",
    "DEFAULT_Z_BIN_WIDTH_A",
    "DEFAULT_THETA_BIN_DEG",
    "DEFAULT_WATER_OH_CUTOFF_A",
    "WATER_MOLAR_MASS_G_PER_MOL",
    "Layer",
    "SurfaceDetectionResult",
    "SurfaceGeometryError",
    "detect_interface_layers",
    "format_detection_summary",
    "WaterTopologyError",
    "detect_water_molecule_indices",
    "get_water_oxygen_indices_array",
]
