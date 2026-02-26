"""Public package interface for `src.structure`."""

from __future__ import annotations

from .utils import DEFAULT_METAL_SYMBOLS
from .utils import DEFAULT_THETA_BIN_DEG
from .utils import DEFAULT_WATER_OH_CUTOFF_A
from .utils import DEFAULT_Z_BIN_WIDTH_A
from .utils import WATER_MOLAR_MASS_G_PER_MOL
from .utils import Layer
from .utils import SurfaceDetectionResult
from .utils import SurfaceGeometryError
from .utils import WaterTopologyError
from .utils import detect_interface_layers
from .utils import detect_water_molecule_indices
from .utils import format_detection_summary
from .utils import get_water_oxygen_indices_array

__all__ = [
    "Layer",
    "SurfaceDetectionResult",
    "SurfaceGeometryError",
    "WaterTopologyError",
    "detect_interface_layers",
    "format_detection_summary",
    "detect_water_molecule_indices",
    "get_water_oxygen_indices_array",
    "DEFAULT_METAL_SYMBOLS",
    "DEFAULT_Z_BIN_WIDTH_A",
    "DEFAULT_THETA_BIN_DEG",
    "DEFAULT_WATER_OH_CUTOFF_A",
    "WATER_MOLAR_MASS_G_PER_MOL",
]
