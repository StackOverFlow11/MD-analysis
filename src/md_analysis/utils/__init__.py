"""Public package interface for `src.utils`."""

from __future__ import annotations

# --- Internal shared helpers (used by water layer, NOT part of public API) ---
# Intentionally excluded from __all__ but explicitly collected here to centralise
# cross-layer internal dependencies.  If a private function is renamed or moved,
# only this file needs updating.
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
from .ClusterUtils import cluster_1d_periodic
from .ClusterUtils import find_largest_gap_periodic
from .ClusterUtils import gap_midpoint_periodic
from .CubeParser import CubeHeader
from .CubeParser import read_cube_header_and_values
from .CubeParser import slab_average_potential_ev
from .CubeParser import plane_avg_phi_z_ev
from .CubeParser import z_coords_ang
from .CubeParser import extract_step_from_cube_filename
from .config import HA_TO_EV
from .config import BOHR_TO_ANG
from .config import DP_A_H3O_W_EV
from .config import MU_HPLUS_G0_EV
from .config import DELTA_E_ZP_EV
from .BaderParser import BaderParseError
from .BaderParser import load_bader_atoms
from .CellParser import CellParseError
from .CellParser import parse_abc_from_md_inp
from .CellParser import parse_abc_from_restart
from .RestartParser import ColvarDef
from .RestartParser import ConstraintInfo
from .RestartParser import LagrangeMultLog
from .RestartParser import SlowGrowthParseError
from .RestartParser import SlowGrowthRestart
from .RestartParser import compute_target_series
from .RestartParser import parse_lagrange_mult_log
from .RestartParser import parse_slowgrowth_restart

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
    "cluster_1d_periodic",
    "find_largest_gap_periodic",
    "gap_midpoint_periodic",
    "CubeHeader",
    "read_cube_header_and_values",
    "slab_average_potential_ev",
    "plane_avg_phi_z_ev",
    "z_coords_ang",
    "extract_step_from_cube_filename",
    "HA_TO_EV",
    "BOHR_TO_ANG",
    "DP_A_H3O_W_EV",
    "MU_HPLUS_G0_EV",
    "DELTA_E_ZP_EV",
    "BaderParseError",
    "load_bader_atoms",
    "CellParseError",
    "parse_abc_from_md_inp",
    "parse_abc_from_restart",
    "ColvarDef",
    "ConstraintInfo",
    "LagrangeMultLog",
    "SlowGrowthParseError",
    "SlowGrowthRestart",
    "compute_target_series",
    "parse_lagrange_mult_log",
    "parse_slowgrowth_restart",
]
