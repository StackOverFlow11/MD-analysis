"""StructureParser — 结构解析子包（周期聚类、金属层识别、水分子拓扑）。"""
from __future__ import annotations

from .ClusterUtils import cluster_1d_periodic, find_largest_gap_periodic, gap_midpoint_periodic
from .LayerParser import (Layer, SurfaceDetectionResult, SurfaceGeometryError,
                          circular_mean_fractional, detect_interface_layers,
                          format_detection_summary, mic_delta_fractional)
from .WaterParser import (WaterTopologyError, detect_water_molecule_indices,
                          get_water_oxygen_indices_array)

__all__ = [
    "cluster_1d_periodic", "find_largest_gap_periodic", "gap_midpoint_periodic",
    "Layer", "SurfaceDetectionResult", "SurfaceGeometryError",
    "circular_mean_fractional", "detect_interface_layers",
    "format_detection_summary", "mic_delta_fractional",
    "WaterTopologyError", "detect_water_molecule_indices", "get_water_oxygen_indices_array",
]
