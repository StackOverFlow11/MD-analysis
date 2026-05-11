"""Charge-potential calibration sub-package.

Provides tools for calibrating the σ → φ relationship from known
data points and predicting electrode potential from surface charge density.
"""

from __future__ import annotations

from ._data import (
    CalibrationData,
    load_calibration_csv,
    load_calibration_json,
    save_calibration_json,
)
from ._mapper import (
    ChargePotentialMapper,
    DifferentialCapacitanceMapper,
    FittingInfo,
    LinearMapper,
    PolynomialMapper,
    SplineMapper,
    create_mapper,
    mapper_from_dict,
)
from ._plot import plot_calibration
from .CalibrationWorkflow import calibrate, convert_reference, predict_potential

__all__ = [
    "CalibrationData",
    "load_calibration_csv",
    "load_calibration_json",
    "save_calibration_json",
    "ChargePotentialMapper",
    "DifferentialCapacitanceMapper",
    "FittingInfo",
    "LinearMapper",
    "PolynomialMapper",
    "SplineMapper",
    "create_mapper",
    "mapper_from_dict",
    "plot_calibration",
    "calibrate",
    "convert_reference",
    "predict_potential",
]
