"""Constrained thermodynamic integration convergence analysis.

Public API
----------
- Data structures: ``models``
- Analysis engines: ``analysis``
- Workflow: ``workflow``
- Integration: ``integration``
- I/O: ``io``
- Plotting: ``plot``
"""

from __future__ import annotations

from .correction import ConstantPotentialCorrection, ConstantPotentialResult
from .models import (
    AutocorrResult,
    BlockAverageResult,
    ConstraintPointInput,
    ConstraintPointReport,
    ConvergenceError,
    GewekeResult,
    InsufficientSamplingError,
    RunningAverageResult,
    TIPointDefinition,
    TIReport,
)

__all__ = [
    # Models
    "AutocorrResult",
    "BlockAverageResult",
    "ConstantPotentialCorrection",
    "ConstantPotentialResult",
    "ConstraintPointInput",
    "ConstraintPointReport",
    "GewekeResult",
    "RunningAverageResult",
    "TIPointDefinition",
    "TIReport",
    # Exceptions
    "ConvergenceError",
    "InsufficientSamplingError",
]
