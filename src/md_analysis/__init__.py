"""Public package interface for `src`."""

from __future__ import annotations

from . import utils
from . import water
from . import electrochemical
from .electrochemical import potential, charge
from .exceptions import MDAnalysisError

__version__ = "0.1.0"

__all__ = ["utils", "water", "electrochemical", "potential", "charge", "__version__", "MDAnalysisError"]

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
