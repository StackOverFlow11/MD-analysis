"""Slow-growth thermodynamic integration analysis."""

from __future__ import annotations

from .SlowGrowth import Slowgrowth, SlowgrowthFull, SlowgrowthSegment
from .SlowGrowthPlot import (
    plot_slowgrowth_quick,
    plot_slowgrowth_publication,
    slowgrowth_analysis,
    write_slowgrowth_csv,
)

__all__ = [
    "Slowgrowth",
    "SlowgrowthFull",
    "SlowgrowthSegment",
    "plot_slowgrowth_quick",
    "plot_slowgrowth_publication",
    "write_slowgrowth_csv",
    "slowgrowth_analysis",
]
