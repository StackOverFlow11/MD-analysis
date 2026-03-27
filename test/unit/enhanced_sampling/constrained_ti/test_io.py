"""Unit tests for constrained_ti.io — discover_ti_points sorting."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from md_analysis.enhanced_sampling.constrained_ti.io import discover_ti_points
from md_analysis.enhanced_sampling.constrained_ti.models import TIPointDefinition


def _make_ti_dirs(tmp_path: Path, xi_values: list[float]) -> None:
    """Create dummy ti_target_<xi>/ dirs with mock restart & log files."""
    for xi in xi_values:
        d = tmp_path / f"ti_target_{xi}"
        d.mkdir()
        (d / "cMD-1.restart").touch()
        (d / "cMD-constraint_force.dat-1.LagrangeMultLog").touch()


class TestDiscoverTiPointsReverse:
    """Tests for the reverse parameter of discover_ti_points."""

    def test_default_ascending(self, tmp_path: Path) -> None:
        _make_ti_dirs(tmp_path, [0.3, -0.5, 1.0])
        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_restart",
            side_effect=lambda d: d / "cMD-1.restart",
        ), patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_log",
            side_effect=lambda d: d / "cMD-constraint_force.dat-1.LagrangeMultLog",
        ):
            points = discover_ti_points(tmp_path, pattern="ti_target")
        xis = [p.xi for p in points]
        assert xis == sorted(xis)

    def test_reverse_descending(self, tmp_path: Path) -> None:
        _make_ti_dirs(tmp_path, [0.3, -0.5, 1.0])
        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_restart",
            side_effect=lambda d: d / "cMD-1.restart",
        ), patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_log",
            side_effect=lambda d: d / "cMD-constraint_force.dat-1.LagrangeMultLog",
        ):
            points = discover_ti_points(tmp_path, pattern="ti_target", reverse=True)
        xis = [p.xi for p in points]
        assert xis == sorted(xis, reverse=True)

    def test_reverse_same_points_as_default(self, tmp_path: Path) -> None:
        _make_ti_dirs(tmp_path, [0.3, -0.5, 1.0])
        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_restart",
            side_effect=lambda d: d / "cMD-1.restart",
        ), patch(
            "md_analysis.enhanced_sampling.constrained_ti.io._find_log",
            side_effect=lambda d: d / "cMD-constraint_force.dat-1.LagrangeMultLog",
        ):
            asc = discover_ti_points(tmp_path, pattern="ti_target")
            desc = discover_ti_points(tmp_path, pattern="ti_target", reverse=True)
        assert set(p.xi for p in asc) == set(p.xi for p in desc)
        assert [p.xi for p in desc] == list(reversed([p.xi for p in asc]))
