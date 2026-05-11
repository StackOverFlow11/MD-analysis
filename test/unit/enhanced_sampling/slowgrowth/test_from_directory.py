"""Tests for SlowgrowthFull.from_directory (parser-driven entry point)."""

from __future__ import annotations

from pathlib import Path

import pytest

from md_analysis.enhanced_sampling._parsers import (
    CP2KParser,
    ParserInferenceError,
)
from md_analysis.enhanced_sampling.slowgrowth.SlowGrowth import SlowgrowthFull


REPO_ROOT = Path(__file__).resolve().parents[4]
SG_ANGLE = REPO_ROOT / "data_example" / "sg" / "angle"


_skip = pytest.mark.skipif(
    not (SG_ANGLE.is_dir()), reason="SG angle fixture missing",
)


@_skip
def test_from_directory_auto_matches_from_paths():
    """from_directory(parser='auto') must produce identical numerical
    output to the legacy from_paths() entry point."""
    restart = SG_ANGLE / "slowgrowth-1.restart"
    log = SG_ANGLE / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
    legacy = SlowgrowthFull.from_paths(str(restart), str(log))
    new = SlowgrowthFull.from_directory(SG_ANGLE)
    assert new.steps.shape == legacy.steps.shape
    assert (new.lagrange_shake == legacy.lagrange_shake).all()
    assert (new.free_energy_au == legacy.free_energy_au).all()
    assert new.timestep_fs == legacy.timestep_fs


@_skip
def test_from_directory_explicit_parser_instance():
    sg = SlowgrowthFull.from_directory(SG_ANGLE, parser=CP2KParser())
    assert sg.steps.size > 0


@_skip
def test_from_directory_explicit_parser_name():
    sg = SlowgrowthFull.from_directory(SG_ANGLE, parser="cp2k")
    assert sg.steps.size > 0


def test_from_directory_unrecognised_raises(tmp_path):
    with pytest.raises(ParserInferenceError):
        SlowgrowthFull.from_directory(tmp_path, parser="auto")
