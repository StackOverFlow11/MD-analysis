"""Tests for utils._frame_discovery shared helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from md_analysis.utils._frame_discovery import (
    FRAME_DIR_STEP_TIME_RE,
    discover_frame_dirs,
    extract_step_time_from_dirname,
)


class TestExtractStepTime:
    def test_basic_bader(self):
        assert extract_step_time_from_dirname("bader_t50_i10") == (10, 50)

    def test_basic_potential(self):
        assert extract_step_time_from_dirname("potential_t1000_i50") == (50, 1000)

    def test_zero_values(self):
        assert extract_step_time_from_dirname("bader_t0_i0") == (0, 0)

    def test_no_match_returns_none(self):
        assert extract_step_time_from_dirname("no_match") is None
        assert extract_step_time_from_dirname("bader_t50") is None  # missing _i
        assert extract_step_time_from_dirname("bader_i10") is None  # missing _t
        assert extract_step_time_from_dirname("") is None

    def test_regex_anywhere_in_name(self):
        # The regex uses `search`, so it matches anywhere in the name
        assert extract_step_time_from_dirname("prefix_t7_i3_suffix") == (3, 7)

    def test_exported_regex(self):
        assert FRAME_DIR_STEP_TIME_RE.search("bader_t50_i10") is not None


class TestDiscoverFrameDirs:
    def test_numeric_sort_same_step(self, tmp_path: Path):
        """When all dirs share i=0, sort by t lexically: 5 < 50 < 200 < 1000."""
        for t in [50, 1000, 200, 5]:
            (tmp_path / f"bader_t{t}_i0").mkdir()
        result = discover_frame_dirs(tmp_path, "bader_t*_i*")
        names = [p.name for p in result]
        assert names == [
            "bader_t5_i0",
            "bader_t50_i0",
            "bader_t200_i0",
            "bader_t1000_i0",
        ]

    def test_numeric_sort_monotonic_t_i(self, tmp_path: Path):
        """Real trajectory dirs: both t and i monotonic."""
        for i, t in [(0, 0), (10, 50), (20, 100), (30, 150)]:
            (tmp_path / f"bader_t{t}_i{i}").mkdir()
        result = discover_frame_dirs(tmp_path, "bader_t*_i*")
        steps = [extract_step_time_from_dirname(p.name)[0] for p in result]
        assert steps == [0, 10, 20, 30]

    def test_drops_non_matching_dirs(self, tmp_path: Path):
        (tmp_path / "bader_t5_i0").mkdir()
        (tmp_path / "no_pattern_here").mkdir()  # skipped
        (tmp_path / "bader_t10_i1").mkdir()
        result = discover_frame_dirs(tmp_path, "*")
        names = [p.name for p in result]
        assert names == ["bader_t5_i0", "bader_t10_i1"]

    def test_drops_files(self, tmp_path: Path):
        (tmp_path / "bader_t5_i0").mkdir()
        (tmp_path / "bader_t10_i1").touch()  # file, not directory
        result = discover_frame_dirs(tmp_path, "bader_t*_i*")
        assert [p.name for p in result] == ["bader_t5_i0"]

    def test_required_true_raises_on_empty(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="_t<int>_i<int>"):
            discover_frame_dirs(tmp_path, "bader_t*_i*", required=True)

    def test_required_false_returns_empty(self, tmp_path: Path):
        assert discover_frame_dirs(tmp_path, "bader_t*_i*", required=False) == []

    def test_potential_pattern(self, tmp_path: Path):
        """Works identically for the potential_t*_i* convention."""
        for i in [0, 10, 20]:
            (tmp_path / f"potential_t{i * 5}_i{i}").mkdir()
        result = discover_frame_dirs(tmp_path, "potential_t*_i*")
        assert len(result) == 3
