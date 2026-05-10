"""Tests for md_analysis.scripts._frame_selector — unified trajectory slicing."""

from __future__ import annotations

from pathlib import Path

import pytest

from md_analysis.scripts._frame_selector import (
    FrameSelection,
    FrameSelectionError,
    iter_selected_frames,
    resolve_single_frame,
)


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _write_xyz_with_time(path: Path, times_fs: list[float]) -> None:
    """Write a minimal CP2K-style XYZ trajectory with given per-frame times.

    Each frame has 2 atoms (Cu at origin, O slightly offset) — enough to be a
    valid XYZ. The comment line uses CP2K format so ``atoms.info['time']`` is
    populated on read.
    """
    lines: list[str] = []
    for step, t in enumerate(times_fs):
        lines.append("2")
        lines.append(f" i = {step:>8d}, time = {t:>12.3f}, E = -1000.0")
        lines.append("Cu  0.0 0.0 0.0")
        lines.append("O   1.0 1.0 1.0")
    path.write_text("\n".join(lines) + "\n")


def _write_xyz_no_time(path: Path, n_frames: int) -> None:
    """Write an XYZ trajectory without time metadata (bare comment line)."""
    lines: list[str] = []
    for _ in range(n_frames):
        lines.append("2")
        lines.append("no-time-here")
        lines.append("Cu  0.0 0.0 0.0")
        lines.append("O   1.0 1.0 1.0")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# FrameSelection validation
# ---------------------------------------------------------------------------


class TestFrameSelectionValidation:
    def test_defaults_are_index_mode(self):
        sel = FrameSelection()
        assert sel.mode == "index"

    def test_index_mode_valid(self):
        sel = FrameSelection(mode="index", frame_start=0, frame_end=10, frame_step=2)
        assert sel.frame_step == 2

    def test_time_mode_requires_all_three(self):
        with pytest.raises(FrameSelectionError, match="Time mode requires"):
            FrameSelection(mode="time", time_start_fs=0.0, time_end_fs=100.0)
        with pytest.raises(FrameSelectionError, match="Time mode requires"):
            FrameSelection(mode="time", time_start_fs=0.0, time_step_fs=5.0)
        with pytest.raises(FrameSelectionError, match="Time mode requires"):
            FrameSelection(mode="time", time_end_fs=100.0, time_step_fs=5.0)

    def test_time_mode_complete_is_valid(self):
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=100.0, time_step_fs=10.0,
        )
        assert sel.is_time_mode if hasattr(sel, "is_time_mode") else True

    def test_invalid_mode_string(self):
        with pytest.raises(FrameSelectionError, match="Invalid mode"):
            FrameSelection(mode="random")  # type: ignore[arg-type]

    def test_frame_step_must_be_positive(self):
        with pytest.raises(FrameSelectionError, match="frame_step"):
            FrameSelection(frame_step=0)
        with pytest.raises(FrameSelectionError, match="frame_step"):
            FrameSelection(frame_step=-1)

    def test_time_step_must_be_positive(self):
        with pytest.raises(FrameSelectionError, match="time_step_fs"):
            FrameSelection(
                mode="time", time_start_fs=0, time_end_fs=100, time_step_fs=0,
            )
        with pytest.raises(FrameSelectionError, match="time_step_fs"):
            FrameSelection(
                mode="time", time_start_fs=0, time_end_fs=100, time_step_fs=-5,
            )

    def test_time_start_must_not_exceed_end(self):
        with pytest.raises(FrameSelectionError, match="time_start_fs"):
            FrameSelection(
                mode="time", time_start_fs=100, time_end_fs=10, time_step_fs=5,
            )


# ---------------------------------------------------------------------------
# iter_selected_frames — index mode
# ---------------------------------------------------------------------------


class TestIterIndexMode:
    def test_all_frames_default(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0])
        sel = FrameSelection()  # defaults
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [0, 1, 2, 3]

    def test_start_end_slicing(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0])
        sel = FrameSelection(frame_start=1, frame_end=3)
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [1, 2]

    def test_step(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0, 20.0])
        sel = FrameSelection(frame_step=2)
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [0, 2, 4]

    def test_empty_range(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0])
        sel = FrameSelection(frame_start=5, frame_end=10)
        result = list(iter_selected_frames(xyz, sel))
        assert result == []


# ---------------------------------------------------------------------------
# iter_selected_frames — time mode
# ---------------------------------------------------------------------------


class TestIterTimeMode:
    def test_exact_step_match(self, tmp_path):
        """dt=5fs trajectory, time_step_fs=10 → every other frame."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0, 20.0])
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=20.0, time_step_fs=10.0,
        )
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [0, 2, 4]

    def test_greedy_step_larger_than_dt(self, tmp_path):
        """time_step_fs=7 with dt=5 → greedy picks frames at t=0, 10, 20."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0, 20.0, 25.0])
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=25.0, time_step_fs=7.0,
        )
        result = list(iter_selected_frames(xyz, sel))
        # target=0 → yield t=0, next=7
        # target=7 → yield t=10, next=17
        # target=17 → yield t=20, next=27
        # target=27 → next frame t=25 < 27, no; next t>25 → break
        times = [float(atoms.info["time"]) for _, atoms in result]
        assert times == [0.0, 10.0, 20.0]

    def test_step_smaller_than_dt(self, tmp_path):
        """time_step_fs=2 with dt=5 → every frame in range."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0])
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=15.0, time_step_fs=2.0,
        )
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [0, 1, 2, 3]

    def test_skip_frames_before_start(self, tmp_path):
        """Frames with time < time_start_fs are skipped."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0, 20.0])
        sel = FrameSelection(
            mode="time", time_start_fs=10.0, time_end_fs=20.0, time_step_fs=5.0,
        )
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [2, 3, 4]

    def test_terminate_at_end(self, tmp_path):
        """Frames with time > time_end_fs break iteration."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0, 20.0])
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=10.0, time_step_fs=5.0,
        )
        result = list(iter_selected_frames(xyz, sel))
        assert [idx for idx, _ in result] == [0, 1, 2]

    def test_missing_time_raises(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_no_time(xyz, n_frames=3)
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=100.0, time_step_fs=10.0,
        )
        with pytest.raises(FrameSelectionError, match="atoms.info\\['time'\\]"):
            list(iter_selected_frames(xyz, sel))


# ---------------------------------------------------------------------------
# resolve_single_frame
# ---------------------------------------------------------------------------


class TestResolveSingleFrame:
    def test_index_mode_basic(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0])
        idx, atoms, warnings = resolve_single_frame(xyz, mode="index", frame=1)
        assert idx == 1
        assert warnings == []
        assert float(atoms.info["time"]) == 5.0

    def test_index_mode_out_of_range(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0])
        with pytest.raises(FrameSelectionError, match="not found"):
            resolve_single_frame(xyz, mode="index", frame=10)

    def test_time_mode_exact_match(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0])
        idx, atoms, warnings = resolve_single_frame(xyz, mode="time", time_fs=10.0)
        assert idx == 2
        assert warnings == []  # exact match → no warning

    def test_time_mode_nearest_neighbor_warns(self, tmp_path):
        """t=7.5 not in trajectory → nearest frame (t=5 or t=10), emit warning."""
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0, 10.0, 15.0])
        idx, atoms, warnings = resolve_single_frame(xyz, mode="time", time_fs=7.5)
        # Both t=5 and t=10 are equidistant; early-break picks the first seen (t=5)
        # or whichever has strictly smaller diff. Implementation detail — we assert
        # a warning was emitted and actual time != requested.
        assert len(warnings) == 1
        assert "nearest" in warnings[0].lower()
        actual_time = float(atoms.info["time"])
        assert actual_time != 7.5

    def test_time_mode_missing_metadata(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_no_time(xyz, n_frames=2)
        with pytest.raises(FrameSelectionError, match="atoms.info\\['time'\\]"):
            resolve_single_frame(xyz, mode="time", time_fs=5.0)

    def test_time_mode_requires_time_fs(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0, 5.0])
        with pytest.raises(FrameSelectionError, match="time_fs"):
            resolve_single_frame(xyz, mode="time", time_fs=None)

    def test_invalid_mode(self, tmp_path):
        xyz = tmp_path / "traj.xyz"
        _write_xyz_with_time(xyz, [0.0])
        with pytest.raises(FrameSelectionError, match="Invalid mode"):
            resolve_single_frame(xyz, mode="random")  # type: ignore[arg-type]
