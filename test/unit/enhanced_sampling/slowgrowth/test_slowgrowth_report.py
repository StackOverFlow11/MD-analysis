"""Tests for slowgrowth_analysis_with_report wrapper + contract dispatch."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot import (
    SlowgrowthAnalysisReport,
    slowgrowth_analysis_with_report,
)

DATA = Path(__file__).resolve().parents[4] / "data_example" / "sg" / "angle"
RESTART = DATA / "slowgrowth-1.restart"
LOG = DATA / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"


@pytest.fixture(scope="module")
def data_paths():
    assert RESTART.is_file(), f"fixture restart missing: {RESTART}"
    assert LOG.is_file(), f"fixture log missing: {LOG}"
    return {"restart": RESTART, "log": LOG}


# ---------------------------------------------------------------------------
# Wrapper level
# ---------------------------------------------------------------------------


class TestSlowgrowthReportWrapper:
    def test_quick_returns_csv_and_quick_png(self, data_paths, tmp_path):
        r = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            output_dir=tmp_path, plot_style="quick",
        )
        assert isinstance(r, SlowgrowthAnalysisReport)
        assert set(r.artifacts.keys()) == {"csv", "quick_png"}

    def test_publication_returns_csv_and_publication_png(
        self, data_paths, tmp_path,
    ):
        r = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            output_dir=tmp_path, plot_style="publication",
        )
        assert set(r.artifacts.keys()) == {"csv", "publication_png"}

    def test_both_returns_csv_and_both_pngs(self, data_paths, tmp_path):
        r = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            output_dir=tmp_path, plot_style="both",
        )
        assert set(r.artifacts.keys()) == {"csv", "quick_png", "publication_png"}

    def test_invalid_plot_style_raises_value_error(
        self, data_paths, tmp_path,
    ):
        with pytest.raises(ValueError):
            slowgrowth_analysis_with_report(
                str(data_paths["restart"]), str(data_paths["log"]),
                output_dir=tmp_path, plot_style="bogus",
            )

    def test_invalid_plot_style_without_inputs_still_raises(
        self, tmp_path,
    ):
        """Input validation — must fail before any file I/O."""
        with pytest.raises(ValueError):
            slowgrowth_analysis_with_report(
                "any.restart", "any.log",
                output_dir=tmp_path, plot_style="bogus",
            )

    def test_missing_restart_raises_file_not_found(self, data_paths, tmp_path):
        with pytest.raises(FileNotFoundError):
            slowgrowth_analysis_with_report(
                "nope.restart", str(data_paths["log"]),
                output_dir=tmp_path, plot_style="quick",
            )

    def test_report_to_dict_is_json_serializable(self, data_paths, tmp_path):
        r = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            output_dir=tmp_path, plot_style="quick",
        )
        json.dumps(r.to_dict())

    def test_reverse_segment_is_reversed_true_and_preserves_barrier_step(
        self, data_paths, tmp_path,
    ):
        """Reverse analysis: barrier_step MUST equal
        ``absolute_steps[peak_idx]`` from the pre-reversal segment, NOT
        ``seg.steps[peak_idx]`` which gets reset to ``[0, N-1]`` by
        ``Slowgrowth.reversed()``.

        We lock the contract by replaying the wrapper's computation and
        asserting equality, and additionally verify that the two
        candidate indices actually differ for this fixture so that a
        regression swapping them WOULD fail this test.
        """
        import numpy as np

        from md_analysis.enhanced_sampling.slowgrowth.SlowGrowth import (
            SlowgrowthFull,
        )
        from md_analysis.utils.constants import HA_TO_EV

        # Forward baseline for is_reversed comparison.
        fwd = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            initial_step=0, final_step=200,
            output_dir=tmp_path, plot_style="quick",
        )
        assert fwd.is_reversed is False

        rev = slowgrowth_analysis_with_report(
            str(data_paths["restart"]), str(data_paths["log"]),
            initial_step=200, final_step=0,
            output_dir=tmp_path, plot_style="quick",
        )
        assert rev.is_reversed is True

        # Replay the wrapper's reverse-branch computation exactly.
        full = SlowgrowthFull.from_paths(
            str(data_paths["restart"]), str(data_paths["log"]),
        )
        pre_rev = full.segment(0, 200)
        absolute_steps = pre_rev.steps[::-1].copy()
        seg = pre_rev.reversed()
        fe_ev = seg.free_energy_au * HA_TO_EV
        peak_idx = int(np.nanargmax(fe_ev))

        expected_barrier_step = int(absolute_steps[peak_idx])
        reset_index_candidate = int(seg.steps[peak_idx])

        # Regression guard: this fixture must actually distinguish the
        # two candidate values, otherwise the equality assertion below
        # would not really lock the contract.
        assert expected_barrier_step != reset_index_candidate, (
            "test fixture degenerate: absolute step equals reset index, "
            "regression into reset-index bug could slip through"
        )

        assert rev.barrier_step == expected_barrier_step

    def test_empty_segment_raises_value_error(self, data_paths, tmp_path):
        with pytest.raises(ValueError):
            slowgrowth_analysis_with_report(
                str(data_paths["restart"]), str(data_paths["log"]),
                initial_step=100, final_step=100,
                output_dir=tmp_path, plot_style="quick",
            )
