"""Integration tests for slow-growth free-energy plotting."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from md_analysis.enhanced_sampling.slowgrowth import (
    SlowgrowthFull,
    plot_slowgrowth_quick,
    plot_slowgrowth_publication,
    slowgrowth_analysis,
    write_slowgrowth_csv,
)
from md_analysis.utils.config import HA_TO_EV

# ---------------------------------------------------------------------------
# Test data paths
# ---------------------------------------------------------------------------

_DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "sg" / "distance_combinedCV"
_RESTART = str(_DATA_DIR / "slowgrowth-1.restart")
_LOG = str(_DATA_DIR / "slowgrowth-constraint_force.dat-1.LagrangeMultLog")


@pytest.fixture
def sg_full() -> SlowgrowthFull:
    return SlowgrowthFull.from_paths(_RESTART, _LOG)


# ---------------------------------------------------------------------------
# SlowgrowthFull basics
# ---------------------------------------------------------------------------

class TestSlowgrowthFullParsing:
    def test_full_has_steps(self, sg_full: SlowgrowthFull):
        assert sg_full.n_steps > 0

    def test_arrays_aligned(self, sg_full: SlowgrowthFull):
        n = sg_full.n_steps
        assert sg_full.steps.shape == (n,)
        assert sg_full.times_fs.shape == (n,)
        assert sg_full.target_au.shape == (n,)
        assert sg_full.lagrange_shake.shape == (n,)
        assert sg_full.free_energy_au.shape == (n,)

    def test_free_energy_starts_at_zero(self, sg_full: SlowgrowthFull):
        assert sg_full.free_energy_au[0] == 0.0


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

class TestWriteCSV:
    def test_csv_columns_and_rows(self, sg_full: SlowgrowthFull, tmp_path: Path):
        csv_path = write_slowgrowth_csv(sg_full, output_dir=tmp_path)
        assert csv_path.exists()

        with csv_path.open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == sg_full.n_steps
        expected_cols = {"step", "time_fs", "target_au", "lagrange_au",
                         "free_energy_au", "free_energy_ev"}
        assert set(rows[0].keys()) == expected_cols

        # Check free_energy_ev is correct conversion
        fe_au = float(rows[-1]["free_energy_au"])
        fe_ev = float(rows[-1]["free_energy_ev"])
        assert abs(fe_ev - fe_au * HA_TO_EV) < 1e-8


# ---------------------------------------------------------------------------
# Quick plot
# ---------------------------------------------------------------------------

class TestQuickPlot:
    def test_quick_plot_creates_png(self, sg_full: SlowgrowthFull, tmp_path: Path):
        png = plot_slowgrowth_quick(sg_full, output_dir=tmp_path)
        assert png.exists()
        assert png.suffix == ".png"
        assert png.stat().st_size > 1000  # non-trivial file


# ---------------------------------------------------------------------------
# Publication plot
# ---------------------------------------------------------------------------

class TestPublicationPlot:
    def test_publication_plot_creates_png(self, sg_full: SlowgrowthFull, tmp_path: Path):
        png = plot_slowgrowth_publication(sg_full, output_dir=tmp_path)
        assert png.exists()
        assert png.suffix == ".png"
        assert png.stat().st_size > 1000


# ---------------------------------------------------------------------------
# Unified entry point — forward
# ---------------------------------------------------------------------------

class TestSlowgrowthAnalysis:
    def test_forward_full_range(self, tmp_path: Path):
        results = slowgrowth_analysis(
            _RESTART, _LOG,
            initial_step=0,
            final_step=None,
            output_dir=tmp_path,
            plot_style="both",
        )
        assert "csv" in results
        assert "quick_png" in results
        assert "publication_png" in results
        for path in results.values():
            assert path.exists()

    def test_forward_sub_segment(self, tmp_path: Path):
        results = slowgrowth_analysis(
            _RESTART, _LOG,
            initial_step=10,
            final_step=100,
            output_dir=tmp_path,
            plot_style="quick",
        )
        assert "csv" in results
        assert "quick_png" in results
        assert "publication_png" not in results

        # Check CSV has the right number of rows
        with results["csv"].open() as f:
            n_data_rows = sum(1 for _ in f) - 1  # minus header
        assert n_data_rows == 90  # [10, 100)

    def test_reversed_segment(self, tmp_path: Path):
        """initial_step > final_step triggers reversal."""
        results = slowgrowth_analysis(
            _RESTART, _LOG,
            initial_step=100,
            final_step=10,
            output_dir=tmp_path,
            plot_style="publication",
        )
        assert "csv" in results
        assert "publication_png" in results
        assert "quick_png" not in results

        with results["csv"].open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 90

        # Free energy should start at 0 after reversal
        assert float(rows[0]["free_energy_au"]) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Preview output (optional manual inspection)
# ---------------------------------------------------------------------------

class TestPreviewOutput:
    """Generate output in test/_tmp_preview/ for manual inspection.

    Not a correctness test — only verifies files are created.
    """

    def test_preview_output(self):
        preview_dir = Path(__file__).resolve().parents[2] / "_tmp_preview"
        preview_dir.mkdir(parents=True, exist_ok=True)

        results = slowgrowth_analysis(
            _RESTART, _LOG,
            output_dir=preview_dir,
            plot_style="both",
        )
        for path in results.values():
            assert path.exists()
