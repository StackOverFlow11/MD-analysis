"""Integration regression tests for constrained TI with arctan extrapolation.

Uses real data from data_example/ti/ti_target_0.302356/ (1970 frames).
"""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "ti" / "ti_target_0.302356"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason=f"Test data not found: {DATA_DIR}"
)


class TestRegressionTiTarget:
    """Regression tests on known-good data (1970 frames, τ_corr ≈ 8.2)."""

    def _run_diagnostics(self, output_dir=None):
        from md_analysis.enhanced_sampling.constrained_ti.workflow import (
            standalone_diagnostics,
        )

        restart = str(DATA_DIR / "cMD-1_1500.restart")
        log = str(DATA_DIR / "cMD-constraint_force.dat-1.LagrangeMultLog")
        return standalone_diagnostics(
            restart, log, output_dir=output_dir
        )

    def test_arctan_reliable_and_matches_acf(self) -> None:
        """#1 (P0): arctan reliable, sem_asymptote ≈ sem_auto, sem_auto ≈ 2.84e-04."""
        result = self._run_diagnostics()
        report = result["report"]

        # ACF SEM should be close to known value
        assert report.autocorr.sem_auto == pytest.approx(2.84e-04, rel=0.05)

        # Arctan should be reliable
        assert report.block_avg.arctan is not None
        assert report.block_avg.arctan.reliable == True  # noqa: E712

        # Arctan asymptote should match ACF SEM (10% tolerance due to
        # dense sampling pulling fit toward low-B region)
        assert report.block_avg.arctan.sem_asymptote == pytest.approx(
            report.autocorr.sem_auto, rel=0.10
        )

    def test_diagnostics_png_generated(self, tmp_path) -> None:
        """#2 (P1): PNG file generated without error."""
        result = self._run_diagnostics(output_dir=tmp_path)
        png_path = result.get("diagnostics_png")
        assert png_path is not None
        assert Path(png_path).exists()
        assert Path(png_path).stat().st_size > 0

    def test_plot_with_arctan_none(self, tmp_path) -> None:
        """#3 (P1): plot works when arctan=None (mock scipy unavailable)."""
        from unittest.mock import patch

        from md_analysis.enhanced_sampling.constrained_ti.models import (
            BlockAverageResult,
        )
        from md_analysis.enhanced_sampling.constrained_ti.plot import (
            plot_point_diagnostics,
        )
        from md_analysis.enhanced_sampling.constrained_ti.workflow import (
            analyze_standalone,
        )

        # Get a real report first
        import numpy as np

        rng = np.random.default_rng(42)
        series = rng.normal(0, 1, 2000)

        # Patch block_average to return arctan=None
        def patched_ba(series, **kwargs):
            real_module = __import__(
                "md_analysis.enhanced_sampling.constrained_ti.analysis.block_average",
                fromlist=["analyze_block_average"],
            )
            real_result = real_module.analyze_block_average(
                series, dense_sampling=False, **{
                    k: v for k, v in kwargs.items() if k != "dense_sampling"
                }
            )
            return BlockAverageResult(
                block_sizes=real_result.block_sizes,
                sem_curve=real_result.sem_curve,
                sem_plateau=real_result.sem_plateau,
                sem_at_max_B=real_result.sem_at_max_B,
                plateau_rtol=real_result.plateau_rtol,
                plateau_reached=real_result.plateau_reached,
                cross_valid_ok=real_result.cross_valid_ok,
                passed=real_result.passed,
                arctan=None,
            )

        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.workflow.analyze_block_average",
            side_effect=patched_ba,
        ):
            report = analyze_standalone(series)

        # Plot should not crash with arctan=None
        png_path = plot_point_diagnostics(report, output_dir=tmp_path)
        assert png_path.exists()
