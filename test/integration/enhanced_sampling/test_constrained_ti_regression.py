"""Integration regression tests for constrained TI with F&P block averaging.

Uses real data from data_example/ti/ti_target_0.302356/ (3112 frames).
"""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "ti" / "ti_target_0.302356"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason=f"Test data not found: {DATA_DIR}"
)


class TestRegressionTiTarget:
    """Regression tests on known-good data (3112 frames, τ_corr ≈ 9.0)."""

    def _run_diagnostics(self, output_dir=None):
        from md_analysis.enhanced_sampling.constrained_ti.workflow import (
            standalone_diagnostics,
        )

        restart = str(DATA_DIR / "cMD-1_1500.restart")
        log = str(DATA_DIR / "cMD-constraint_force.dat-1.LagrangeMultLog")
        return standalone_diagnostics(
            restart, log, output_dir=output_dir
        )

    def test_fp_plateau_and_matches_acf(self) -> None:
        """F&P plateau detected, SEM_block ≈ SEM_auto."""
        result = self._run_diagnostics()
        report = result["report"]

        # ACF SEM should be close to known value
        assert report.autocorr.sem_auto == pytest.approx(2.42e-04, rel=0.05)

        # F&P plateau should be detected
        assert report.block_avg.plateau_reached == True  # noqa: E712
        assert report.block_avg.plateau_block_size is not None

        # SEM_block should agree with SEM_auto (20% tolerance)
        assert report.block_avg.plateau_sem == pytest.approx(
            report.autocorr.sem_auto, rel=0.20
        )

        # sem_final should be the plateau value
        assert report.sem_final == pytest.approx(
            report.block_avg.plateau_sem, rel=1e-10
        )

    def test_diagnostics_png_generated(self, tmp_path) -> None:
        """PNG file generated without error."""
        result = self._run_diagnostics(output_dir=tmp_path)
        png_path = result.get("diagnostics_png")
        assert png_path is not None
        assert Path(png_path).exists()
        assert Path(png_path).stat().st_size > 0

    def test_plot_with_no_plateau(self, tmp_path) -> None:
        """Plot works when plateau is not reached (mock)."""
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

        import numpy as np

        rng = np.random.default_rng(42)
        series = rng.normal(0, 1, 2000)

        def patched_ba(series, **kwargs):
            from md_analysis.enhanced_sampling.constrained_ti.analysis.block_average import (
                analyze_block_average as real_ba,
            )
            real = real_ba(series, **kwargs)
            return BlockAverageResult(
                block_sizes=real.block_sizes,
                sem_curve=real.sem_curve,
                delta_sem=real.delta_sem,
                n_total=real.n_total,
                plateau_index=None,
                plateau_sem=real.plateau_sem,
                plateau_delta=real.plateau_delta,
                plateau_block_size=None,
                plateau_reached=False,
                passed=None,
            )

        with patch(
            "md_analysis.enhanced_sampling.constrained_ti.workflow.analyze_block_average",
            side_effect=patched_ba,
        ):
            report = analyze_standalone(series)

        png_path = plot_point_diagnostics(report, output_dir=tmp_path)
        assert png_path.exists()
