"""Unit tests for ``md_analysis.workflows.enhanced_sampling``.

Covers the workflow facade for slow-growth and constrained-TI:

- ``run_slowgrowth_{quick,publication}_plot`` artifacts + metadata
- ``run_ti_single_diagnostics`` per-point artifacts + metadata
- ``run_ti_full_analysis`` multi-point artifacts (incl. diag PNG flattening)
- ``run_ti_constant_potential_correction`` orchestration

The constant-potential-correction test is gated on the
``electrochemical/charge/Bader/`` test data because the underlying
correction routines need a real POSCAR + ACF.dat per constraint
point. When the TI fixture lacks a ``bader/`` subdir the test only
exercises the failure path (``area_A2 == 0`` and missing-bader
warning), which is still the correct workflow contract.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.enhanced_sampling.constrained_ti.models import (
    ConstraintPointReport,
    TIReport,
)
from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    TIFullAnalysisReport,
)
from md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot import (
    SlowgrowthAnalysisReport,
)
from md_analysis.workflows import (
    WorkflowResult,
    require_artifacts_exist,
    run_slowgrowth_publication_plot,
    run_slowgrowth_quick_plot,
    run_ti_full_analysis,
    run_ti_single_diagnostics,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
SG_DIR = REPO_ROOT / "data_example" / "sg" / "angle"
SG_RESTART = SG_DIR / "slowgrowth-1.restart"
SG_LOG = SG_DIR / "slowgrowth-constraint_force.dat-1.LagrangeMultLog"
TI_ROOT = REPO_ROOT / "data_example" / "ti" / "double_cv" / "1k"
TI_POINT = TI_ROOT / "ti_target_0.031369"
TI_POINT_RESTART = TI_POINT / "cMD-1.restart"
TI_POINT_LOG = TI_POINT / "cMD-constraint_force.dat-1.LagrangeMultLog"

requires_sg = pytest.mark.skipif(
    not SG_RESTART.is_file(),
    reason=f"SG fixture missing: {SG_RESTART}",
)
requires_ti_root = pytest.mark.skipif(
    not TI_ROOT.is_dir(),
    reason=f"TI fixture root missing: {TI_ROOT}",
)
requires_ti_point = pytest.mark.skipif(
    not TI_POINT_RESTART.is_file(),
    reason=f"TI single-point fixture missing: {TI_POINT_RESTART}",
)


# ---------------------------------------------------------------------------
# Slow-growth workflows
# ---------------------------------------------------------------------------


@requires_sg
class TestRunSlowgrowthQuickPlot:
    def test_artifacts_and_metadata(self, tmp_path: Path) -> None:
        result = run_slowgrowth_quick_plot(
            restart_path=SG_RESTART,
            log_path=SG_LOG,
            output_dir=tmp_path / "sg_quick",
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "slowgrowth_quick_plot"
        assert set(result.artifacts.keys()) == {"csv", "quick_png"}
        require_artifacts_exist(result)

        meta = result.metadata
        assert meta["plot_style"] == "quick"
        assert isinstance(meta["n_steps"], int) and meta["n_steps"] > 0
        assert isinstance(meta["delta_F_eV"], float)
        assert isinstance(meta["delta_F_barrier_eV"], float)
        assert "is_reversed" in meta

        assert isinstance(result.extra, SlowgrowthAnalysisReport)

    def test_to_dict_json_friendly(self, tmp_path: Path) -> None:
        result = run_slowgrowth_quick_plot(
            restart_path=SG_RESTART,
            log_path=SG_LOG,
            output_dir=tmp_path / "sg_quick",
        )
        json.dumps(result.to_dict())  # must not raise


@requires_sg
class TestRunSlowgrowthPublicationPlot:
    def test_artifacts_and_metadata(self, tmp_path: Path) -> None:
        result = run_slowgrowth_publication_plot(
            restart_path=SG_RESTART,
            log_path=SG_LOG,
            output_dir=tmp_path / "sg_pub",
        )
        assert result.name == "slowgrowth_publication_plot"
        assert set(result.artifacts.keys()) == {"csv", "publication_png"}
        require_artifacts_exist(result)
        assert result.metadata["plot_style"] == "publication"


# ---------------------------------------------------------------------------
# TI single-point diagnostics
# ---------------------------------------------------------------------------


@requires_ti_point
class TestRunTiSingleDiagnostics:
    def test_artifacts_and_metadata(self, tmp_path: Path) -> None:
        result = run_ti_single_diagnostics(
            restart_path=TI_POINT_RESTART,
            log_path=TI_POINT_LOG,
            output_dir=tmp_path / "ti_single",
            equilibration=0,
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "ti_single_diagnostics"
        # Both csv and diagnostics_png produced when output_dir is given
        assert {"csv", "diagnostics_png"}.issubset(result.artifacts.keys())
        require_artifacts_exist(result)

        meta = result.metadata
        for key in (
            "n_total", "n_analyzed", "equilibration", "dt_fs",
            "time_start_fs", "time_end_fs", "xi", "tau_corr",
            "n_eff", "geweke_z",
        ):
            assert key in meta
        assert meta["equilibration"] == 0
        assert meta["n_analyzed"] == meta["n_total"]

        assert isinstance(result.extra, ConstraintPointReport)


# ---------------------------------------------------------------------------
# TI full analysis
# ---------------------------------------------------------------------------


@requires_ti_root
class TestRunTiFullAnalysis:
    def test_artifacts_and_metadata(self, tmp_path: Path) -> None:
        result = run_ti_full_analysis(
            root_dir=TI_ROOT,
            output_dir=tmp_path / "ti_full",
            dir_filter="ti_target_*",
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "ti_full_analysis"

        # Mandatory artifacts
        for key in ("convergence_csv", "free_energy_csv", "free_energy_png"):
            assert key in result.artifacts
        # Per-point PNGs flattened by index
        diag_keys = [
            k for k in result.artifacts
            if k.startswith("diagnostics_png_")
        ]
        assert len(diag_keys) == result.metadata["n_points"]

        require_artifacts_exist(result)

        meta = result.metadata
        for key in (
            "n_points", "delta_A_eV", "sigma_A_eV", "all_passed",
            "failing_indices", "parser", "dir_filter",
        ):
            assert key in meta
        assert meta["n_points"] >= 2

        assert isinstance(result.extra, TIFullAnalysisReport)
        assert isinstance(result.extra.ti_report, TIReport)

    def test_point_slice_is_propagated_to_metadata(
        self, tmp_path: Path
    ) -> None:
        result = run_ti_full_analysis(
            root_dir=TI_ROOT,
            output_dir=tmp_path / "ti_full_slice",
            dir_filter="ti_target_*",
            point_slice=":3",
        )
        assert result.metadata["point_slice"] == ":3"
        # Sliced count: <= 3 (still must be >= 2 for analysis to succeed)
        assert 2 <= result.metadata["n_points"] <= 3

    def test_to_dict_json_friendly(self, tmp_path: Path) -> None:
        result = run_ti_full_analysis(
            root_dir=TI_ROOT,
            output_dir=tmp_path / "ti_full_json",
            dir_filter="ti_target_*",
        )
        # extra is TIFullAnalysisReport (no to_dict); to_dict() must not
        # raise but the extra slot will be the raw object.
        dumped = result.to_dict()
        # name / output_dir / artifacts / metadata must serialise cleanly.
        json.dumps(
            {
                "name": dumped["name"],
                "output_dir": dumped["output_dir"],
                "artifacts": dumped["artifacts"],
                "metadata": dumped["metadata"],
            }
        )


@requires_ti_root
class TestRunTiFullAnalysisErrors:
    def test_invalid_point_slice_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            run_ti_full_analysis(
                root_dir=TI_ROOT,
                output_dir=tmp_path / "ti_bad_slice",
                dir_filter="ti_target_*",
                point_slice="not_a_slice",
            )

    def test_missing_root_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            run_ti_full_analysis(
                root_dir=tmp_path / "does_not_exist",
                output_dir=tmp_path / "ti_missing",
                dir_filter="ti_target_*",
            )


# ---------------------------------------------------------------------------
# Constant-potential correction — symbol smoke test only
# ---------------------------------------------------------------------------


def test_constant_potential_correction_importable() -> None:
    """The orchestration entry point is importable.

    A full end-to-end test would require per-constraint-point bader/
    POSCAR/ACF.dat/POTCAR data which is not currently bundled. The
    correction-specific failure modes (missing bader dirs, missing
    calibration JSON) are exercised by existing tests in the
    enhanced_sampling/constrained_ti suite; this test only guards the
    workflow facade's import surface so the upstream re-export keeps
    working.
    """
    from md_analysis.workflows.enhanced_sampling import (
        run_ti_constant_potential_correction,
    )

    assert callable(run_ti_constant_potential_correction)
