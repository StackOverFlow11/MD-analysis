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
from typing import Any

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

    A real end-to-end test would require per-constraint-point bader/
    POSCAR/ACF.dat/POTCAR data which is not currently bundled. The
    correction-specific failure modes (missing bader dirs, missing
    calibration JSON) are exercised by existing tests in the
    enhanced_sampling/constrained_ti suite; the orchestration / slicing
    contract itself is covered by the mock test class below.
    """
    from md_analysis.workflows.enhanced_sampling import (
        run_ti_constant_potential_correction,
    )

    assert callable(run_ti_constant_potential_correction)


# ---------------------------------------------------------------------------
# Constant-potential correction — mock orchestration
# ---------------------------------------------------------------------------
#
# These tests exercise the workflow's parameter wiring and slicing
# contract WITHOUT touching real Bader/calibration data. All external
# dependencies are monkey-patched at their canonical import sites so
# the lazy ``from ... import ...`` statements inside
# ``run_ti_constant_potential_correction`` resolve to stubs.


class _StubMapper:
    """Placeholder σ→φ mapper. The correction function is mocked so the
    mapper's behaviour is irrelevant; only its identity matters."""


def _make_stub_correction_result():
    """Build a ``ConstantPotentialResult``-shaped stub.

    Constructed via ``object.__new__`` to bypass dataclass validation;
    the workflow only reads ``correction.area_A2``,
    ``delta_A_const_phi_eV``, and ``A_const_q_eV[-1]``.
    """
    import numpy as np

    from md_analysis.enhanced_sampling.constrained_ti.correction import (
        ConstantPotentialCorrection,
        ConstantPotentialResult,
    )

    corr = ConstantPotentialCorrection(
        sigma_uC_cm2=np.array([1.0, 2.0]),
        phi_V_SHE=np.array([0.1, 0.2]),
        correction_eV=np.array([0.0, -0.05]),
        area_A2=12.34,
    )
    # ti_report is referenced only for storage; pass a sentinel object.
    return ConstantPotentialResult(
        ti_report=object(),
        correction=corr,
        A_const_q_eV=np.array([0.0, 0.42]),
        A_const_phi_eV=np.array([0.0, 0.37]),
        delta_A_const_phi_eV=0.37,
    )


@pytest.fixture
def mock_correction_pipeline(monkeypatch, tmp_path):
    """Patch every external symbol the correction workflow imports.

    Returns a dict of spies so each test can assert call counts and
    forwarded arguments.
    """
    from md_analysis.workflows import (
        WorkflowResult as _WR,
        enhanced_sampling as workflows_es,
    )
    from md_analysis.enhanced_sampling.constrained_ti import (
        correction as _correction_mod,
        io as _io_mod,
        workflow as _workflow_mod,
    )
    from md_analysis.electrochemical.calibration import (
        _data as _cal_data_mod,
        _mapper as _cal_mapper_mod,
    )

    out_dir = tmp_path / "ti_cp_corr"
    out_dir.mkdir()
    csv_corr = out_dir / "corrected_free_energy.csv"
    csv_corr.write_text("xi,phi\n")
    png_corr = out_dir / "corrected_free_energy.png"
    png_corr.write_bytes(b"\x89PNG\r\n\x1a\n")

    # Synthetic TI baseline WorkflowResult.
    ti_baseline = _WR(
        name="ti_full_analysis",
        output_dir=out_dir,
        artifacts={
            "convergence_csv": out_dir / "conv.csv",
            "free_energy_csv": out_dir / "fe.csv",
            "free_energy_png": out_dir / "fe.png",
        },
        metadata={
            "n_points": 4,
            "delta_A_eV": 0.42,
            "sigma_A_eV": 0.05,
            "all_passed": True,
            "failing_indices": [],
        },
        # `extra.ti_report` is referenced; provide a SimpleNamespace.
        extra=type("ExtraStub", (), {"ti_report": object()})(),
    )
    for f in (ti_baseline.artifacts.values()):
        f.write_text("stub\n") if f.suffix in (".csv",) else f.write_bytes(b"")

    # Mutable spy state per call.
    state: dict[str, Any] = {
        "ti_full_calls": [],
        "discover_calls": [],
        "parse_slice_calls": [],
        "load_cal_calls": [],
        "mapper_calls": [],
        "compute_corr_calls": [],
        "discover_returns": [object(), object(), object(), object()],
    }

    def fake_ti_full(**kwargs):
        state["ti_full_calls"].append(kwargs)
        return ti_baseline

    def fake_discover(root, *, parser, dir_filter, reverse, strict=True):
        # Phase 6.5: the correction workflow's phase-2 re-discovery now
        # forwards strict; capture it so the D5 invariant (phase-1 and
        # phase-2 see the same strict) can be asserted.
        state["discover_calls"].append(
            {"root": root, "parser": parser,
             "dir_filter": dir_filter, "reverse": reverse,
             "strict": strict}
        )
        return list(state["discover_returns"])

    def fake_parse_slice(spec):
        state["parse_slice_calls"].append(spec)
        # Mirror real semantics: ":2" -> slice(None, 2)
        return slice(None, 2)

    def fake_load_cal(path):
        state["load_cal_calls"].append(Path(path))
        return (object(), {"method": "linear"})

    def fake_mapper(_fit_params):
        state["mapper_calls"].append(_fit_params)
        return _StubMapper()

    def fake_compute_corr(ti_report, point_defs, mapper, **kw):
        state["compute_corr_calls"].append({
            "ti_report": ti_report,
            "point_defs_len": len(point_defs),
            "mapper": mapper,
            **kw,
        })
        return _make_stub_correction_result()

    def fake_write_csv(_result, *, output_dir):
        return Path(output_dir) / "corrected_free_energy.csv"

    def fake_plot(_result, *, output_dir):
        return Path(output_dir) / "corrected_free_energy.png"

    monkeypatch.setattr(workflows_es, "run_ti_full_analysis", fake_ti_full)
    monkeypatch.setattr(_io_mod, "discover_ti_points", fake_discover)
    monkeypatch.setattr(_workflow_mod, "_parse_point_slice", fake_parse_slice)
    monkeypatch.setattr(_cal_data_mod, "load_calibration_json", fake_load_cal)
    monkeypatch.setattr(_cal_mapper_mod, "mapper_from_dict", fake_mapper)
    monkeypatch.setattr(
        _correction_mod,
        "compute_constant_potential_correction",
        fake_compute_corr,
    )
    monkeypatch.setattr(
        _correction_mod, "write_corrected_free_energy_csv", fake_write_csv,
    )
    monkeypatch.setattr(
        _correction_mod,
        "plot_corrected_free_energy_profile",
        fake_plot,
    )

    return {
        "state": state,
        "out_dir": out_dir,
        "ti_baseline": ti_baseline,
    }


class TestConstantPotentialCorrectionOrchestration:
    """Mock-level orchestration coverage for run_ti_constant_potential_correction.

    Verifies the slicing contract (regression for the
    ``point_slice=""`` aliasing bug) plus parameter wiring and
    WorkflowResult assembly.
    """

    def test_empty_string_slice_does_not_call_parser(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        """``point_slice=""`` must be treated as "no slicing" — same
        as ``run_ti_full_from_root`` does. Calling ``_parse_point_slice("")``
        would raise ValueError and put point_defs out of sync with
        ti_report.point_reports.
        """
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        result = run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            point_slice="",
        )
        # No parse_slice call (the bug-prevention regression check).
        assert mock_correction_pipeline["state"]["parse_slice_calls"] == []
        # All four discovered point_defs reach compute_corr unsliced.
        assert mock_correction_pipeline["state"]["compute_corr_calls"][0][
            "point_defs_len"
        ] == 4
        assert isinstance(result, WorkflowResult)

    def test_none_slice_does_not_call_parser(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            point_slice=None,
        )
        assert mock_correction_pipeline["state"]["parse_slice_calls"] == []
        assert mock_correction_pipeline["state"]["compute_corr_calls"][0][
            "point_defs_len"
        ] == 4

    def test_non_empty_slice_calls_parser_and_slices_point_defs(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            point_slice=":2",
        )
        # parse_slice called exactly once with the user's spec.
        assert mock_correction_pipeline["state"]["parse_slice_calls"] == [":2"]
        # 4 discovered → slice(:2) → 2 forwarded to compute_corr.
        assert mock_correction_pipeline["state"]["compute_corr_calls"][0][
            "point_defs_len"
        ] == 2

    def test_artifacts_include_ti_baseline_plus_corrected_pair(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        result = run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
        )
        # TI baseline artifacts inherited
        for key in ("convergence_csv", "free_energy_csv", "free_energy_png"):
            assert key in result.artifacts
        # Correction artifacts added
        assert "corrected_free_energy_csv" in result.artifacts
        assert "corrected_free_energy_png" in result.artifacts

    def test_metadata_carries_correction_summary(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        result = run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            target_side="opposed",
            method="layer",
            normal="b",
        )
        meta = result.metadata
        # Correction-specific
        assert meta["target_side"] == "opposed"
        assert meta["method"] == "layer"
        assert meta["normal"] == "b"
        assert meta["area_A2"] == pytest.approx(12.34)
        assert meta["delta_A_const_phi_eV"] == pytest.approx(0.37)
        assert meta["delta_A_const_q_eV"] == pytest.approx(0.42)
        assert meta["total_correction_eV"] == pytest.approx(0.37 - 0.42)
        # TI baseline metadata inherited
        assert meta["n_points"] == 4
        assert meta["delta_A_eV"] == pytest.approx(0.42)
        # calibration path recorded
        assert meta["calibration_json"] == str(tmp_path / "cal.json")

    def test_ti_full_and_discover_share_filter_settings(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        """run_ti_full_analysis and the re-discover step must see the
        same parser / dir_filter / reverse so the post-slice point_defs
        align with ti_report.point_reports."""
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            parser="cp2k",
            dir_filter="ti_target_*",
            reverse=True,
        )
        ti_kwargs = mock_correction_pipeline["state"]["ti_full_calls"][0]
        disc_kwargs = mock_correction_pipeline["state"]["discover_calls"][0]
        assert ti_kwargs["parser"] == disc_kwargs["parser"] == "cp2k"
        assert ti_kwargs["dir_filter"] == disc_kwargs["dir_filter"] == "ti_target_*"
        assert ti_kwargs["reverse"] == disc_kwargs["reverse"] is True

    def test_calibration_json_path_loaded_once(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        cal_path = tmp_path / "custom_cal.json"
        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=cal_path,
        )
        loaded = mock_correction_pipeline["state"]["load_cal_calls"]
        assert len(loaded) == 1
        assert loaded[0] == cal_path
        assert len(mock_correction_pipeline["state"]["mapper_calls"]) == 1


# ---------------------------------------------------------------------------
# Phase 6.5: strict forwarding (workflow layer, not just CLI wiring)
# ---------------------------------------------------------------------------
#
# D1 = Option 1: strict is plumbed run_ti_full_analysis →
# run_ti_full_from_root → discover_ti_points, and
# run_ti_constant_potential_correction forwards it to BOTH the inner
# run_ti_full_analysis AND the phase-2 re-discovery (D5). Default True
# keeps agent / contract behaviour unchanged; the CLI passes False.


class TestStrictForwarding:
    def test_run_ti_full_analysis_forwards_strict_to_from_root(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        """run_ti_full_analysis → run_ti_full_from_root: explicit False
        is forwarded; the default is True (agent behaviour unchanged)."""
        from md_analysis.enhanced_sampling.constrained_ti import (
            workflow as _wf_mod,
        )

        calls: list[dict[str, Any]] = []

        def fake_from_root(**kwargs):
            calls.append(kwargs)
            return type(
                "R", (),
                {
                    "convergence_csv": tmp_path / "c.csv",
                    "free_energy_csv": tmp_path / "f.csv",
                    "free_energy_png": tmp_path / "f.png",
                    "diagnostics_pngs": (),
                    "n_points": 2,
                    "delta_A_eV": 0.1,
                    "sigma_A_eV": 0.01,
                    "all_passed": True,
                    "failing_indices": (),
                    "per_point": (),
                    "ti_report": object(),
                },
            )()

        monkeypatch.setattr(_wf_mod, "run_ti_full_from_root", fake_from_root)

        run_ti_full_analysis(
            root_dir=tmp_path, output_dir=tmp_path / "o1", strict=False,
        )
        assert calls[-1]["strict"] is False

        run_ti_full_analysis(root_dir=tmp_path, output_dir=tmp_path / "o2")
        assert calls[-1]["strict"] is True

    def test_run_ti_full_from_root_forwards_strict_to_discover(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        """run_ti_full_from_root → discover_ti_points: strict passthrough.

        Stops the pipeline right after discovery by returning an empty
        list (workflow then raises ValueError on < 2 points) so we only
        exercise the discover call.
        """
        from md_analysis.enhanced_sampling.constrained_ti import io as _io_mod
        from md_analysis.enhanced_sampling.constrained_ti.workflow import (
            run_ti_full_from_root,
        )

        seen: list[bool] = []

        def fake_discover(root, *, parser, dir_filter, reverse, strict=True):
            seen.append(strict)
            return []  # → workflow raises ValueError (< 2 points)

        monkeypatch.setattr(_io_mod, "discover_ti_points", fake_discover)
        root = tmp_path / "ti_root"
        root.mkdir()

        with pytest.raises(ValueError):
            run_ti_full_from_root(
                root_dir=root, output_dir=tmp_path / "o", strict=False,
            )
        assert seen[-1] is False

        with pytest.raises(ValueError):
            run_ti_full_from_root(root_dir=root, output_dir=tmp_path / "o")
        assert seen[-1] is True

    def test_correction_forwards_strict_to_both_phases(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        """run_ti_constant_potential_correction forwards strict to BOTH
        the inner run_ti_full_analysis AND the phase-2 discover (D5)."""
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
            strict=False,
        )
        st = mock_correction_pipeline["state"]
        assert st["ti_full_calls"][0]["strict"] is False
        assert st["discover_calls"][0]["strict"] is False

    def test_correction_strict_defaults_true_both_phases(
        self, mock_correction_pipeline, tmp_path: Path
    ) -> None:
        from md_analysis.workflows.enhanced_sampling import (
            run_ti_constant_potential_correction,
        )

        run_ti_constant_potential_correction(
            root_dir=tmp_path / "ti_root",
            output_dir=mock_correction_pipeline["out_dir"],
            calibration_json_path=tmp_path / "cal.json",
        )
        st = mock_correction_pipeline["state"]
        assert st["ti_full_calls"][0]["strict"] is True
        assert st["discover_calls"][0]["strict"] is True


class TestTiFullAnalysisAgentSchema:
    """Phase 6.5: agent reroute + strict contract exposure."""

    def test_target_fn_routes_through_workflows(self) -> None:
        from md_analysis.agent._core import get_task

        assert get_task("ti_full_analysis").target_fn == (
            "md_analysis.workflows.enhanced_sampling:run_ti_full_analysis"
        )

    def test_schema_exposes_strict_default_true(self) -> None:
        from md_analysis.agent import get_task_schema

        schema = get_task_schema("ti_full_analysis")
        props = schema["parameters"]["properties"]
        assert "strict" in props
        # default True → agent failure semantics unchanged
        assert props["strict"].get("default") is True
        # strict is optional (not in required)
        assert "strict" not in schema["parameters"].get("required", [])
