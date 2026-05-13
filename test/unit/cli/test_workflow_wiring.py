"""Monkey-patched wiring tests for the CLI commands migrated to
``md_analysis.workflows`` during the entrance refactor.

These tests confirm three properties at the CLI execute() boundary
without exercising any scientific code:

1. The workflow facade is looked up at the expected module path
   (i.e. the ``lazy_import`` target is the new workflows path, not
   the legacy backend).
2. The CLI passes through the right kwargs from its ctx dict —
   especially the ``OUTDIR_RESOLVED.parent`` adapter for the two
   charge commands, the ``mode/frame/time_fs`` mapping built by
   ``_single_frame_workflow_kwargs`` for the scripts single
   commands, and the artifact-key reads used by the print loops.
3. The workflow's ``WorkflowResult.artifacts`` keys the CLI reads
   actually match the workflow contract (so a future refactor
   that renames an artifact key would fail this test, not the
   integration suite).

No real fixtures, no disk I/O. Each test patches the *module-level*
``lazy_import`` in the CLI module so ``execute()`` receives a stub
workflow that records its kwargs and returns a hand-crafted
``WorkflowResult``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from md_analysis.cli import _charge, _enhanced_sampling, _scripts
from md_analysis.cli._params import K
from md_analysis.workflows import WorkflowResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _LazyImportStub:
    """Replaces ``lazy_import`` so the CLI receives a captured stub.

    Records every (module_path, name) lookup and the kwargs the CLI
    forwards to the resolved callable, then returns a caller-supplied
    workflow stub.
    """

    def __init__(self, stub: Any) -> None:
        self.lookups: list[tuple[str, str]] = []
        self.calls: list[dict[str, Any]] = []
        self._stub = stub

    def __call__(self, module_path: str, name: str) -> Any:
        self.lookups.append((module_path, name))

        def _wrapped(**kwargs: Any) -> Any:
            self.calls.append(kwargs)
            return self._stub(**kwargs)

        return _wrapped


def _make_workflow_result(
    *,
    name: str,
    output_dir: Path,
    artifacts: dict[str, Path] | None = None,
) -> WorkflowResult:
    return WorkflowResult(
        name=name,
        output_dir=output_dir,
        artifacts=artifacts or {},
        metadata={},
    )


# ---------------------------------------------------------------------------
# _charge.py: OUTDIR_RESOLVED.parent adapter
# ---------------------------------------------------------------------------


class TestTrackedChargeCmd:
    """CLI 225 forwards OUTDIR_RESOLVED.parent (not OUTDIR_RESOLVED)
    so the workflow's internal ``/tracked`` append yields the
    canonical ``<base>/electrochemical/charge/tracked/`` layout."""

    def test_passes_parent_dir_and_reads_tracked_csv_artifact(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        resolved_outdir = (
            tmp_path / "analysis" / "electrochemical" / "charge" / "tracked"
        )
        resolved_outdir.mkdir(parents=True)
        expected_csv = resolved_outdir / "atomic_charges_tracked_xyz.csv"
        expected_csv.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="tracked_charge",
                output_dir=resolved_outdir,
                artifacts={"tracked_charge_csv": expected_csv},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_charge, "lazy_import", stub)

        ctx = {
            K.OUTDIR_RESOLVED: str(resolved_outdir),
            K.ROOT_DIR: ".",
            K.ATOM_INDICES_XYZ: [0, 1, 2],
            K.DIR_PATTERN: "bader_t*_i*",
            K.FRAME_START: None,
            K.FRAME_END: None,
            K.FRAME_STEP: None,
        }
        cmd = _charge.TrackedChargeCmd("225", "Tracked Atom Charges")
        cmd.execute(ctx)

        # workflow facade resolved through the canonical workflows path
        assert stub.lookups == [
            ("md_analysis.workflows.charge", "run_tracked_charge"),
        ]
        # exactly one call with the OUTDIR_RESOLVED.parent adapter
        assert len(stub.calls) == 1
        call = stub.calls[0]
        assert call["output_dir"] == Path(resolved_outdir).parent
        # ctx fields forwarded straight through
        assert call["root_dir"] == "."
        assert call["atom_indices_xyz"] == [0, 1, 2]
        assert call["dir_pattern"] == "bader_t*_i*"


class TestCounterionChargeCmd:
    """CLI 226 mirrors the parent-dir adapter for the
    ``/counterion_tracking`` workflow segment."""

    def test_passes_parent_dir_and_reads_counterion_csv_artifact(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        resolved_outdir = (
            tmp_path
            / "analysis"
            / "electrochemical"
            / "charge"
            / "counterion_tracking"
        )
        resolved_outdir.mkdir(parents=True)
        expected_csv = resolved_outdir / "counterion_charges.csv"
        expected_csv.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="counterion_charge",
                output_dir=resolved_outdir,
                artifacts={"counterion_charge_csv": expected_csv},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_charge, "lazy_import", stub)

        ctx = {
            K.OUTDIR_RESOLVED: str(resolved_outdir),
            K.ROOT_DIR: ".",
            K.METAL_ELEMENTS: None,
            K.NORMAL: "c",
            K.LAYER_TOL: 0.5,
            K.DIR_PATTERN: "bader_t*_i*",
            K.FRAME_START: None,
            K.FRAME_END: None,
            K.FRAME_STEP: None,
        }
        cmd = _charge.CounterionChargeCmd("226", "Counterion Charges")
        cmd.execute(ctx)

        assert stub.lookups == [
            ("md_analysis.workflows.charge", "run_counterion_charge"),
        ]
        assert len(stub.calls) == 1
        call = stub.calls[0]
        assert call["output_dir"] == Path(resolved_outdir).parent
        # forwarding spot-checks
        assert call["root_dir"] == "."
        assert call["normal"] == "c"
        assert call["layer_tol_A"] == 0.5


# ---------------------------------------------------------------------------
# _enhanced_sampling.py: workflow name dispatch
# ---------------------------------------------------------------------------


class TestSlowGrowthPlotCmds:
    """CLI 301/302 must dispatch to the matching SG workflow name."""

    @pytest.mark.parametrize(
        "cmd_cls, expected_name",
        [
            (_enhanced_sampling.SGQuickPlotCmd, "run_slowgrowth_quick_plot"),
            (
                _enhanced_sampling.SGPublicationPlotCmd,
                "run_slowgrowth_publication_plot",
            ),
        ],
    )
    def test_dispatches_to_correct_workflow(
        self,
        monkeypatch,
        tmp_path: Path,
        cmd_cls: type,
        expected_name: str,
    ) -> None:
        outdir = tmp_path / "out"
        outdir.mkdir()

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="slowgrowth_stub", output_dir=outdir,
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_enhanced_sampling, "lazy_import", stub)

        ctx = {
            K.RESTART_PATH: "stub.restart",
            K.LOG_PATH: "stub.LagrangeMultLog",
            K.OUTDIR_RESOLVED: str(outdir),
            K.INITIAL_STEP: 0,
            K.FINAL_STEP: None,
            K.COLVAR_ID: None,
        }
        cmd = cmd_cls("3xx", "SG plot")
        cmd.execute(ctx)

        # Module path is the workflows enhanced_sampling facade
        assert len(stub.lookups) == 1
        module_path, fn_name = stub.lookups[0]
        assert module_path == "md_analysis.workflows.enhanced_sampling"
        assert fn_name == expected_name
        # All six required kwargs forwarded
        call = stub.calls[0]
        assert call["restart_path"] == "stub.restart"
        assert call["log_path"] == "stub.LagrangeMultLog"
        assert call["output_dir"] == str(outdir)
        assert call["initial_step"] == 0
        assert call["final_step"] is None
        assert call["colvar_id"] is None

    def test_base_class_without_workflow_name_raises(self) -> None:
        """The base class has _workflow_name = "" and must reject
        direct instantiation/execute() to catch regressions where
        a subclass forgot to set it."""
        base = _enhanced_sampling._SlowgrowthPlotCmd("3yy", "base")
        with pytest.raises(RuntimeError, match="_workflow_name"):
            base.execute({})


# ---------------------------------------------------------------------------
# _scripts.py: single-frame workflow kwarg mapping + artifact reads
# ---------------------------------------------------------------------------


class TestSingleFrameWorkflowKwargs:
    """The new ``_single_frame_workflow_kwargs`` helper maps ctx's
    FRAME_MODE + FRAME/SINGLE_TIME_FS to the workflow signature."""

    def test_index_mode_emits_mode_and_frame(self) -> None:
        ctx = {K.FRAME_MODE: "index", K.FRAME: 7}
        kwargs = _scripts._single_frame_workflow_kwargs(ctx)
        assert kwargs == {"mode": "index", "frame": 7}

    def test_time_mode_emits_mode_and_time_fs(self) -> None:
        ctx = {K.FRAME_MODE: "time", K.SINGLE_TIME_FS: 250.0}
        kwargs = _scripts._single_frame_workflow_kwargs(ctx)
        assert kwargs == {"mode": "time", "time_fs": 250.0}


class TestBaderSingleCmd:
    """CLI 411 forwards XYZ + cell_abc + single-frame mode kwargs to
    run_bader_single and prints from result.artifacts['workdir']."""

    def test_dispatches_to_run_bader_single(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        workdir = tmp_path / "out" / "bader"
        workdir.mkdir(parents=True)

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="bader_single",
                output_dir=tmp_path / "out",
                artifacts={"workdir": workdir},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_scripts, "lazy_import", stub)

        ctx = {
            K.XYZ: "stub.xyz",
            K.CELL_ABC: (10.0, 10.0, 30.0),
            K.OUTDIR: str(tmp_path / "out"),
            K.WORKDIR_NAME: "bader",
            K.SCRIPT_PATH: None,
            K.GEN_POTCAR: False,
            K.FRAME_MODE: "index",
            K.FRAME: 0,
        }
        cmd = _scripts.BaderSingleCmd("411", "Bader Single")
        cmd.execute(ctx)

        # workflow path
        assert stub.lookups == [
            ("md_analysis.workflows.scripts", "run_bader_single"),
        ]
        # kwarg map includes single-frame helper output
        call = stub.calls[0]
        assert call["xyz_path"] == "stub.xyz"
        assert call["cell_abc"] == (10.0, 10.0, 30.0)
        assert call["output_dir"] == str(tmp_path / "out")
        assert call["workdir_name"] == "bader"
        assert call["script_path"] is None
        assert call["generate_potcar"] is False
        assert call["mode"] == "index"
        assert call["frame"] == 0
        # ``time_fs`` MUST be absent in index mode (regression guard
        # against accidentally forwarding the unused branch).
        assert "time_fs" not in call


class TestBaderBatchCmd:
    """CLI 412 forwards frame slicing kwargs and reads workdir_<i>
    keys from result.artifacts."""

    def test_dispatches_and_reads_workdir_indexed_artifacts(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        wd0 = out / "bader_t000_i000"
        wd1 = out / "bader_t100_i100"
        wd0.mkdir()
        wd1.mkdir()

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="bader_batch",
                output_dir=out,
                artifacts={
                    "workdir_0": wd0,
                    "workdir_1": wd1,
                    # Non-workdir keys must be ignored by the CLI's
                    # filter expression.
                    "summary_csv": out / "summary.csv",
                },
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_scripts, "lazy_import", stub)

        ctx = {
            K.XYZ: "stub.xyz",
            K.CELL_ABC: (10.0, 10.0, 30.0),
            K.OUTDIR: str(out),
            K.SCRIPT_PATH: None,
            K.GEN_POTCAR: False,
            K.FRAME_MODE: "index",
            K.FRAME_START: 0,
            K.FRAME_END: 10,
            K.FRAME_STEP: 1,
            K.TIME_START_FS: None,
            K.TIME_END_FS: None,
            K.TIME_STEP_FS: None,
        }
        cmd = _scripts.BaderBatchCmd("412", "Bader Batch")
        cmd.execute(ctx)

        assert stub.lookups == [
            ("md_analysis.workflows.scripts", "run_bader_batch"),
        ]
        call = stub.calls[0]
        # Batch frame slicing forwarded
        assert call["mode"] == "index"
        assert call["frame_start"] == 0
        assert call["frame_end"] == 10
        assert call["frame_step"] == 1
        # In index mode the time_* kwargs are explicitly None (not
        # missing) because the workflow signature lists them.
        assert call["time_start_fs"] is None
        assert call["time_end_fs"] is None
        assert call["time_step_fs"] is None


class TestTISingleCmd:
    """CLI 421 forwards colvar_id (TI workflow accepts it; the batch
    workflow does NOT, which is why CLI 422 was left on the legacy
    backend in Phase 4)."""

    def test_dispatches_to_run_ti_single_with_colvar_id(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        workdir = tmp_path / "out" / "ti_target_0.5"
        workdir.mkdir(parents=True)

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="ti_single",
                output_dir=tmp_path / "out",
                artifacts={"workdir": workdir},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_scripts, "lazy_import", stub)

        ctx = {
            K.INP_PATH: "sg.inp",
            K.XYZ: "sg.xyz",
            K.RESTART_PATH: "sg.restart",
            K.TARGET_AU: 0.5,
            K.OUTDIR: str(tmp_path / "out"),
            K.STEPS: 10000,
            K.COLVAR_ID: 2,
            K.WORKDIR_NAME: None,
            K.SCRIPT_PATH: None,
        }
        cmd = _scripts.TISingleCmd("421", "TI Single")
        cmd.execute(ctx)

        assert stub.lookups == [
            ("md_analysis.workflows.scripts", "run_ti_single"),
        ]
        call = stub.calls[0]
        assert call["target_au"] == 0.5
        assert call["steps"] == 10000
        assert call["colvar_id"] == 2
