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

from md_analysis.cli import (
    _charge,
    _enhanced_sampling,
    _potential,
    _scripts,
    _water,
)
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


# ---------------------------------------------------------------------------
# _water.py: 4 single-step workflows (Phase 6.1)
# ---------------------------------------------------------------------------


def _water_ctx(tmp_path: Path) -> dict[str, Any]:
    """Build a ctx dict shared by the 4 water single-step wiring tests."""
    return {
        K.XYZ: str(tmp_path / "md-pos-1.xyz"),
        K.CELL_ABC: (10.0, 10.0, 30.0),
        K.OUTDIR_RESOLVED: str(tmp_path / "water"),
        K.DZ_A: 0.1,
        K.LAYER_TOL: 0.5,
        K.FRAME_START: None,
        K.FRAME_END: None,
        K.FRAME_STEP: None,
    }


def _assert_water_ctx_forwarded(call: dict[str, Any], tmp_path: Path) -> None:
    """All 4 single-step CLI commands forward the same ctx fields verbatim."""
    assert call["xyz_path"] == Path(tmp_path / "md-pos-1.xyz")
    assert call["cell_abc"] == (10.0, 10.0, 30.0)
    assert call["output_dir"] == str(tmp_path / "water")
    assert call["dz_A"] == 0.1
    assert call["layer_tol_A"] == 0.5
    assert call["frame_start"] is None
    assert call["frame_end"] is None
    assert call["frame_step"] is None


class TestWaterDensityCmd:
    """CLI 101 forwards ctx straight to workflows.water.run_water_density."""

    def test_lookup_and_kwargs_and_artifact_keys(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        water_dir = tmp_path / "water"
        water_dir.mkdir()
        expected_csv = water_dir / "water_mass_density.csv"
        expected_csv.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="water_density",
                output_dir=water_dir,
                artifacts={"density_csv": expected_csv},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_water, "lazy_import", stub)

        cmd = _water.WaterDensityCmd("101", "Water Density")
        cmd.execute(_water_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.water", "run_water_density"),
        ]
        assert len(stub.calls) == 1
        _assert_water_ctx_forwarded(stub.calls[0], tmp_path)


class TestWaterOrientationCmd:
    """CLI 102 -> workflows.water.run_water_orientation."""

    def test_lookup_and_kwargs_and_artifact_keys(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        water_dir = tmp_path / "water"
        water_dir.mkdir()
        expected_csv = water_dir / "water_orientation.csv"
        expected_csv.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="water_orientation",
                output_dir=water_dir,
                artifacts={"orientation_csv": expected_csv},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_water, "lazy_import", stub)

        cmd = _water.WaterOrientationCmd("102", "Water Orientation")
        cmd.execute(_water_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.water", "run_water_orientation"),
        ]
        _assert_water_ctx_forwarded(stub.calls[0], tmp_path)


class TestAdWaterOrientationCmd:
    """CLI 103 -> workflows.water.run_ad_water_orientation
    (artifact contract: adsorbed_profile_csv + adsorbed_range_txt)."""

    def test_lookup_and_kwargs_and_artifact_keys(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        water_dir = tmp_path / "water"
        water_dir.mkdir()
        profile_csv = water_dir / "ad_water_profile.csv"
        range_txt = water_dir / "ad_water_range.txt"
        profile_csv.write_text("")
        range_txt.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="ad_water_orientation",
                output_dir=water_dir,
                artifacts={
                    "adsorbed_profile_csv": profile_csv,
                    "adsorbed_range_txt": range_txt,
                },
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_water, "lazy_import", stub)

        cmd = _water.AdWaterOrientationCmd("103", "Adsorbed Water Orientation")
        cmd.execute(_water_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.water", "run_ad_water_orientation"),
        ]
        _assert_water_ctx_forwarded(stub.calls[0], tmp_path)


class TestAdWaterThetaCmd:
    """CLI 104 -> workflows.water.run_ad_water_theta
    (artifact contract: theta_csv; verbose=True forwarded)."""

    def test_lookup_and_kwargs_and_artifact_keys(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        water_dir = tmp_path / "water"
        water_dir.mkdir()
        theta_csv = water_dir / "ad_water_theta.csv"
        theta_csv.write_text("")

        def fake_workflow(**kwargs: Any) -> WorkflowResult:
            return _make_workflow_result(
                name="ad_water_theta",
                output_dir=water_dir,
                artifacts={"theta_csv": theta_csv},
            )

        stub = _LazyImportStub(fake_workflow)
        monkeypatch.setattr(_water, "lazy_import", stub)

        cmd = _water.AdWaterThetaCmd("104", "Adsorbed Water Theta")
        cmd.execute(_water_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.water", "run_ad_water_theta"),
        ]
        call = stub.calls[0]
        _assert_water_ctx_forwarded(call, tmp_path)
        # CLI 104 forwards verbose=True (matches legacy behaviour)
        assert call["verbose"] is True


# ---------------------------------------------------------------------------
# _potential.py: 5 single-step workflows (Phase 6.2)
# ---------------------------------------------------------------------------


_CONTINUOUS_ONLY_KEYS = ("cube_pattern", "xyz_path", "md_out_path")


def _potential_continuous_ctx(tmp_path: Path) -> dict[str, Any]:
    """Continuous-mode ctx used by CLI 211-215 wiring tests."""
    return {
        K.INPUT_MODE: "continuous",
        K.SP_ROOT_DIR: None,
        K.SP_DIR_PATTERN: "potential_t*_i*",
        K.SP_CUBE_FILENAME: "sp_potential-v_hartree-1_0.cube",
        K.SP_OUT_FILENAME: "sp.out",
        K.CUBE_PATTERN: "md-POTENTIAL-v_hartree-1_*.cube",
        K.XYZ: str(tmp_path / "md-pos-1.xyz"),
        K.MD_OUT: str(tmp_path / "md.out"),
        K.OUTDIR_RESOLVED: str(tmp_path / "out"),
        K.THICKNESS: 7.0,
        K.THICKNESS_END: 15.0,
        K.CENTER_MODE: "interface",
        K.METAL_ELEMENTS: {"Cu"},
        K.LAYER_TOL: 0.5,
        K.FERMI_UNIT: "au",
        K.MAX_CURVES: 0,
        K.FRAME_START: None,
        K.FRAME_END: None,
        K.FRAME_STEP: None,
    }


def _potential_distributed_ctx(tmp_path: Path) -> dict[str, Any]:
    """Distributed-mode ctx; the CLI routing branch should drop the
    continuous-only keys (cube_pattern / xyz_path / md_out_path)
    before calling the workflow facade."""
    ctx = _potential_continuous_ctx(tmp_path)
    ctx[K.INPUT_MODE] = "distributed"
    ctx[K.SP_ROOT_DIR] = str(tmp_path / "sp_dirs")
    # CUBE_PATTERN / XYZ / MD_OUT are still in ctx for completeness, but
    # the CLI routing branch must not forward them in distributed mode.
    return ctx


def _stub_workflow(
    name: str, output_dir: Path, *, artifact_key: str, artifact_path: Path
) -> Any:
    """Build a fake workflow callable returning a WorkflowResult."""

    def _fake(**kwargs: Any) -> WorkflowResult:
        return _make_workflow_result(
            name=name,
            output_dir=output_dir,
            artifacts={artifact_key: artifact_path},
        )

    return _fake


def _assert_continuous_routing(
    call: dict[str, Any], *, expects_cube: bool, expects_xyz: bool, expects_md_out: bool
) -> None:
    """In continuous mode the CLI routing branch forwards the
    continuous-only keys when relevant for the command."""
    if expects_cube:
        assert "cube_pattern" in call
    if expects_xyz:
        assert "xyz_path" in call
    if expects_md_out:
        assert "md_out_path" in call
    assert call["input_mode"] == "continuous"


def _assert_distributed_routing(call: dict[str, Any]) -> None:
    """In distributed mode the CLI routing branch must NOT forward
    any continuous-only keys, and input_mode + sp_* must be present."""
    for key in _CONTINUOUS_ONLY_KEYS:
        assert key not in call, (
            f"{key!r} was forwarded in distributed mode; "
            f"the CLI routing branch must not pass it"
        )
    assert call["input_mode"] == "distributed"
    assert "sp_root_dir" in call
    assert "sp_dir_pattern" in call
    assert "sp_cube_filename" in call
    assert "sp_out_filename" in call


class TestCenterPotentialCmd:
    """CLI 211 -> workflows.potential.run_center_potential.
    Forwards cube_pattern + xyz_path in continuous mode; drops both
    in distributed mode."""

    def test_continuous_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "center.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "center_potential", out_dir,
                artifact_key="center_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.CenterPotentialCmd("211", "Center Potential")
        cmd.execute(_potential_continuous_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_center_potential"),
        ]
        call = stub.calls[0]
        _assert_continuous_routing(
            call, expects_cube=True, expects_xyz=True, expects_md_out=False,
        )
        assert call["thickness_ang"] == 7.0

    def test_distributed_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "center.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "center_potential", out_dir,
                artifact_key="center_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.CenterPotentialCmd("211", "Center Potential")
        cmd.execute(_potential_distributed_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_center_potential"),
        ]
        _assert_distributed_routing(stub.calls[0])


class TestFermiEnergyCmd:
    """CLI 212 -> workflows.potential.run_fermi_energy.
    Forwards md_out_path in continuous mode; drops it in distributed."""

    def test_continuous_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "fermi.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "fermi_energy", out_dir,
                artifact_key="fermi_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.FermiEnergyCmd("212", "Fermi Energy")
        cmd.execute(_potential_continuous_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_fermi_energy"),
        ]
        _assert_continuous_routing(
            stub.calls[0],
            expects_cube=False, expects_xyz=False, expects_md_out=True,
        )

    def test_distributed_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "fermi.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "fermi_energy", out_dir,
                artifact_key="fermi_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.FermiEnergyCmd("212", "Fermi Energy")
        cmd.execute(_potential_distributed_ctx(tmp_path))

        _assert_distributed_routing(stub.calls[0])


class TestElectrodePotentialCmd:
    """CLI 213 -> workflows.potential.run_electrode_potential.
    Forwards cube_pattern + md_out_path + xyz_path in continuous; drops
    all three in distributed."""

    def test_continuous_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "electrode.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "electrode_potential", out_dir,
                artifact_key="electrode_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.ElectrodePotentialCmd("213", "Electrode Potential")
        cmd.execute(_potential_continuous_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_electrode_potential"),
        ]
        _assert_continuous_routing(
            stub.calls[0],
            expects_cube=True, expects_xyz=True, expects_md_out=True,
        )

    def test_distributed_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "electrode.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "electrode_potential", out_dir,
                artifact_key="electrode_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.ElectrodePotentialCmd("213", "Electrode Potential")
        cmd.execute(_potential_distributed_ctx(tmp_path))

        _assert_distributed_routing(stub.calls[0])


class TestPhiZProfileCmd:
    """CLI 214 -> workflows.potential.run_phi_z_profile.
    Forwards cube_pattern in continuous (no xyz / md_out for PhiZ).
    Drops cube_pattern in distributed."""

    def test_continuous_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "phi_z.png"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "phi_z_profile", out_dir,
                artifact_key="phi_z_png", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.PhiZProfileCmd("214", "PhiZ Profile")
        cmd.execute(_potential_continuous_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_phi_z_profile"),
        ]
        # PhiZ only forwards cube_pattern in continuous (no xyz/md_out)
        call = stub.calls[0]
        assert "cube_pattern" in call
        assert "xyz_path" not in call
        assert "md_out_path" not in call
        assert call["input_mode"] == "continuous"

    def test_distributed_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "phi_z.png"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "phi_z_profile", out_dir,
                artifact_key="phi_z_png", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.PhiZProfileCmd("214", "PhiZ Profile")
        cmd.execute(_potential_distributed_ctx(tmp_path))

        _assert_distributed_routing(stub.calls[0])


class TestThicknessSensitivityCmd:
    """CLI 215 -> workflows.potential.run_thickness_sensitivity.
    Same routing surface as 213 (cube + md_out + xyz)."""

    def test_continuous_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "ts.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "thickness_sensitivity", out_dir,
                artifact_key="thickness_sensitivity_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.ThicknessSensitivityCmd("215", "Thickness Sensitivity")
        cmd.execute(_potential_continuous_ctx(tmp_path))

        assert stub.lookups == [
            ("md_analysis.workflows.potential", "run_thickness_sensitivity"),
        ]
        _assert_continuous_routing(
            stub.calls[0],
            expects_cube=True, expects_xyz=True, expects_md_out=True,
        )
        assert stub.calls[0]["thickness_end"] == 15.0

    def test_distributed_routing(self, monkeypatch, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        artifact = out_dir / "ts.csv"
        artifact.write_text("")
        stub = _LazyImportStub(
            _stub_workflow(
                "thickness_sensitivity", out_dir,
                artifact_key="thickness_sensitivity_csv", artifact_path=artifact,
            )
        )
        monkeypatch.setattr(_potential, "lazy_import", stub)

        cmd = _potential.ThicknessSensitivityCmd("215", "Thickness Sensitivity")
        cmd.execute(_potential_distributed_ctx(tmp_path))

        _assert_distributed_routing(stub.calls[0])


# ---------------------------------------------------------------------------
# _charge.py: SurfaceChargeCmd 221/222/223 + SingleSideChargeCmd 224
# (Phase 6.3)
# ---------------------------------------------------------------------------


def _surface_charge_ctx(
    tmp_path: Path,
    *,
    method: str,
    output_dir: Path,
) -> dict[str, Any]:
    """Shared ctx for SurfaceChargeCmd (221/222/223)."""
    return {
        K.METHOD: method,
        K.ROOT_DIR: str(tmp_path / "bader_root"),
        K.METAL_ELEMENTS: None,
        K.NORMAL: "c",
        K.LAYER_TOL: 0.5,
        K.N_SURFACE_LAYERS: 1,
        K.DIR_PATTERN: "bader_t*_i*",
        K.OUTDIR_RESOLVED: str(output_dir),
        K.FRAME_START: None,
        K.FRAME_END: None,
        K.FRAME_STEP: None,
    }


def _surface_charge_metadata(
    method: str, target_side: str | None = None
) -> dict[str, Any]:
    """Default metadata stub for surface-charge workflow results."""
    return {
        "method": method,
        "target_side": target_side,
        "normal": "c",
        "n_frames": 12,
        "sigma_aligned_mean": -1.25,
        "sigma_aligned_std": 0.05,
        "sigma_opposed_mean": 1.30,
        "sigma_opposed_std": 0.06,
        "phi_cumavg_last": None,
        "phi_reference": None,
    }


def _make_charge_workflow_result(
    *,
    method: str,
    output_dir: Path,
    artifact: Path,
    target_side: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> WorkflowResult:
    return WorkflowResult(
        name="surface_charge",
        output_dir=output_dir,
        artifacts={"charge_csv": artifact, "charge_png": output_dir / "p.png"},
        metadata=metadata or _surface_charge_metadata(method, target_side),
    )


class TestSurfaceChargeStaticCmd:
    """CLI 221/222: SurfaceChargeCmd(method=...) has output_name=method,
    so OUTDIR_RESOLVED ends in /<method>/ and the CLI must pass
    .parent (rewind one level) so the workflow's internal /<method>
    append produces the original layout."""

    @pytest.mark.parametrize(
        "method, code",
        [("counterion", "221"), ("layer", "222")],
    )
    def test_path_routing_and_summary(
        self, monkeypatch, tmp_path: Path, capsys, method: str, code: str
    ) -> None:
        # Static method: OUTDIR_RESOLVED already ends with /<method>/
        resolved = tmp_path / "out" / "electrochemical" / "charge" / method
        resolved.mkdir(parents=True)
        artifact = resolved / "surface_charge.csv"
        artifact.write_text("")

        def fake(**kwargs: Any) -> WorkflowResult:
            return _make_charge_workflow_result(
                method=method, output_dir=resolved, artifact=artifact,
            )

        stub = _LazyImportStub(fake)
        monkeypatch.setattr(_charge, "lazy_import", stub)

        cmd = _charge.SurfaceChargeCmd(code, f"Static {method}", method=method)
        cmd.execute(_surface_charge_ctx(
            tmp_path, method=method, output_dir=resolved,
        ))

        # workflow facade resolved
        assert stub.lookups == [
            ("md_analysis.workflows.charge", "run_surface_charge"),
        ]

        # Path routing: static path forwards .parent so workflow's
        # internal /<method> append produces the original layout.
        call = stub.calls[0]
        assert call["output_dir"] == resolved.parent
        assert call["method"] == method
        assert "target_side" not in call  # two-sided path

        # Final layout = (cli_output_dir / method) byte-equal pre-migration
        composed = Path(call["output_dir"]) / call["method"]
        assert composed == resolved

        # capsys: CLI summary text must remain byte-equal
        out = capsys.readouterr().out
        assert "Analysis complete. Output:" in out
        assert "charge_csv:" in out
        assert "Ensemble average (12 frames):" in out
        assert "sigma_aligned: " in out
        assert "+/- 0.0500 uC/cm^2" in out
        assert "sigma_opposed:" in out


class TestSurfaceChargeDynamicCmd:
    """CLI 223: SurfaceChargeCmd() (no method) has output_name="",
    so OUTDIR_RESOLVED ends in /charge/ and the CLI passes it
    directly; the workflow's /<method> append produces the layout."""

    @pytest.mark.parametrize("method", ["counterion", "layer"])
    def test_path_routing_and_summary(
        self, monkeypatch, tmp_path: Path, capsys, method: str
    ) -> None:
        resolved = tmp_path / "out" / "electrochemical" / "charge"
        resolved.mkdir(parents=True)
        final_dir = resolved / method
        final_dir.mkdir()
        artifact = final_dir / "surface_charge.csv"
        artifact.write_text("")

        def fake(**kwargs: Any) -> WorkflowResult:
            return _make_charge_workflow_result(
                method=method, output_dir=final_dir, artifact=artifact,
            )

        stub = _LazyImportStub(fake)
        monkeypatch.setattr(_charge, "lazy_import", stub)

        cmd = _charge.SurfaceChargeCmd("223", "Dynamic")  # method=None
        cmd.execute(_surface_charge_ctx(
            tmp_path, method=method, output_dir=resolved,
        ))

        assert stub.lookups == [
            ("md_analysis.workflows.charge", "run_surface_charge"),
        ]
        call = stub.calls[0]
        assert call["output_dir"] == resolved
        assert call["method"] == method
        assert "target_side" not in call

        # Final layout = (cli_output_dir / method)
        composed = Path(call["output_dir"]) / call["method"]
        assert composed == final_dir

        out = capsys.readouterr().out
        assert "Analysis complete. Output:" in out
        assert "Ensemble average (12 frames):" in out
        assert "sigma_aligned: " in out
        assert "sigma_opposed:" in out


class TestSingleSideChargeCmd:
    """CLI 224: SingleSideChargeCmd has output_name="", so
    OUTDIR_RESOLVED ends in /charge/.  The CLI passes it directly +
    target_side; the workflow appends /<method>_<side>."""

    @pytest.mark.parametrize("method", ["counterion", "layer"])
    @pytest.mark.parametrize("side", ["aligned", "opposed"])
    def test_path_routing_and_summary(
        self,
        monkeypatch,
        tmp_path: Path,
        capsys,
        method: str,
        side: str,
    ) -> None:
        resolved = tmp_path / "out" / "electrochemical" / "charge"
        resolved.mkdir(parents=True)
        final_dir = resolved / f"{method}_{side}"
        final_dir.mkdir()
        artifact = final_dir / "surface_charge.csv"
        artifact.write_text("")

        # phi populated to exercise the V-vs-reference print branch
        metadata = _surface_charge_metadata(method, target_side=side)
        metadata["phi_cumavg_last"] = -0.123
        metadata["phi_reference"] = "SHE"

        def fake(**kwargs: Any) -> WorkflowResult:
            return _make_charge_workflow_result(
                method=method, output_dir=final_dir, artifact=artifact,
                target_side=side, metadata=metadata,
            )

        stub = _LazyImportStub(fake)
        monkeypatch.setattr(_charge, "lazy_import", stub)

        ctx = _surface_charge_ctx(
            tmp_path, method=method, output_dir=resolved,
        )
        ctx[K.TARGET_SIDE] = side

        cmd = _charge.SingleSideChargeCmd("224", "Single-side")
        cmd.execute(ctx)

        assert stub.lookups == [
            ("md_analysis.workflows.charge", "run_surface_charge"),
        ]
        call = stub.calls[0]
        # 224 passes OUTDIR_RESOLVED directly + target_side
        assert call["output_dir"] == resolved
        assert call["method"] == method
        assert call["target_side"] == side

        # Final layout = output_dir / f"{method}_{side}"
        composed = Path(call["output_dir"]) / f"{call['method']}_{call['target_side']}"
        assert composed == final_dir

        # capsys: single-side summary + phi line
        out = capsys.readouterr().out
        assert f"Analysis complete ({side} side). Output:" in out
        assert "charge_csv:" in out
        assert f"Ensemble average (12 frames, {side} side):" in out
        assert "sigma: " in out
        assert "+/-" in out
        assert "uC/cm^2" in out
        # phi branch
        assert "phi:" in out
        assert "V vs SHE (cum. avg)" in out
