"""Unit tests for ``md_analysis.workflows.composite.run_interface_analysis``.

Covers orchestration contract: layout, parameter forwarding,
artifact merging, collision detection, and metadata aggregation.

Leaf workflows (``run_water_three_panel`` / ``run_potential_full``)
are monkey-patched so this test does not depend on the
``data_example/potential/`` fixture (which is missing the
``md.inp`` file expected by the legacy integration test path —
this is a pre-existing baseline gap documented in
``context4agent/requirements/refactor_repair_plan.md`` and is
out of scope for the entrance refactor).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from md_analysis.workflows import (
    WorkflowResult,
    require_artifacts_exist,
    run_interface_analysis,
)


# ---------------------------------------------------------------------------
# Mock plumbing
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_leaf_workflows(monkeypatch, tmp_path):
    """Patch the two leaf workflows at their composite-side import sites.

    The composite imports ``run_water_three_panel`` / ``run_potential_full``
    via ``from .water import ...`` / ``from .potential import ...``, so
    the canonical patch site is
    ``md_analysis.workflows.composite.{run_water_three_panel,
    run_potential_full}``.

    Returns a spy dict the tests can introspect.
    """
    from md_analysis.workflows import composite as composite_mod

    state: dict[str, Any] = {
        "water_calls": [],
        "potential_calls": [],
    }

    def fake_water(*, xyz_path, md_inp_path, cell_abc, output_dir,
                   frame_start, frame_end, frame_step, verbose):
        state["water_calls"].append(
            {
                "xyz_path": Path(xyz_path),
                "md_inp_path": md_inp_path,
                "cell_abc": cell_abc,
                "output_dir": Path(output_dir),
                "frame_start": frame_start,
                "frame_end": frame_end,
                "frame_step": frame_step,
                "verbose": verbose,
            }
        )
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        csv_path = Path(output_dir) / "density.csv"
        png_path = Path(output_dir) / "water.png"
        csv_path.write_text("z,rho\n")
        png_path.write_bytes(b"\x89PNG")
        return WorkflowResult(
            name="water_three_panel",
            output_dir=Path(output_dir),
            artifacts={
                "density_csv": csv_path,
                "plot_png": png_path,
            },
            metadata={
                "n_water_artifacts_logical": 2,
            },
        )

    def fake_potential(**kwargs):
        state["potential_calls"].append(kwargs)
        outdir = Path(kwargs["output_dir"])
        outdir.mkdir(parents=True, exist_ok=True)
        elec = outdir / "electrode"
        elec.mkdir(parents=True, exist_ok=True)
        u_csv = elec / "electrode.csv"
        u_csv.write_text("t,U\n")
        return WorkflowResult(
            name="potential_full",
            output_dir=outdir,
            artifacts={"electrode_csv": u_csv},
            metadata={
                "input_mode": kwargs.get("input_mode", "continuous"),
                "has_fermi": True,
                "sub_analyses": ["electrode"],
            },
        )

    monkeypatch.setattr(composite_mod, "run_water_three_panel", fake_water)
    monkeypatch.setattr(composite_mod, "run_potential_full", fake_potential)
    return state


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


class TestInterfaceAnalysisLayout:
    def test_water_uses_canonical_subdir(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        out = tmp_path / "interface"
        run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            md_inp_path=tmp_path / "md.inp",
            output_dir=out,
            md_out_path=tmp_path / "md.out",
        )
        assert (
            mock_leaf_workflows["water_calls"][0]["output_dir"]
            == out / "water"
        )

    def test_potential_uses_canonical_subdir(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        out = tmp_path / "interface"
        run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=out,
            md_out_path=tmp_path / "md.out",
        )
        assert (
            mock_leaf_workflows["potential_calls"][0]["output_dir"]
            == out / "electrochemical" / "potential"
        )

    def test_output_dir_is_created(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        out = tmp_path / "fresh" / "nested" / "interface"
        result = run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=out,
        )
        assert out.is_dir()
        assert result.output_dir == out


# ---------------------------------------------------------------------------
# Parameter forwarding
# ---------------------------------------------------------------------------


class TestParameterForwarding:
    def test_shared_frame_slicing_reaches_both_leaves(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
            frame_start=10,
            frame_end=100,
            frame_step=2,
        )
        for which in ("water_calls", "potential_calls"):
            kwargs = mock_leaf_workflows[which][0]
            assert kwargs["frame_start"] == 10
            assert kwargs["frame_end"] == 100
            assert kwargs["frame_step"] == 2

    def test_potential_kwargs_forwarded(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
            md_out_path=tmp_path / "md.out",
            cube_pattern="custom-*.cube",
            thickness_ang=4.5,
            center_mode="cell",
            layer_tol_ang=0.8,
            fermi_unit="ev",
            compute_u=False,
            compute_phi_z=False,
            max_curves=5,
            thickness_end=12.0,
            metal_elements={"Cu", "Ag"},
        )
        pk = mock_leaf_workflows["potential_calls"][0]
        assert pk["cube_pattern"] == "custom-*.cube"
        assert pk["thickness_ang"] == 4.5
        assert pk["center_mode"] == "cell"
        assert pk["layer_tol_ang"] == 0.8
        assert pk["fermi_unit"] == "ev"
        assert pk["compute_u"] is False
        assert pk["compute_phi_z"] is False
        assert pk["max_curves"] == 5
        assert pk["thickness_end"] == 12.0
        assert pk["metal_elements"] == {"Cu", "Ag"}
        # xyz_path is shared with water so the potential leaf can do its
        # own interface detection in interface center_mode.
        assert pk["xyz_path"] == tmp_path / "md.xyz"

    def test_water_specific_args_do_not_leak_into_potential(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        """cell_abc and md_inp_path are water-only; the potential
        leaf must not receive them."""
        run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            md_inp_path=tmp_path / "md.inp",
            cell_abc=(10.0, 10.0, 30.0),
            output_dir=tmp_path / "out",
        )
        pk = mock_leaf_workflows["potential_calls"][0]
        assert "cell_abc" not in pk
        assert "md_inp_path" not in pk

    def test_xyz_path_normalised_to_path(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        run_interface_analysis(
            xyz_path=str(tmp_path / "md.xyz"),
            md_out_path=str(tmp_path / "md.out"),
            output_dir=str(tmp_path / "out"),
        )
        assert isinstance(
            mock_leaf_workflows["water_calls"][0]["xyz_path"], Path
        )
        assert isinstance(
            mock_leaf_workflows["potential_calls"][0]["md_out_path"], Path
        )


# ---------------------------------------------------------------------------
# Result assembly
# ---------------------------------------------------------------------------


class TestResultAssembly:
    def test_returns_workflow_result_with_merged_artifacts(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        result = run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "interface_analysis"
        # Water artifacts
        assert "density_csv" in result.artifacts
        assert "plot_png" in result.artifacts
        # Potential artifacts
        assert "electrode_csv" in result.artifacts
        require_artifacts_exist(result)

    def test_metadata_summarises_subworkflows(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        result = run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
            frame_start=0,
            frame_end=10,
        )
        meta = result.metadata
        assert meta["sub_workflows"] == ["water_three_panel", "potential_full"]
        assert meta["water_n_artifacts"] == 2
        assert meta["potential_n_artifacts"] == 1
        assert meta["n_artifacts"] == 3
        assert meta["water_output_dir"].endswith("water")
        assert meta["potential_output_dir"].endswith(
            "electrochemical/potential"
        )
        # Potential metadata is forwarded
        assert meta["potential_input_mode"] == "continuous"
        assert meta["potential_has_fermi"] is True
        assert meta["potential_sub_analyses"] == ["electrode"]
        # Shared slicing is recorded
        assert meta["frame_start"] == 0
        assert meta["frame_end"] == 10

    def test_to_dict_json_friendly(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        result = run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
        )
        json.dumps(result.to_dict())

    def test_no_extra_when_leaves_have_none(
        self, mock_leaf_workflows, tmp_path: Path
    ) -> None:
        """Water + potential leaf workflows do not populate ``extra``;
        the composite must not synthesise one either."""
        result = run_interface_analysis(
            xyz_path=tmp_path / "md.xyz",
            output_dir=tmp_path / "out",
        )
        assert result.extra is None


# ---------------------------------------------------------------------------
# Artifact collision detection
# ---------------------------------------------------------------------------


class TestArtifactCollision:
    def test_collision_raises_runtime_error(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        """If both leaves emit the same artifact key the composite
        refuses to silently overwrite."""
        from md_analysis.workflows import composite as composite_mod

        def water_with_collision(**kwargs):
            outdir = Path(kwargs["output_dir"])
            outdir.mkdir(parents=True, exist_ok=True)
            f = outdir / "x.csv"
            f.write_text("")
            return WorkflowResult(
                name="water_three_panel",
                output_dir=outdir,
                artifacts={"clash": f},
            )

        def potential_with_collision(**kwargs):
            outdir = Path(kwargs["output_dir"])
            outdir.mkdir(parents=True, exist_ok=True)
            f = outdir / "y.csv"
            f.write_text("")
            return WorkflowResult(
                name="potential_full",
                output_dir=outdir,
                artifacts={"clash": f},
                metadata={"sub_analyses": []},
            )

        monkeypatch.setattr(
            composite_mod, "run_water_three_panel", water_with_collision
        )
        monkeypatch.setattr(
            composite_mod, "run_potential_full", potential_with_collision
        )

        with pytest.raises(RuntimeError, match="clash"):
            run_interface_analysis(
                xyz_path=tmp_path / "md.xyz",
                output_dir=tmp_path / "out",
            )
