"""Batch 3: contract-backed charge tasks.

Covers:

- report dataclasses + ``run_*_with_report`` wrappers in ``main.py``
- ``*_with_report`` cores in ``AtomCharges.py``
- agent dispatch for ``charge_surface`` / ``charge_tracked`` /
  ``charge_counterion`` (schema, happy path, error classification)

The charge fixture reuses the single-frame POSCAR / ACF.dat / POTCAR
from ``data_example/bader/single_frame`` — we symlink/copy it into
two ``bader_t*_i*`` subdirs of a ``tmp_path``-based root so the sorted
frame discovery has something real to chew on.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.electrochemical.charge.Bader.AtomCharges import (
    CounterionChargeResult,
    TrackedChargeResult,
    counterion_charge_analysis_with_report,
    tracked_atom_charge_analysis_with_report,
)
from md_analysis.main import (
    ChargeSurfaceReport,
    CounterionChargeReport,
    TrackedChargeReport,
    run_charge_analysis_with_report,
    run_counterion_charge_analysis_with_report,
    run_tracked_charge_analysis_with_report,
)

DATA_DIR = (
    Path(__file__).resolve().parents[3]
    / "data_example" / "bader" / "single_frame"
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_frame_dir(frame_dir: Path) -> None:
    """Populate a frame dir with POSCAR (carrying an md_analysis::v1
    IndexMap comment), ACF.dat, and POTCAR.  The test-data POSCAR lacks
    the IndexMap comment, so we rewrite it through IndexMapper here —
    the mapping is effectively identity because we treat the existing
    POSCAR atom order as the XYZ order."""
    from ase.io import read as ase_read

    from md_analysis.scripts.utils.IndexMapper import (
        compute_index_map,
        write_poscar_with_map,
    )

    frame_dir.mkdir(parents=True, exist_ok=True)
    atoms = ase_read(DATA_DIR / "POSCAR", format="vasp")
    imap = compute_index_map(atoms, frame=int(frame_dir.name.rsplit("_i", 1)[-1]))
    write_poscar_with_map(atoms, frame_dir / "POSCAR", imap)
    shutil.copy2(DATA_DIR / "ACF.dat", frame_dir / "ACF.dat")
    shutil.copy2(DATA_DIR / "POTCAR", frame_dir / "POTCAR")


@pytest.fixture
def bader_root(tmp_path):
    """Build a minimal Bader frame layout: two bader_t*_i* subdirs."""
    root = tmp_path / "bader_root"
    for t, i in [(100, 1), (200, 2)]:
        _make_frame_dir(root / f"bader_t{t}_i{i}")
    return {"root": root, "outdir": tmp_path / "out", "tmp": tmp_path}


# ---------------------------------------------------------------------------
# Wrapper-level — tracked / counterion *_with_report
# ---------------------------------------------------------------------------


class TestTrackedAnalysisWithReport:
    def test_success_returns_result(self, bader_root):
        out = bader_root["outdir"] / "tracked"
        out.mkdir(parents=True)
        r = tracked_atom_charge_analysis_with_report(
            bader_root["root"],
            atom_indices_xyz=[0, 1],
            output_dir=out,
        )
        assert isinstance(r, TrackedChargeResult)
        assert r.csv_path.is_file()
        assert r.png_path.is_file()
        assert r.n_frames == 2
        assert r.n_atoms_tracked == 2
        assert r.atom_indices_xyz == (0, 1)


class TestCounterionAnalysisWithReport:
    def test_success_zero_detections_reports_png_as_none(self, bader_root):
        """Stock fixture (pure Cu/Ag/O/H, no ions) → zero counterions.

        Invariant: :func:`plot_counterion_charges` does not write a PNG
        when nothing was detected, so the result's ``png_path`` must
        be ``None`` — never a bogus path to a missing file.
        """
        out = bader_root["outdir"] / "ci"
        out.mkdir(parents=True)
        r = counterion_charge_analysis_with_report(
            bader_root["root"], output_dir=out,
        )
        assert isinstance(r, CounterionChargeResult)
        assert r.csv_path.is_file()
        assert r.summary_path.is_file()
        assert r.n_frames == 2
        assert r.n_unique_counterions == 0
        assert r.png_path is None


# ---------------------------------------------------------------------------
# Wrapper-level — main.py run_*_with_report
# ---------------------------------------------------------------------------


class TestRunChargeAnalysisWithReport:
    def test_success(self, bader_root):
        r = run_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
            method="counterion",
        )
        assert isinstance(r, ChargeSurfaceReport)
        assert r.charge_csv.is_file()
        assert r.charge_png.is_file()
        assert r.n_frames == 2
        assert r.method == "counterion"
        assert isinstance(r.sigma_aligned_mean, float)

    def test_report_to_dict_is_json_serializable(self, bader_root):
        r = run_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
            method="counterion",
        )
        json.dumps(r.to_dict())


class TestRunTrackedChargeAnalysisWithReport:
    def test_success(self, bader_root):
        r = run_tracked_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
            atom_indices_xyz=[0, 1, 2],
        )
        assert isinstance(r, TrackedChargeReport)
        assert r.tracked_charge_csv.is_file()
        assert r.tracked_charge_png.is_file()
        assert r.n_frames == 2
        assert r.n_atoms_tracked == 3
        assert list(r.atom_indices_xyz) == [0, 1, 2]

    def test_report_to_dict_is_json_serializable(self, bader_root):
        r = run_tracked_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
            atom_indices_xyz=[0, 1],
        )
        json.dumps(r.to_dict())


class TestRunCounterionChargeAnalysisWithReport:
    def test_success_zero_detections_reports_png_as_none(self, bader_root):
        r = run_counterion_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
        )
        assert isinstance(r, CounterionChargeReport)
        assert r.counterion_charge_csv.is_file()
        assert r.counterion_summary_csv.is_file()
        assert r.n_frames == 2
        assert r.n_unique_counterions == 0
        # Propagated nullability from CounterionChargeResult — when the
        # underlying detection was empty the report must explicitly
        # mark the PNG as absent (None) rather than carry a bogus Path.
        assert r.counterion_charge_png is None

    def test_report_to_dict_is_json_serializable(self, bader_root):
        r = run_counterion_charge_analysis_with_report(
            output_dir=bader_root["outdir"],
            root_dir=bader_root["root"],
        )
        payload = r.to_dict()
        # Serialises without custom encoder
        json.dumps(payload)
        # None is preserved through to_dict()
        assert payload["counterion_charge_png"] is None


# ---------------------------------------------------------------------------
# charge_surface schema + dispatch
# ---------------------------------------------------------------------------


class TestChargeSurfaceSchema:
    def test_schema_keys(self):
        s = get_task_schema("charge_surface")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("charge_surface")
        assert s["parameters"]["required"] == ["output_dir"]

    def test_method_enum(self):
        s = get_task_schema("charge_surface")
        assert set(
            s["parameters"]["properties"]["method"].get("enum", [])
        ) == {"counterion", "layer"}

    def test_normal_enum(self):
        s = get_task_schema("charge_surface")
        assert set(
            s["parameters"]["properties"]["normal"].get("enum", [])
        ) == {"a", "b", "c"}

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["charge_surface"].contract is not None


class TestChargeSurfaceDispatch:
    def test_success(self, bader_root):
        r = dispatch("charge_surface", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "method": "counterion",
        })
        assert r.success, r.errors
        assert "charge_csv" in r.outputs
        assert "charge_png" in r.outputs
        for k in (
            "n_frames", "method",
            "sigma_aligned_mean", "sigma_aligned_std",
            "sigma_opposed_mean", "sigma_opposed_std",
            "phi_cumavg_last", "phi_reference",
        ):
            assert k in r.summary
        # summary is JSON-serializable
        json.dumps(r.summary)

    def test_missing_root_dir(self, tmp_path):
        r = dispatch("charge_surface", {
            "output_dir": str(tmp_path / "out"),
            "root_dir": str(tmp_path / "nope"),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_missing_poscar_in_frame(self, bader_root):
        # Corrupt a frame by deleting its POSCAR
        frame = bader_root["root"] / "bader_t100_i1"
        (frame / "POSCAR").unlink()
        r = dispatch("charge_surface", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_method(self, bader_root):
        r = dispatch("charge_surface", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "method": "bogus",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_invalid_normal(self, bader_root):
        r = dispatch("charge_surface", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "normal": "z",
        })
        assert not r.success
        assert r.error_type == "validation"


# ---------------------------------------------------------------------------
# charge_tracked schema + dispatch
# ---------------------------------------------------------------------------


class TestChargeTrackedSchema:
    def test_schema_keys(self):
        s = get_task_schema("charge_tracked")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("charge_tracked")
        assert set(s["parameters"]["required"]) == {
            "output_dir", "atom_indices_xyz",
        }

    def test_atom_indices_schema(self):
        s = get_task_schema("charge_tracked")
        ai = s["parameters"]["properties"]["atom_indices_xyz"]
        assert ai["type"] == "array"
        assert ai["items"]["type"] == "integer"
        assert ai["items"]["minimum"] == 0
        assert ai["minItems"] == 1

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["charge_tracked"].contract is not None


class TestChargeTrackedDispatch:
    def test_success(self, bader_root):
        r = dispatch("charge_tracked", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "atom_indices_xyz": [0, 1],
        })
        assert r.success, r.errors
        assert "tracked_charge_csv" in r.outputs
        assert "tracked_charge_png" in r.outputs
        assert r.summary["n_frames"] == 2
        assert r.summary["n_atoms_tracked"] == 2
        assert r.summary["atom_indices_xyz"] == [0, 1]
        json.dumps(r.summary)

    def test_missing_root_dir(self, tmp_path):
        r = dispatch("charge_tracked", {
            "output_dir": str(tmp_path / "out"),
            "root_dir": str(tmp_path / "nope"),
            "atom_indices_xyz": [0],
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_missing_poscar_in_frame(self, bader_root):
        frame = bader_root["root"] / "bader_t100_i1"
        (frame / "POSCAR").unlink()
        r = dispatch("charge_tracked", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "atom_indices_xyz": [0],
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_empty_atom_list(self, bader_root):
        r = dispatch("charge_tracked", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "atom_indices_xyz": [],
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_negative_index(self, bader_root):
        r = dispatch("charge_tracked", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "atom_indices_xyz": [0, -1],
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_out_of_bounds_index(self, bader_root):
        r = dispatch("charge_tracked", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "atom_indices_xyz": [99999],
        })
        assert not r.success
        assert r.error_type == "validation"


# ---------------------------------------------------------------------------
# charge_counterion schema + dispatch
# ---------------------------------------------------------------------------


class TestChargeCounterionSchema:
    def test_schema_keys(self):
        s = get_task_schema("charge_counterion")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required(self):
        s = get_task_schema("charge_counterion")
        assert s["parameters"]["required"] == ["output_dir"]

    def test_normal_enum(self):
        s = get_task_schema("charge_counterion")
        assert set(
            s["parameters"]["properties"]["normal"].get("enum", [])
        ) == {"a", "b", "c"}

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["charge_counterion"].contract is not None


class TestChargeCounterionDispatch:
    def test_success_zero_detections_omits_png_from_outputs(self, bader_root):
        """End-to-end invariant: when no counterions are detected, the
        handler MUST NOT surface ``counterion_charge_png`` in
        ``TaskResult.outputs`` — the file doesn't exist on disk and
        ``outputs`` is reserved for real artifacts."""
        r = dispatch("charge_counterion", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
        })
        assert r.success, r.errors
        # CSV artifacts always produced
        assert "counterion_charge_csv" in r.outputs
        assert "counterion_summary_csv" in r.outputs
        # PNG not produced for empty detection → not in outputs
        assert "counterion_charge_png" not in r.outputs, (
            "TaskResult.outputs must not advertise a non-existent PNG; "
            "the contract treats counterion_charge_png as nullable."
        )
        # Every declared output path actually exists on disk
        for k, path in r.outputs.items():
            assert Path(path).is_file(), (
                f"outputs[{k!r}] points to a non-existent file: {path}"
            )
        # Summary always populated
        assert r.summary["n_frames"] == 2
        assert r.summary["n_unique_counterions"] == 0
        json.dumps(r.summary)

    def test_missing_root_dir(self, tmp_path):
        r = dispatch("charge_counterion", {
            "output_dir": str(tmp_path / "out"),
            "root_dir": str(tmp_path / "nope"),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_missing_poscar_in_frame(self, bader_root):
        frame = bader_root["root"] / "bader_t100_i1"
        (frame / "POSCAR").unlink()
        r = dispatch("charge_counterion", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_normal(self, bader_root):
        r = dispatch("charge_counterion", {
            "output_dir": str(bader_root["outdir"]),
            "root_dir": str(bader_root["root"]),
            "normal": "z",
        })
        assert not r.success
        assert r.error_type == "validation"
