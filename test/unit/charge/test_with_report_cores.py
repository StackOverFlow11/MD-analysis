"""Tests for the *low-level* ``*_with_report`` wrappers in
``electrochemical.charge.Bader.AtomCharges``.

These cover the per-module wrappers that survive the entrance
refactor (see ``context4agent/requirements/refactor_repair_plan.md``)
and are independent of the ``md_analysis.main`` agent-facing
wrappers — which were removed in the Phase 3 charge legacy
cleanup.  They will continue to be exercised until
``workflows/charge.py`` fully subsumes them and the underlying
``*_with_report`` API is retired in a later phase.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.electrochemical.charge.Bader.AtomCharges import (
    CounterionChargeResult,
    TrackedChargeResult,
    counterion_charge_analysis_with_report,
    tracked_atom_charge_analysis_with_report,
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
