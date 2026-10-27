"""Integration tests for center_slab_potential_analysis and related functions.

Uses ``data_example/potential/dense/`` as input data (continuous-mode
fixture; the legacy root-level layout was moved into ``dense/`` during
the utils/engines refactor).
Can also be run as a standalone script: python test/integration/potential/test_center_potential.py
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from md_analysis.electrochemical.potential import (
    center_slab_potential_analysis,
    fermi_energy_analysis,
    electrode_potential_analysis,
)

# Resolve data directory relative to this file
_DATA_DIR = (
    Path(__file__).resolve().parents[3]
    / "data_example" / "potential" / "dense"
)

# Skip all tests if data directory doesn't exist
pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"data_example/potential/dense/ not found at {_DATA_DIR}",
)


class TestCenterSlabPotentialAnalysis:

    def test_basic_cell_center(self, tmp_path: Path) -> None:
        """Run with center_mode='cell' (no xyz needed)."""
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            csv_path = center_slab_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="cell",
            )
        finally:
            os.chdir(old_cwd)

        assert csv_path.exists()
        data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        assert data.size > 0
        assert "step" in data.dtype.names
        assert "phi_center_ev" in data.dtype.names

    def test_interface_center(self, tmp_path: Path) -> None:
        """Run with center_mode='interface' using xyz trajectory."""
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            csv_path = center_slab_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="interface",
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
            )
        finally:
            os.chdir(old_cwd)

        assert csv_path.exists()
        # Should also produce slab_center_and_interfaces.csv
        slab_csv = tmp_path / "slab_center_and_interfaces.csv"
        assert slab_csv.exists()


class TestFermiEnergyAnalysis:

    def test_basic(self, tmp_path: Path) -> None:
        """End-to-end Fermi-energy analysis.

        Phase 5 Commit 3 pins CSV numeric equivalence with the typed
        ``engines.cp2k.read_fermi_series`` facade output (+ ``HA_TO_EV``):
        ``step`` / ``time_fs`` / ``fermi_raw`` / ``fermi_ev`` must agree
        byte-for-byte with the in-memory ``FermiRecord`` stream that
        ``read_fermi_series`` produces from the same ``md.out``.

        PNG output is checked for existence only (byte-equality is too
        brittle to enforce).
        """
        from md_analysis.engines.cp2k import read_fermi_series
        from md_analysis.utils.constants import HA_TO_EV

        csv_path = fermi_energy_analysis(
            _DATA_DIR / "md.out",
            output_dir=tmp_path,
            fermi_unit="au",
        )
        assert csv_path.exists()
        png_path = csv_path.with_suffix(".png")
        assert png_path.exists(), f"expected PNG at {png_path}"

        data = np.genfromtxt(
            csv_path, delimiter=",", names=True,
            dtype=None, encoding="utf-8",
        )
        assert data.size > 0
        for col in ("step", "time_fs", "fermi_raw", "fermi_ev"):
            assert col in data.dtype.names, f"missing column {col}"

        expected = read_fermi_series(_DATA_DIR / "md.out")
        assert len(expected) == data.size, (
            f"row count mismatch: CSV={data.size} vs facade={len(expected)}"
        )

        for i, rec in enumerate(expected):
            assert int(data["step"][i]) == rec.step
            # CP2K md.out fixture provides time_fs for every row; the
            # ``time_fs is None`` branch is exercised by unit tests
            # (FermiRecord.from_legacy_dict).
            csv_time = float(data["time_fs"][i])
            if rec.time_fs is None:
                assert np.isnan(csv_time)
            else:
                assert csv_time == pytest.approx(rec.time_fs, abs=1e-12)
            assert float(data["fermi_raw"][i]) == pytest.approx(
                rec.fermi_raw, abs=1e-12,
            )
            assert float(data["fermi_ev"][i]) == pytest.approx(
                rec.fermi_raw * HA_TO_EV, abs=1e-12,
            )


class TestElectrodePotentialAnalysis:

    def test_full_pipeline(self, tmp_path: Path) -> None:
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            u_csv = electrode_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                _DATA_DIR / "md.out",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="interface",
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
                fermi_unit="au",
            )
        finally:
            os.chdir(old_cwd)

        assert u_csv.exists()
        # Should also have center_potential.csv and fermi_energy.csv
        assert (tmp_path / "center_potential.csv").exists()
        assert (tmp_path / "fermi_energy.csv").exists()


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    _PREVIEW = Path(__file__).resolve().parents[2] / "_tmp_preview"
    _PREVIEW.mkdir(exist_ok=True)

    old_cwd = os.getcwd()
    try:
        os.chdir(_DATA_DIR)
        csv_path = electrode_potential_analysis(
            "md-POTENTIAL-v_hartree-1_*.cube",
            _DATA_DIR / "md.out",
            output_dir=_PREVIEW / "potential",
            thickness_ang=7.0,
            center_mode="interface",
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
        )
    finally:
        os.chdir(old_cwd)

    print(f"Electrode potential CSV: {csv_path}")
    sys.exit(0)
