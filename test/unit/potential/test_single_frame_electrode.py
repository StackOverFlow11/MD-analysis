"""Regression test: electrode_potential_analysis merge logic with single-frame CSVs.

np.genfromtxt(..., names=True) returns a 0-d structured array for a single data
row.  The np.atleast_1d() guard added in the P0 fix ensures the downstream
iteration does not crash.
"""

from __future__ import annotations

import numpy as np
import pytest

from md_analysis.utils.config import DP_A_H3O_W_EV, MU_HPLUS_G0_EV, DELTA_E_ZP_EV


def _merge_csvs(center_csv_path, fermi_csv_path):
    """Reproduce the merge logic from electrode_potential_analysis (lines 522-538)."""
    center_data = np.atleast_1d(
        np.genfromtxt(center_csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )
    fermi_data = np.atleast_1d(
        np.genfromtxt(fermi_csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )

    fermi_by_step = {int(r["step"]): float(r["fermi_ev"]) for r in fermi_data}

    u_rows: list[dict] = []
    for row in center_data:
        s = int(row["step"])
        phi_ev = float(row["phi_center_ev"])
        if s not in fermi_by_step:
            continue
        ef_ev = fermi_by_step[s]
        u_v = -ef_ev + phi_ev + DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV
        u_rows.append({"step": s, "U_vs_SHE_V": float(u_v)})

    return u_rows


class TestSingleFrameElectrode:
    """Verify that single-row CSVs do not crash the merge logic."""

    def test_single_row_merge(self, tmp_path):
        center_csv = tmp_path / "center_potential.csv"
        center_csv.write_text("step,phi_center_ev,phi_center_cumavg_ev\n100,-2.5,-2.5\n")

        fermi_csv = tmp_path / "fermi_energy.csv"
        fermi_csv.write_text("step,time_fs,fermi_raw,fermi_ev,fermi_cumavg_ev\n100,50.0,-0.18,-4.9,-4.9\n")

        rows = _merge_csvs(center_csv, fermi_csv)

        assert len(rows) == 1
        assert rows[0]["step"] == 100
        assert isinstance(rows[0]["U_vs_SHE_V"], float)

    def test_no_matching_steps(self, tmp_path):
        center_csv = tmp_path / "center_potential.csv"
        center_csv.write_text("step,phi_center_ev,phi_center_cumavg_ev\n100,-2.5,-2.5\n")

        fermi_csv = tmp_path / "fermi_energy.csv"
        fermi_csv.write_text("step,time_fs,fermi_raw,fermi_ev,fermi_cumavg_ev\n200,100.0,-0.18,-4.9,-4.9\n")

        rows = _merge_csvs(center_csv, fermi_csv)
        assert len(rows) == 0

    def test_multi_row_still_works(self, tmp_path):
        center_csv = tmp_path / "center_potential.csv"
        center_csv.write_text(
            "step,phi_center_ev,phi_center_cumavg_ev\n"
            "100,-2.5,-2.5\n"
            "200,-2.6,-2.55\n"
        )

        fermi_csv = tmp_path / "fermi_energy.csv"
        fermi_csv.write_text(
            "step,time_fs,fermi_raw,fermi_ev,fermi_cumavg_ev\n"
            "100,50.0,-0.18,-4.9,-4.9\n"
            "200,100.0,-0.19,-5.0,-4.95\n"
        )

        rows = _merge_csvs(center_csv, fermi_csv)
        assert len(rows) == 2
