"""Tests for ``CellSpec`` model + ``engines.cp2k.read_cell`` facade.

Phase 4 Commit 2 (D7 Layer 1).  Covers:

  - ``CellSpec.abc_ang`` derived formula (orthogonal + non-orthogonal)
  - ``CellSpec.is_orthorhombic`` boundary behaviour (tol = 1e-6 A)
  - ``read_cell`` suffix-based dispatch
    (``.restart`` / ``.restart.bak-N`` / ``md.inp``)
  - Numerical agreement with the existing
    ``parse_abc_from_restart`` test baseline.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.engines import CellSpec, read_cell

REPO_ROOT = Path(__file__).resolve().parents[3]
SG_DATA = REPO_ROOT / "data_example" / "sg"


# =========================================================================
# CellSpec.abc_ang
# =========================================================================


class TestAbcAng:
    """Derived ``abc_ang`` matches the norms of the three lattice rows."""

    def test_orthogonal_returns_diagonal(self) -> None:
        m = np.diag([10.0, 20.0, 30.0])
        spec = CellSpec(cell_matrix_ang=m)
        np.testing.assert_allclose(spec.abc_ang, (10.0, 20.0, 30.0), atol=1e-12)

    def test_non_orthogonal_returns_row_norms_not_diagonal(self) -> None:
        # Non-orthogonal cell: row 0 = (3, 4, 0) has norm 5, NOT 3.
        m = np.array(
            [
                [3.0, 4.0, 0.0],
                [0.0, 5.0, 0.0],
                [0.0, 0.0, 7.0],
            ]
        )
        spec = CellSpec(cell_matrix_ang=m)
        np.testing.assert_allclose(spec.abc_ang, (5.0, 5.0, 7.0), atol=1e-12)


# =========================================================================
# CellSpec.is_orthorhombic
# =========================================================================


class TestIsOrthorhombic:
    """Off-diagonal tolerance is 1e-6 A (matches parse_abc_from_restart)."""

    def test_diagonal_matrix(self) -> None:
        spec = CellSpec(cell_matrix_ang=np.diag([10.0, 10.0, 30.0]))
        assert spec.is_orthorhombic is True

    def test_off_diagonal_above_tolerance_is_false(self) -> None:
        m = np.diag([10.0, 10.0, 30.0]).astype(float)
        m[0, 1] = 1e-5  # 10x tolerance
        spec = CellSpec(cell_matrix_ang=m)
        assert spec.is_orthorhombic is False

    def test_off_diagonal_within_tolerance_is_true(self) -> None:
        m = np.diag([10.0, 10.0, 30.0]).astype(float)
        m[0, 1] = 1e-7  # well below tolerance
        spec = CellSpec(cell_matrix_ang=m)
        assert spec.is_orthorhombic is True


# =========================================================================
# read_cell — restart dispatch
# =========================================================================


class TestReadCellFromRestart:
    """Suffix dispatch for ``.restart`` files (with and without bak)."""

    def test_angle_restart(self) -> None:
        spec = read_cell(SG_DATA / "angle" / "slowgrowth-1.restart")
        # Baseline matches test_slowgrowth_parser.py:88 for the same file.
        assert spec.abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)
        assert spec.is_orthorhombic is True
        np.testing.assert_allclose(
            spec.cell_matrix_ang,
            np.diag([10.2239, 10.2239, 26.422]),
            atol=1e-4,
        )

    def test_distance_restart_bak(self) -> None:
        spec = read_cell(
            SG_DATA / "distance" / "slowgrowth-1.restart.bak-1"
        )
        assert spec.abc_ang == pytest.approx((10.2239, 10.2239, 26.422), rel=1e-4)
        assert spec.is_orthorhombic is True


# =========================================================================
# read_cell — md.inp dispatch
# =========================================================================


class TestReadCellFromMdInp:
    """Suffix dispatch for non-restart files (md.inp / *.inp)."""

    def test_md_inp_orthogonal(self, tmp_path: Path) -> None:
        md_inp = tmp_path / "md.inp"
        md_inp.write_text(
            "&CELL\n"
            "  ABC [angstrom] 12.5 12.5 28.0\n"
            "&END CELL\n",
            encoding="utf-8",
        )
        spec = read_cell(md_inp)
        assert spec.abc_ang == pytest.approx((12.5, 12.5, 28.0), abs=1e-12)
        assert spec.is_orthorhombic is True
        np.testing.assert_allclose(
            spec.cell_matrix_ang,
            np.diag([12.5, 12.5, 28.0]),
            atol=1e-12,
        )
