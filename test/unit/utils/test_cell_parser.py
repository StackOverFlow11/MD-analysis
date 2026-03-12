"""Tests for md_analysis.utils.CellParser."""

from __future__ import annotations

import pytest

from md_analysis.utils.RestartParser.CellParser import CellParseError, parse_abc_from_md_inp, parse_abc_from_restart


_VALID_RESTART = """\
&FORCE_EVAL
  &SUBSYS
    &CELL
      A     1.0223900000000002E+01    0.0000000000000000E+00    0.0000000000000000E+00
      B     0.0000000000000000E+00    1.0223900000000002E+01    0.0000000000000000E+00
      C     0.0000000000000000E+00    0.0000000000000000E+00    2.6422000000000001E+01
      PERIODIC  XYZ
      MULTIPLE_UNIT_CELL  1 1 1
    &END CELL
    &COORD
      Cu  0.0  0.0  0.0
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""


class TestParseAbcFromRestart:
    """Tests for parse_abc_from_restart."""

    def test_parse_valid_restart(self, tmp_path):
        """Correct diagonal elements are extracted."""
        f = tmp_path / "md.restart"
        f.write_text(_VALID_RESTART)
        a, b, c = parse_abc_from_restart(f)
        assert a == pytest.approx(10.2239, rel=1e-5)
        assert b == pytest.approx(10.2239, rel=1e-5)
        assert c == pytest.approx(26.422, rel=1e-4)

    def test_parse_missing_cell_raises(self, tmp_path):
        """CellParseError if no &CELL block exists."""
        f = tmp_path / "bad.restart"
        f.write_text("&MOTION\n  STEPS 1000\n&END MOTION\n")
        with pytest.raises(CellParseError, match="No &CELL"):
            parse_abc_from_restart(f)

    def test_parse_missing_vector_raises(self, tmp_path):
        """CellParseError if A/B/C vector line is missing."""
        content = (
            "&CELL\n"
            "  A  10.0  0.0  0.0\n"
            "  B   0.0  10.0  0.0\n"
            "  PERIODIC XYZ\n"
            "&END CELL\n"
        )
        f = tmp_path / "missing_c.restart"
        f.write_text(content)
        with pytest.raises(CellParseError, match="Vector C not found"):
            parse_abc_from_restart(f)

    def test_parse_non_orthogonal_raises(self, tmp_path):
        """CellParseError if off-diagonal elements are nonzero."""
        content = (
            "&CELL\n"
            "  A  10.0  5.0  0.0\n"
            "  B   0.0  10.0  0.0\n"
            "  C   0.0   0.0  20.0\n"
            "&END CELL\n"
        )
        f = tmp_path / "non_ortho.restart"
        f.write_text(content)
        with pytest.raises(CellParseError, match="Non-orthogonal"):
            parse_abc_from_restart(f)

    def test_parse_with_surrounding_sections(self, tmp_path):
        """Only &CELL block is extracted; other sections are ignored."""
        content = (
            "&MOTION\n  STEPS 1000\n&END MOTION\n"
            "&FORCE_EVAL\n  &SUBSYS\n"
            "    &CELL\n"
            "      A  5.0  0.0  0.0\n"
            "      B  0.0  8.0  0.0\n"
            "      C  0.0  0.0  12.5\n"
            "    &END CELL\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n"
        )
        f = tmp_path / "complex.restart"
        f.write_text(content)
        a, b, c = parse_abc_from_restart(f)
        assert a == pytest.approx(5.0)
        assert b == pytest.approx(8.0)
        assert c == pytest.approx(12.5)


class TestParseAbcFromMdInp:
    """Tests for parse_abc_from_md_inp."""

    def test_parse_valid_md_inp(self, tmp_path):
        """Correct ABC values are extracted."""
        f = tmp_path / "md.inp"
        f.write_text(
            "&SUBSYS\n"
            "  &CELL\n"
            "    ABC [angstrom] 10.2239 10.2239 26.422\n"
            "  &END CELL\n"
            "&END SUBSYS\n"
        )
        a, b, c = parse_abc_from_md_inp(f)
        assert a == pytest.approx(10.2239)
        assert b == pytest.approx(10.2239)
        assert c == pytest.approx(26.422)

    def test_parse_missing_abc_raises(self, tmp_path):
        """CellParseError if no ABC line exists."""
        f = tmp_path / "bad.inp"
        f.write_text("&GLOBAL\n  PROJECT test\n&END GLOBAL\n")
        with pytest.raises(CellParseError, match="Cannot find"):
            parse_abc_from_md_inp(f)

    def test_parse_scientific_notation(self, tmp_path):
        """Fortran-style scientific notation is parsed correctly."""
        f = tmp_path / "md.inp"
        f.write_text("  ABC [angstrom] 1.0E+01 2.0E+01 3.0E+01\n")
        a, b, c = parse_abc_from_md_inp(f)
        assert a == pytest.approx(10.0)
        assert b == pytest.approx(20.0)
        assert c == pytest.approx(30.0)
