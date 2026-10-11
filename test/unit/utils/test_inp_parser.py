"""Tests for ``md_analysis.utils.formats.cp2k.inp`` (refactor steps 1-3).

Anchor groups follow the refactor plan §4 (step-1 group): preprocessing
workarounds (W3/W1/W2), lazy singletons, xmlspec pass-through, exception
wrapping, two-form section normalisation, W1 x line-continuation, W2 unit
anchors, W3 case-insensitivity, comment cleaning, stabilise/self-check,
alias normalisation, ``@`` directive gate, value types and the
``set_cell_abc`` vector-form no-op.

Step-2 group (see §4 / §3.2): CELL ``ABC`` unit three-form anchors,
``cell._abc_ang_from_cell_tree`` vector-form anchors, and the md_inp
fourth-message anchor (GLOBAL-only fixture, migrated verbatim from
``test_cell_parser.py``).

Step-3 group (see §4 / §3.2): the rewritten fixed-atoms W1 anchor
(formerly ``test_slowgrowth_parser.py::test_fixed_atoms_range_expansion``,
which directly tested the deleted regex helper), the disordered-LIST
sorted-tuple anchor through ``parse_colvar_restart``, and the
non-orthogonal dual-path anchor (``parse_abc_from_restart`` and
``parse_colvar_restart`` on the same fixture).

Real-data anchors use ``data_example/sg/distance_combinedCV`` (kept
untouched); all synthetic anchors are inline strings, so no synthetic data
files are needed under ``test/unit/data/cp2k_inp/``.
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

import pytest

from md_analysis.utils.formats.cp2k import inp as inp_module
from md_analysis.utils.formats.cp2k.cell import (
    CellParseError,
    _abc_ang_from_cell_tree,
    parse_abc_from_md_inp,
    parse_abc_from_restart,
)
from md_analysis.utils.formats.cp2k.colvar import parse_colvar_restart
from md_analysis.utils.formats.cp2k.inp import (
    InpGenerateError,
    InpParseError,
    generate_inp,
    get_keyword,
    get_section,
    iter_collectives,
    parse_inp,
    set_cell_abc,
    set_keyword,
)

DATA = Path(__file__).resolve().parents[3] / "data_example" / "sg"
SG_INP = DATA / "distance_combinedCV" / "sg.inp"
SG_RESTART = DATA / "distance_combinedCV" / "slowgrowth-1.restart"


def _constraint_text(*keyword_lines: str) -> str:
    """Minimal MOTION/CONSTRAINT/COLLECTIVE input carrying *keyword_lines*."""
    lines = ["&MOTION", "  &CONSTRAINT", "    &COLLECTIVE", "      COLVAR 1"]
    lines += [f"      {kw}" for kw in keyword_lines]
    lines += ["    &END COLLECTIVE", "  &END CONSTRAINT", "&END MOTION", ""]
    return "\n".join(lines)


def _fixed_atoms_text(list_line: str) -> str:
    return (
        "&MOTION\n"
        "  &CONSTRAINT\n"
        "    &FIXED_ATOMS\n"
        f"      {list_line}\n"
        "    &END FIXED_ATOMS\n"
        "  &END CONSTRAINT\n"
        "&END MOTION\n"
    )


# ---------------------------------------------------------------------------
# Lazy singletons (cit + pint must stay off the import path)
# ---------------------------------------------------------------------------


class TestLazyImports:
    """cit and pint are imported on first parse only, never at module import."""

    def test_module_import_does_not_load_cit_or_pint(self):
        """Hermetic subprocess check: neither cit nor pint in sys.modules."""
        code = (
            "import sys\n"
            "import md_analysis.utils.formats.cp2k.inp\n"
            "assert 'cp2k_input_tools' not in sys.modules, 'cit imported eagerly'\n"
            "assert 'pint' not in sys.modules, 'pint imported eagerly'\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr

    def test_parser_and_ureg_singletons(self):
        assert inp_module._get_parser(None) is inp_module._get_parser(None)
        assert inp_module._get_generator(None) is inp_module._get_generator(None)
        assert inp_module._get_ureg() is inp_module._get_ureg()


class TestXmlspecPassThrough:
    """xmlspec parameterisation of parser and generator."""

    def test_explicit_xmlspec_parses_identically(self):
        from cp2k_input_tools import DEFAULT_CP2K_INPUT_XML

        text = "&GLOBAL\n  PROJECT_NAME x\n&END GLOBAL\n"
        default = parse_inp(text)
        explicit = parse_inp(text, xmlspec=DEFAULT_CP2K_INPUT_XML)
        assert explicit == default == {"global": {"project_name": "x"}}

    def test_xmlspec_keys_cached_separately(self):
        from cp2k_input_tools import DEFAULT_CP2K_INPUT_XML

        assert inp_module._get_parser(DEFAULT_CP2K_INPUT_XML) is not inp_module._get_parser(None)
        assert inp_module._get_generator(DEFAULT_CP2K_INPUT_XML) is not inp_module._get_generator(None)

    def test_generate_with_xmlspec_round_trips(self):
        from cp2k_input_tools import DEFAULT_CP2K_INPUT_XML

        tree = parse_inp(_constraint_text("TARGET 1.5"), xmlspec=DEFAULT_CP2K_INPUT_XML)
        text = generate_inp(tree, xmlspec=DEFAULT_CP2K_INPUT_XML)
        assert parse_inp(text, xmlspec=DEFAULT_CP2K_INPUT_XML) == tree


# ---------------------------------------------------------------------------
# Input dispatch and exception wrapping
# ---------------------------------------------------------------------------


class TestInputDispatchAndErrors:
    def test_pathlike_and_plain_str_path(self, tmp_path):
        f = tmp_path / "job.inp"
        f.write_text("&GLOBAL\n  PROJECT_NAME viapath\n&END GLOBAL\n", encoding="utf-8")
        assert parse_inp(f)["global"]["project_name"] == "viapath"
        assert parse_inp(str(f))["global"]["project_name"] == "viapath"

    def test_text_with_newline_parsed_inline(self, tmp_path):
        # a str containing a newline is text, never touched as a path
        text = "&GLOBAL\n  PROJECT_NAME inline\n&END GLOBAL\n"
        assert parse_inp(text)["global"]["project_name"] == "inline"

    def test_str_without_newline_missing_file_raises_with_hint(self):
        with pytest.raises(FileNotFoundError, match="no such file"):
            parse_inp("definitely/not/a/file.inp")

    def test_unsupported_input_type_raises_typeerror(self):
        with pytest.raises(TypeError, match="path"):
            parse_inp(42)

    def test_cit_error_wrapped_with_line_info(self):
        with pytest.raises(InpParseError, match="line 2") as excinfo:
            parse_inp("&GLOBAL\n  BOGUS_KEYWORD 1\n&END GLOBAL\n")
        assert type(excinfo.value.__cause__).__name__ == "InvalidNameError"
        assert "cannot parse CP2K input" in str(excinfo.value)

    def test_inpparse_error_is_mdanalysis_error(self):
        assert issubclass(InpParseError, Exception)
        from md_analysis.exceptions import MDAnalysisError

        assert issubclass(InpParseError, MDAnalysisError)
        assert issubclass(InpGenerateError, MDAnalysisError)


# ---------------------------------------------------------------------------
# W3: comma normalisation on LIST lines
# ---------------------------------------------------------------------------


class TestListCommaNormalisation:
    def test_list_commas_case_insensitive_with_range(self):
        """`List 1,2,7..10` -> integer tuple (W3 + W1 chain)."""
        tree = parse_inp(_fixed_atoms_text("List 1,2,7..10"))
        assert tree["motion"]["constraint"]["fixed_atoms"]["list"] == (1, 2, 7, 8, 9, 10)

    def test_only_list_lines_touched(self):
        raw = "  OTHER 1,2\n  LIST 1,2\n  listx 1,2\n"
        assert inp_module._normalize_list_commas(raw) == (
            "  OTHER 1,2\n  LIST 1 2\n  listx 1,2\n"
        )

    def test_crlf_normalised_first(self):
        text = "&GLOBAL\r\n  PROJECT_NAME crlf\r\n&END GLOBAL\r\n"
        assert parse_inp(text)["global"]["project_name"] == "crlf"


# ---------------------------------------------------------------------------
# W1: integer range expansion
# ---------------------------------------------------------------------------


class TestRangeExpansion:
    def test_range_expanded(self):
        tree = parse_inp(_fixed_atoms_text("LIST 2..4"))
        assert tree["motion"]["constraint"]["fixed_atoms"]["list"] == (2, 3, 4)

    def test_range_with_line_continuation(self):
        """W1 expansion is inline, so `\\` continuations keep working."""
        tree = parse_inp(_fixed_atoms_text("LIST 2..4 \\\n           5"))
        assert tree["motion"]["constraint"]["fixed_atoms"]["list"] == (2, 3, 4, 5)

    def test_no_range_inside_numbers(self):
        assert inp_module._expand_ranges("EPS 1.2..3.4 X2..4") == "EPS 1.2..3.4 X2..4"

    def test_decreasing_range_left_alone(self):
        assert inp_module._expand_ranges("5..3") == "5..3"


# ---------------------------------------------------------------------------
# W2: internal-unit conversion on TARGET / TARGET_GROWTH
# ---------------------------------------------------------------------------

_W2_LITERALS = [
    # (keyword line, tree key, expected value) — plan §4 synthetic anchors,
    # fs/ps map to atomic_unit_of_time (1 fs = 41.34137333518245 aut)
    ("TARGET_GROWTH [fs^-1*angstrom] 1E-3", "target_growth", 4.5710289043883664e-05),
    ("TARGET_GROWTH [fs^-1*angstrom] -0.5E-3", "target_growth", -2.2855144521941832e-05),
    ("TARGET_GROWTH [fs^-1*deg] 1E-2", "target_growth", 4.2217495723806265e-06),
    ("TARGET [angstrom] 1.0", "target", 1.8897261246222794),
    ("TARGET [deg] 1.0", "target", 0.017453292519943295),
]


class TestInternalUnitConversion:
    @pytest.mark.parametrize(("line", "key", "expected"), _W2_LITERALS)
    def test_synthetic_literals(self, line, key, expected):
        tree = parse_inp(_constraint_text(line))
        assert tree["motion"]["constraint"]["collective"][key] == pytest.approx(
            expected, rel=1e-12
        )

    def test_real_data_cross_check_against_restart(self):
        """sg.inp tree TARGET_GROWTH (converted) matches the .restart value
        (internal units, bare) to rel=1e-6 — direction/scale of the
        fs -> aut conversion is anchored by real data."""
        inp_tree = parse_inp(SG_INP)
        restart_tree = parse_inp(SG_RESTART)
        converted = inp_tree["motion"]["constraint"]["collective"]["target_growth"]
        reference = restart_tree["motion"]["constraint"]["collective"]["target_growth"]
        assert converted == pytest.approx(reference, rel=1e-6)

    def test_timestep_fs_lands_on_fs(self):
        """Reverse direction anchor: TIMESTEP (MD block, default unit fs)
        keeps [fs] values at ~1.0 — a swapped fs/aut mapping would give 41.3."""
        text = "&MOTION\n  &MD\n    TIMESTEP 1.0 [fs]\n  &END MD\n&END MOTION\n"
        assert parse_inp(text)["motion"]["md"]["timestep"] == pytest.approx(1.0)

    def test_bare_internal_value_untouched(self):
        tree = parse_inp(_constraint_text("TARGET 1.8897261246222794"))
        assert tree["motion"]["constraint"]["collective"]["target"] == 1.8897261246222794

    def test_unknown_unit_fails_fast(self):
        """Units W2 cannot translate are left in place, so cit rejects the
        unit token with its own line-numbered error (wrapped)."""
        with pytest.raises(InpParseError, match="line 5"):
            parse_inp(_constraint_text("TARGET [furlong] 1.0"))


# ---------------------------------------------------------------------------
# Comment cleaning (cit issue #111 guard)
# ---------------------------------------------------------------------------


class TestCommentCleaning:
    def test_inline_comment_with_quote_truncated(self):
        clean = inp_module._clean_comments
        assert clean("KEY val # don't") == "KEY val"
        assert clean("KEY val ! l'enfant") == "KEY val"
        assert clean('KEY val # say "hi"') == "KEY val"

    def test_inline_comment_without_quote_byte_untouched(self):
        clean = inp_module._clean_comments
        assert clean("KEY val # plain, comment") == "KEY val # plain, comment"
        assert clean("KEY val # trailing \\") == "KEY val # trailing \\"

    def test_whole_line_comment_byte_untouched(self):
        clean = inp_module._clean_comments
        assert clean("# whole line, untouched") == "# whole line, untouched"
        assert clean("  ! bang line it's") == "  ! bang line it's"

    def test_quote_inside_value_preserved(self):
        clean = inp_module._clean_comments
        assert clean('KEY "a#b" # note') == 'KEY "a#b" # note'

    def test_continuation_joined_comment_cleaned(self):
        """A comment line continuing a `\\` line is inline after cit's join,
        so it must be truncated (and the continuation marker kept)."""
        clean = inp_module._clean_comments
        assert clean("KEY \\\n# joined comment isn't lost") == "KEY \\\n"

    def test_quote_comments_parse_end_to_end(self):
        text = (
            "&GLOBAL\n"
            "  PROJECT_NAME H2O # don't touch this\n"
            "  RUN_TYPE MD\n"
            "  # whole line comment with quote: it's fine\n"
            "  ! bang comment with \"quote\" fine too\n"
            "&END GLOBAL\n"
        )
        tree = parse_inp(text)
        assert tree["global"]["project_name"] == "H2O"

    def test_quoted_hash_in_value_parses(self):
        text = '&GLOBAL\n  PROJECT_NAME "ab#cd" ! trailing\n&END GLOBAL\n'
        assert parse_inp(text)["global"]["project_name"] == "ab#cd"


# ---------------------------------------------------------------------------
# G1: @ preprocessor directive gate
# ---------------------------------------------------------------------------


class TestAtDirectiveGate:
    def test_at_if_rejected(self):
        with pytest.raises(InpParseError, match="preprocessor"):
            parse_inp("&GLOBAL\n@IF (${GO} == 2)\n&END GLOBAL\n")

    def test_at_set_rejected_with_line_number(self):
        with pytest.raises(InpParseError, match="line 3.*preprocessor"):
            parse_inp("&GLOBAL\n  RUN_TYPE MD\n  @SET GO 2\n&END GLOBAL\n")

    def test_commented_directives_not_gated(self):
        text = (
            "&GLOBAL\n"
            "  RUN_TYPE MD\n"
            "  # @SET GO 2\n"
            "  ! @IF (${X} == 2)\n"
            "&END GLOBAL\n"
        )
        assert parse_inp(text)["global"]["run_type"] == "md"

    def test_at_inside_value_not_gated(self):
        text = "&GLOBAL\n  PROJECT a@b\n&END GLOBAL\n"
        assert parse_inp(text)["global"]["project_name"] == "a@b"

    def test_continuation_line_not_gated(self):
        """A `\\`-continued line starting with `@` is a value continuation,
        not a directive (cit joins it mid-line); the gate must not fire."""
        inp_module._gate_at_directives("KEY value \\\n@SET x 1")  # no raise
        with pytest.raises(InpParseError, match="preprocessor"):
            inp_module._gate_at_directives("KEY value\n  @SET x 1")


# ---------------------------------------------------------------------------
# Stabilise + post-generation self-check (generator LIFO compensation)
# ---------------------------------------------------------------------------


class TestStabilizeAndSelfCheck:
    def test_real_fixture_nested_colvar_order_preserved(self):
        """distance_combinedCV/sg.inp: generate -> re-parse keeps the nested
        COMBINE_COLVAR atom pairs [(270,229),(228,229)] (semantic order)."""
        tree = parse_inp(SG_INP)
        text = generate_inp(tree)  # self-check must pass silently
        reparsed = parse_inp(text)
        nested = reparsed["force_eval"]["subsys"]["colvar"]["combine_colvar"]["colvar"]
        atoms = [block["distance"]["atoms"] for block in nested]
        assert atoms == [(270, 229), (228, 229)]

    def test_skipping_stabilize_raises_generate_error(self, monkeypatch):
        monkeypatch.setattr(inp_module, "_stabilize_tree", lambda node: node)
        tree = parse_inp(SG_INP)
        with pytest.raises(InpGenerateError, match="round-trip.*colvar"):
            generate_inp(tree)

    def test_repeated_keyword_order_preserved(self):
        text = (
            "&FORCE_EVAL\n"
            "  &SUBSYS\n"
            "    &KIND Pt\n"
            "      BASIS_SET ORB first-orb\n"
            "      BASIS_SET AUX1 second-aux\n"
            "      POTENTIAL GTH-PBE-q10\n"
            "    &END KIND\n"
            "  &END SUBSYS\n"
            "&END FORCE_EVAL\n"
        )
        tree = parse_inp(text)
        basis = tree["force_eval"]["subsys"]["kind"]["Pt"]["basis_set"]
        assert basis == [("ORB", "first-orb"), ("AUX1", "second-aux")]
        reparsed = parse_inp(generate_inp(tree))
        assert reparsed["force_eval"]["subsys"]["kind"]["Pt"]["basis_set"] == basis

    def test_repeated_sections_must_use_list_containers(self):
        dual = (
            "&MOTION\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n      COLVAR 1\n    &END COLLECTIVE\n"
            "    &COLLECTIVE\n      COLVAR 2\n    &END COLLECTIVE\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        tree = parse_inp(dual)
        blocks = tree["motion"]["constraint"]["collective"]
        assert isinstance(blocks, list) and len(blocks) == 2
        generate_inp(tree)  # list container OK

        tree["motion"]["constraint"]["collective"] = tuple(blocks)
        with pytest.raises(InpGenerateError):
            generate_inp(tree)


# ---------------------------------------------------------------------------
# Alias normalisation (canonical key names)
# ---------------------------------------------------------------------------


class TestAliasNormalisation:
    def test_project_alias_normalised_to_canonical_key(self):
        tree = parse_inp("&GLOBAL\n  PROJECT via-alias\n&END GLOBAL\n")
        assert tree["global"] == {"project_name": "via-alias"}
        assert "project" not in tree["global"]

    def test_set_canonical_key_regenerates_without_repetition(self):
        tree = parse_inp("&GLOBAL\n  PROJECT old\n&END GLOBAL\n")
        set_keyword(tree, "global.project_name", "new")
        reparsed = parse_inp(generate_inp(tree))
        # a single canonical key; alias key never re-appears (which would
        # make a re-parse raise NameRepetitionError later on)
        assert reparsed["global"] == {"project_name": "new"}


# ---------------------------------------------------------------------------
# Tree path primitives (dict/list two-form normalisation)
# ---------------------------------------------------------------------------


class TestTreePrimitives:
    def test_get_keyword_dict_form(self):
        tree = parse_inp(SG_INP)
        assert get_keyword(tree, "motion.md.timestep") == 1.0
        assert get_keyword(tree, "force_eval.subsys.kind.O.basis_set") == "DZVP-MOLOPT-SR-GTH"

    def test_get_keyword_list_form_descends_by_param(self):
        tree = {
            "force_eval": {
                "subsys": {
                    "kind": [
                        {"_": "O", "basis_set": "X"},
                        {"_": "Pt", "basis_set": "Y"},
                    ]
                }
            }
        }
        assert get_keyword(tree, "force_eval.subsys.kind.Pt.basis_set") == "Y"

    def test_get_section_index_path(self):
        tree = {"a": {"b": [{"_": "x"}, {"_": "y"}]}}
        assert get_section(tree, "a.b[1]._") == "y"

    def test_get_section_missing_raises_keyerror(self):
        with pytest.raises(KeyError, match="cell"):
            get_section({"global": {}}, "force_eval.subsys.cell")

    def test_set_keyword_dict_and_list_form(self):
        tree = {"force_eval": {"subsys": {"kind": [{"_": "Pt"}]}}}
        set_keyword(tree, "force_eval.subsys.kind.Pt.potential", "GTH")
        assert tree["force_eval"]["subsys"]["kind"][0]["potential"] == "GTH"
        set_keyword(tree, "force_eval.subsys.kind.Pt.potential", "GTH-PBE")
        assert tree["force_eval"]["subsys"]["kind"][0]["potential"] == "GTH-PBE"

    def test_set_keyword_requires_parent_section(self):
        with pytest.raises(KeyError):
            set_keyword({"global": {}}, "motion.md.timestep", 1.0)

    def test_iter_collectives_list_and_dict_forms(self):
        dual = (
            "&MOTION\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n      COLVAR 1\n    &END COLLECTIVE\n"
            "    &COLLECTIVE\n      COLVAR 2\n    &END COLLECTIVE\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        assert [cid for cid, _ in iter_collectives(parse_inp(dual))] == [1, 2]
        single = parse_inp(_constraint_text("TARGET 1.5"))
        assert [cid for cid, _ in iter_collectives(single)] == [1]


# ---------------------------------------------------------------------------
# Value types on the simplified tree
# ---------------------------------------------------------------------------


class TestValueTypes:
    def test_native_types(self):
        tree = parse_inp(
            "&MOTION\n"
            "  &MD\n    STEPS 100\n    TIMESTEP 0.5\n  &END MD\n"
            "  &CONSTRAINT\n"
            "    &COLLECTIVE\n"
            "      COLVAR 1\n      INTERMOLECULAR .TRUE.\n    &END COLLECTIVE\n"
            "  &END CONSTRAINT\n"
            "&END MOTION\n"
        )
        md = tree["motion"]["md"]
        assert type(md["steps"]) is int
        assert type(md["timestep"]) is float
        collective = tree["motion"]["constraint"]["collective"]
        assert collective["intermolecular"] is True

    def test_multi_value_keyword_is_tuple(self):
        tree = parse_inp(
            "&FORCE_EVAL\n  &SUBSYS\n    &CELL\n"
            "      ABC [angstrom] 10.0 11.0 12.0\n"
            "    &END CELL\n  &END SUBSYS\n&END FORCE_EVAL\n"
        )
        assert tree["force_eval"]["subsys"]["cell"]["abc"] == (10.0, 11.0, 12.0)


# ---------------------------------------------------------------------------
# set_cell_abc vector-form no-op
# ---------------------------------------------------------------------------


class TestSetCellAbc:
    def test_abc_written_on_abc_form_template(self):
        tree = parse_inp(SG_INP)
        set_cell_abc(tree, 11.0, 12.0, 13.0)
        assert tree["force_eval"]["subsys"]["cell"]["abc"] == (11.0, 12.0, 13.0)

    def test_vector_form_template_warns_and_noop(self, caplog):
        vector_tree = {
            "force_eval": {"subsys": {"cell": {"a": (1, 0, 0), "b": (0, 1, 0), "c": (0, 0, 1)}}}
        }
        with caplog.at_level(logging.WARNING, logger="md_analysis.utils.formats.cp2k.inp"):
            set_cell_abc(vector_tree, 11.0, 12.0, 13.0)
        assert any("ABC" in record.message for record in caplog.records)
        assert "abc" not in vector_tree["force_eval"]["subsys"]["cell"]

    def test_missing_cell_section_warns_and_noop(self, caplog):
        tree = {"global": {"project_name": "x"}}
        with caplog.at_level(logging.WARNING, logger="md_analysis.utils.formats.cp2k.inp"):
            set_cell_abc(tree, 11.0, 12.0, 13.0)
        assert any("ABC" in record.message for record in caplog.records)
        assert tree == {"global": {"project_name": "x"}}


# ---------------------------------------------------------------------------
# Step-2 anchors: cell.py on the cit tree path (plan §4 / §3.2)
# ---------------------------------------------------------------------------


def _cell_inp_text(cell_line: str) -> str:
    """Minimal FORCE_EVAL/SUBSYS/CELL input carrying one *cell_line*."""
    return (
        "&FORCE_EVAL\n"
        "  &SUBSYS\n"
        "    &CELL\n"
        f"      {cell_line}\n"
        "    &END CELL\n"
        "  &END SUBSYS\n"
        "&END FORCE_EVAL\n"
    )


def _abc_of(cell_line: str) -> tuple:
    tree = parse_inp(_cell_inp_text(cell_line))
    return get_keyword(tree, "force_eval.subsys.cell.abc")


class TestCellAbcUnitForms:
    """Step 2: CELL ``ABC`` lands on the same angstrom triple for the
    ``[angstrom]`` / ``[bohr]`` / bare-token forms (the XML default unit of
    ``ABC`` is angstrom; cit converts any explicit unit token)."""

    EXPECTED = (10.2239, 10.2239, 26.422)

    def test_angstrom_token(self):
        abc = _abc_of("ABC [angstrom] 10.2239 10.2239 26.422")
        assert abc == pytest.approx(self.EXPECTED, rel=1e-12)

    def test_bare_token_defaults_to_angstrom(self):
        abc = _abc_of("ABC 10.2239 10.2239 26.422")
        assert abc == pytest.approx(self.EXPECTED, rel=1e-12)

    def test_bohr_token_converted_to_same_values(self):
        # bohr values = EXPECTED / 0.5291772109039775 (cit bohr->angstrom),
        # full float precision so the round-trip is exact to rel=1e-12.
        abc = _abc_of("ABC [bohr] 19.320370925525722 19.320370925525722 49.930343664769865")
        assert abc == pytest.approx(self.EXPECTED, rel=1e-12)


class TestAbcAngFromCellTree:
    """Step 2: ``_abc_ang_from_cell_tree`` takes ``cell.{a,b,c}`` vector
    diagonals (restart-style tree, no ``abc`` keyword), keeps the
    orthogonality policy (off-diagonal tolerance 1e-6) and the legacy
    error messages (end-to-end coverage lives in test_cell_parser.py)."""

    @staticmethod
    def _tree(**overrides) -> dict:
        cell = {
            "a": (10.2239, 0.0, 0.0),
            "b": (0.0, 10.2239, 0.0),
            "c": (0.0, 0.0, 26.422),
        }
        cell.update(overrides)
        return {"force_eval": {"subsys": {"cell": cell}}}

    def test_parsed_vector_form_tree_diagonals(self):
        text = (
            "&FORCE_EVAL\n  &SUBSYS\n    &CELL\n"
            "      A  10.2239   0.0    0.0\n"
            "      B   0.0   10.2239   0.0\n"
            "      C   0.0     0.0   26.422\n"
            "    &END CELL\n  &END SUBSYS\n&END FORCE_EVAL\n"
        )
        tree = parse_inp(text)
        assert "abc" not in tree["force_eval"]["subsys"]["cell"]  # vector form
        assert _abc_ang_from_cell_tree(tree, "x.restart") == (
            10.2239, 10.2239, 26.422,
        )

    def test_off_diagonal_below_tolerance_passes(self):
        tree = self._tree(b=(1e-9, 10.2239, 0.0))
        assert _abc_ang_from_cell_tree(tree, "x.restart") == (
            10.2239, 10.2239, 26.422,
        )

    def test_missing_cell_section_raises(self):
        with pytest.raises(CellParseError, match="No &CELL"):
            _abc_ang_from_cell_tree({"global": {}}, "x.restart")

    def test_missing_vector_raises(self):
        tree = {"force_eval": {"subsys": {"cell": {"a": (1, 0, 0), "b": (0, 1, 0)}}}}
        with pytest.raises(CellParseError, match="Vector C not found"):
            _abc_ang_from_cell_tree(tree, "x.restart")

    def test_non_orthogonal_raises(self):
        tree = self._tree(a=(10.0, 5.0, 0.0))
        with pytest.raises(CellParseError, match="Non-orthogonal"):
            _abc_ang_from_cell_tree(tree, "x.restart")


class TestMdInpMissingAbc:
    """Step 2: md_inp fourth message — parse succeeds but the tree lacks
    ``force_eval.subsys.cell.abc`` (GLOBAL-only fixture needs no
    legalisation; migrated verbatim from test_cell_parser.py)."""

    def test_global_only_inp_raises_cannot_find(self, tmp_path):
        f = tmp_path / "bad.inp"
        f.write_text("&GLOBAL\n  PROJECT test\n&END GLOBAL\n", encoding="utf-8")
        with pytest.raises(CellParseError, match="Cannot find"):
            parse_abc_from_md_inp(f)


# ---------------------------------------------------------------------------
# Step-3 anchors: colvar.py on the cit tree path (plan §4 / §3.2)
# ---------------------------------------------------------------------------

_ORTHO_CELL = (
    "      A  10.0  0.0  0.0\n"
    "      B   0.0 10.0  0.0\n"
    "      C   0.0  0.0 20.0"
)

_NON_ORTHO_CELL = (
    "      A  10.0  5.0  0.0\n"
    "      B   0.0 10.0  0.0\n"
    "      C   0.0  0.0 20.0"
)


def _colvar_restart_text(cell_lines: str, constraint_extra: str = "") -> str:
    """Complete synthetic .restart-shaped input for the colvar read path."""
    return (
        "&GLOBAL\n  PROJECT_NAME colvar\n&END GLOBAL\n"
        "&FORCE_EVAL\n  &SUBSYS\n    &CELL\n"
        f"{cell_lines}\n"
        "    &END CELL\n  &END SUBSYS\n&END FORCE_EVAL\n"
        "&MOTION\n  &MD\n"
        "    STEPS 100\n    TIMESTEP 1.0\n"
        "    STEP_START_VAL 0\n    TIME_START_VAL 0.0\n"
        "  &END MD\n  &CONSTRAINT\n"
        "    &COLLECTIVE\n      COLVAR 1\n      TARGET 1.5\n"
        "      TARGET_GROWTH 0.001\n    &END COLLECTIVE\n"
        f"{constraint_extra}"
        "  &END CONSTRAINT\n&END MOTION\n"
    )


class TestFixedAtomsW1Expansion:
    """Step 3: rewrite of ``test_slowgrowth_parser.py::
    test_fixed_atoms_range_expansion``, which directly tested the deleted
    regex helper ``_parse_fixed_atoms_list``.  Range expansion now lives
    in the inp.py pipeline (W1), so the anchor parses a complete legal
    input and checks the tree tuple (the legacy fixture, a bare
    ``&CONSTRAINT`` fragment, was legalised by wrapping it in ``&MOTION``).
    """

    def test_range_expanded_in_tree(self):
        tree = parse_inp(_fixed_atoms_text("LIST 1 3..5 10"))
        assert tree["motion"]["constraint"]["fixed_atoms"]["list"] == (1, 3, 4, 5, 10)


class TestFixedAtomsSortedThroughColvar:
    """Step 3: ``parse_colvar_restart`` sorts the W1-expanded LIST into an
    ascending tuple (explicit sorted -> tuple contract, plan §3.2)."""

    def test_disordered_list_sorted(self, tmp_path):
        f = tmp_path / "disordered.restart"
        f.write_text(
            _colvar_restart_text(
                _ORTHO_CELL,
                constraint_extra=(
                    "    &FIXED_ATOMS\n      LIST 7 2 5..6\n    &END FIXED_ATOMS\n"
                ),
            ),
            encoding="utf-8",
        )
        raw = parse_colvar_restart(f)
        assert raw.fixed_atom_indices == (2, 5, 6, 7)


class TestColvarNonOrthogonalDualAnchor:
    """Step 3: the same non-orthogonal fixture fails on both cell read
    paths with ``CellParseError`` "Non-orthogonal" — the shared
    orthogonality policy of ``_abc_ang_from_cell_tree`` is not lost on
    the colvar path (plan §3.2, Q6)."""

    def test_both_paths_raise_non_orthogonal(self, tmp_path):
        f = tmp_path / "non_ortho.restart"
        f.write_text(_colvar_restart_text(_NON_ORTHO_CELL), encoding="utf-8")
        with pytest.raises(CellParseError, match="Non-orthogonal"):
            parse_abc_from_restart(f)
        with pytest.raises(CellParseError, match="Non-orthogonal"):
            parse_colvar_restart(f)
