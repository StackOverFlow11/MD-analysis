"""Shared pytest fixtures and helpers for MD Analysis tests."""

from __future__ import annotations

import pytest

from md_analysis.utils.formats.cp2k.cell import parse_abc_from_md_inp
from md_analysis.utils.formats.cp2k.inp import get_keyword, parse_inp

__all__ = ["parse_abc_from_md_inp"]


@pytest.fixture
def assert_inp_tree():
    """Tree-assertion helper for CP2K inp text (refactor plan §4).

    Usage::

        assert_inp_tree(text, {"force_eval.subsys.cell.abc": (5.0, 5.0, 10.0)})

    The text is parsed once through the cit pipeline
    (``md_analysis.utils.formats.cp2k.inp.parse_inp``) and every dotted
    path (``.`` / ``[n]`` syntax of ``inp.get_keyword``) is compared
    against the expected value.  Assertions are made on the tree, not
    the generated text, so keyword order / casing / comment loss in the
    generator output never makes a test brittle.
    """

    def _assert(text_or_path, expectations: dict) -> None:
        tree = parse_inp(text_or_path)
        for path, expected in expectations.items():
            actual = get_keyword(tree, path)
            assert actual == expected, (
                f"tree path {path!r}: expected {expected!r}, "
                f"got {actual!r}"
            )

    return _assert
