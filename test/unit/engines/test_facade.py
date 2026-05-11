"""Smoke tests for the ``md_analysis.engines`` package-level facade.

These tests pin the public symbol surface so that future refactors do
not silently break ``from md_analysis.engines import ...`` consumers.
"""

from __future__ import annotations

import pytest


def test_public_symbols_importable_from_package_root() -> None:
    """All names listed in ``engines.__all__`` resolve at the package root."""
    from md_analysis.engines import (
        CP2KParser,
        ConstraintMDParser,
        ParserInferenceError,
        PotentialFrame,
        get_parser,
        infer_parser,
        register_parser,
        resolve_parser,
    )

    # Every symbol is non-None so e.g. typos in the __init__ re-export
    # would surface here even if the import line itself succeeded
    # (which it would for a stale alias).
    assert ConstraintMDParser is not None
    assert ParserInferenceError is not None
    assert CP2KParser is not None
    assert PotentialFrame is not None
    assert callable(register_parser)
    assert callable(get_parser)
    assert callable(infer_parser)
    assert callable(resolve_parser)


def test_default_registration_makes_cp2k_resolvable() -> None:
    """Importing the package must auto-register the CP2K parser."""
    from md_analysis.engines import CP2KParser, get_parser

    parser = get_parser("cp2k")
    assert isinstance(parser, CP2KParser)
    assert parser.name == "cp2k"


def test_unknown_parser_name_raises_parser_inference_error() -> None:
    """``get_parser`` on an unregistered name raises ``ParserInferenceError``."""
    from md_analysis.engines import ParserInferenceError, get_parser

    with pytest.raises(ParserInferenceError):
        get_parser("definitely-not-a-real-engine")


def test_vasp_placeholder_is_not_auto_registered() -> None:
    """VASP must NOT be in the default registry — calling ``infer_parser``
    on a directory the CP2K parser doesn't recognise should raise rather
    than silently return the VASP placeholder."""
    from md_analysis.engines import ParserInferenceError, get_parser

    with pytest.raises(ParserInferenceError):
        get_parser("vasp")


def test_explicit_vasp_construction_raises_not_implemented(tmp_path) -> None:
    """The VASP placeholder implements the Protocol shape but every
    method must raise ``NotImplementedError`` until Phase 8 fills it in."""
    from md_analysis.engines.vasp import VASPParser

    parser = VASPParser()
    assert parser.name == "vasp"

    with pytest.raises(NotImplementedError):
        parser.is_constraint_directory(tmp_path)

    with pytest.raises(NotImplementedError):
        parser.parse_metadata(tmp_path)

    with pytest.raises(NotImplementedError):
        parser.parse_lambda_series(tmp_path)
