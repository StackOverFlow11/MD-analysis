"""Tests for the engine-agnostic parser layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.engines import (
    CP2KParser,
    ConstraintMDParser,
    ParserInferenceError,
    get_parser,
    infer_parser,
    register_parser,
    resolve_parser,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
TI_POINT = (
    REPO_ROOT / "data_example" / "ti" / "double_cv" / "1k" / "ti_target_0.031369"
)
SG_DIR = REPO_ROOT / "data_example" / "sg" / "angle"


# ---------------------------------------------------------------------------
# Protocol conformance
# ---------------------------------------------------------------------------


def test_cp2k_parser_satisfies_protocol():
    parser = CP2KParser()
    assert isinstance(parser, ConstraintMDParser)
    assert parser.name == "cp2k"


# ---------------------------------------------------------------------------
# CP2KParser on real fixtures
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_cp2k_parser_recognises_ti_point():
    parser = CP2KParser()
    assert parser.is_constraint_directory(TI_POINT)


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_cp2k_parse_metadata_extracts_target():
    parser = CP2KParser()
    meta = parser.parse_metadata(TI_POINT)
    # Directory name encodes target = 0.031369; restart should agree
    assert meta.colvars.primary.target_au == pytest.approx(0.031369, abs=1e-5)
    assert meta.timestep_fs > 0


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_cp2k_parse_lambda_series_returns_array():
    parser = CP2KParser()
    log = parser.parse_lambda_series(TI_POINT)
    assert log.n_steps > 0
    assert isinstance(log.collective_shake, np.ndarray)
    assert log.collective_shake.shape == (log.n_steps,)


@pytest.mark.skipif(
    not SG_DIR.exists(), reason="SG fixture missing"
)
def test_cp2k_parser_recognises_sg_directory():
    parser = CP2KParser()
    assert parser.is_constraint_directory(SG_DIR)


# ---------------------------------------------------------------------------
# Negative cases
# ---------------------------------------------------------------------------


def test_cp2k_parser_rejects_empty_dir(tmp_path):
    parser = CP2KParser()
    assert parser.is_constraint_directory(tmp_path) is False


def test_cp2k_parser_rejects_nonexistent_dir(tmp_path):
    parser = CP2KParser()
    assert parser.is_constraint_directory(tmp_path / "missing") is False


def test_cp2k_parser_skips_bak_and_wfn_restarts(tmp_path):
    """A directory with only .bak / RESTART.wfn restart files should be
    treated as missing — the primary restart file is required."""
    (tmp_path / "cMD-1.restart.bak").write_text("")
    (tmp_path / "cMD-1-RESTART.wfn.restart").write_text("")
    (tmp_path / "cMD.LagrangeMultLog").write_text("")
    parser = CP2KParser()
    assert parser.is_constraint_directory(tmp_path) is False


# ---------------------------------------------------------------------------
# Registry + sniffer
# ---------------------------------------------------------------------------


def test_get_parser_returns_registered_cp2k():
    parser = get_parser("cp2k")
    assert isinstance(parser, CP2KParser)


def test_get_parser_case_insensitive():
    assert isinstance(get_parser("CP2K"), CP2KParser)
    assert isinstance(get_parser("Cp2k"), CP2KParser)


def test_get_parser_unknown_raises():
    with pytest.raises(ParserInferenceError):
        get_parser("definitely_not_a_real_engine")


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_infer_parser_finds_cp2k_on_real_data():
    parser = infer_parser(TI_POINT)
    assert isinstance(parser, CP2KParser)


def test_infer_parser_raises_on_unrecognised(tmp_path):
    with pytest.raises(ParserInferenceError):
        infer_parser(tmp_path)


def test_register_parser_then_infer(tmp_path):
    """Registering a fake parser whose marker file exists should let
    infer_parser pick it up."""

    class FakeParser:
        name = "fake_engine"

        def is_constraint_directory(self, directory):
            return (directory / "FAKE_MARKER").exists()

        def parse_metadata(self, directory):  # pragma: no cover
            raise NotImplementedError

        def parse_lambda_series(self, directory):  # pragma: no cover
            raise NotImplementedError

    (tmp_path / "FAKE_MARKER").write_text("")
    register_parser("fake_engine", FakeParser)
    try:
        # Empty tmp_path still doesn't have CP2K files; fake one wins
        # because cp2k.is_constraint_directory returns False here
        parser = infer_parser(tmp_path)
        assert parser.name == "fake_engine"
    finally:
        # Clean up registry to avoid bleed between tests
        from md_analysis.engines.protocols import _REGISTRY
        _REGISTRY.pop("fake_engine", None)


# ---------------------------------------------------------------------------
# resolve_parser
# ---------------------------------------------------------------------------


def test_resolve_parser_rejects_auto_sentinel():
    with pytest.raises(ValueError, match="auto"):
        resolve_parser("auto")


def test_resolve_parser_passes_through_instance():
    p = CP2KParser()
    assert resolve_parser(p) is p


def test_resolve_parser_resolves_by_name():
    assert isinstance(resolve_parser("cp2k"), CP2KParser)
