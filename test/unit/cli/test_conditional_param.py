"""Tests for ConditionalParam wrapper and Potential command sp_* params."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from md_analysis.cli._params import (
    K,
    ConditionalParam,
    StrParam,
    sp_cube_filename,
    sp_dir_pattern,
    sp_out_filename,
    sp_root_dir,
)
from md_analysis.cli._prompt import set_input_source


@pytest.fixture(autouse=True)
def _restore_input():
    yield
    set_input_source(input)


class TestConditionalParam:
    def test_collects_when_predicate_true(self):
        inner = StrParam("k", "label", default="fallback")
        cp = ConditionalParam(inner, lambda ctx: True)
        set_input_source(lambda _: "user_value")
        ctx: dict = {}
        cp.collect(ctx)
        assert ctx["k"] == "user_value"

    def test_applies_default_when_predicate_false(self):
        inner = StrParam("k", "label", default="fallback")
        cp = ConditionalParam(inner, lambda ctx: False)
        # No input source set — if prompted the test would fail
        set_input_source(lambda _: pytest.fail("should not prompt"))
        ctx: dict = {}
        cp.collect(ctx)
        assert ctx["k"] == "fallback"

    def test_apply_default_bypasses_predicate(self):
        """apply_default should always apply the inner default, not prompt."""
        inner = StrParam("k", "label", default="fallback")
        cp = ConditionalParam(inner, lambda ctx: True)
        set_input_source(lambda _: pytest.fail("should not prompt"))
        ctx: dict = {}
        cp.apply_default(ctx)
        assert ctx["k"] == "fallback"

    def test_predicate_reads_ctx(self):
        """Predicate should be re-evaluated against current ctx each collect call."""
        inner = StrParam("k", "label", default="fallback")
        cp = ConditionalParam(inner, lambda ctx: ctx.get("mode") == "on")

        set_input_source(lambda _: "prompted")
        ctx_on: dict = {"mode": "on"}
        cp.collect(ctx_on)
        assert ctx_on["k"] == "prompted"

        set_input_source(lambda _: pytest.fail("should not prompt"))
        ctx_off: dict = {"mode": "off"}
        cp.collect(ctx_off)
        assert ctx_off["k"] == "fallback"


class TestPotentialSpParams:
    """sp_* params should only prompt in distributed mode."""

    def test_continuous_mode_silent(self):
        """In continuous mode all four sp_* ConditionalParams apply defaults silently."""
        set_input_source(lambda _: pytest.fail("should not prompt"))
        ctx: dict = {K.INPUT_MODE: "continuous"}
        for p in (sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename):
            assert isinstance(p, ConditionalParam)
            p.collect(ctx)
        assert ctx[K.SP_ROOT_DIR] == "."
        assert ctx[K.SP_DIR_PATTERN] == "potential_t*_i*"
        assert ctx[K.SP_CUBE_FILENAME] == "sp_potential-v_hartree-1_0.cube"
        assert ctx[K.SP_OUT_FILENAME] == "sp.out"

    def test_distributed_mode_prompts(self):
        """In distributed mode all four sp_* ConditionalParams prompt the user."""
        responses = iter(["/data/sp", "custom_t*", "custom.cube", "custom.out"])
        set_input_source(lambda _: next(responses))
        ctx: dict = {K.INPUT_MODE: "distributed"}
        for p in (sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename):
            p.collect(ctx)
        assert ctx[K.SP_ROOT_DIR] == "/data/sp"
        assert ctx[K.SP_DIR_PATTERN] == "custom_t*"
        assert ctx[K.SP_CUBE_FILENAME] == "custom.cube"
        assert ctx[K.SP_OUT_FILENAME] == "custom.out"
