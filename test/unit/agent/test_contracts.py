"""Tests for the Tools-layer contract infrastructure."""

from __future__ import annotations

import importlib

import pytest

from md_analysis.agent import get_task_schema, list_tasks
from md_analysis.agent._contracts import (
    ExceptionMapping,
    FieldSpec,
    TaskContract,
)
from md_analysis.agent._dispatch import _coerce_params_from_contract


# ---------------------------------------------------------------------------
# Schema generation — public API shape is preserved for contract tasks
# ---------------------------------------------------------------------------


class TestSchemaPublicShape:
    """Contract tasks must keep {name, description, parameters} shape."""

    def test_ti_gen_batch_schema_keys(self):
        s = get_task_schema("ti_gen_batch")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_ti_gen_batch_parameters_object(self):
        s = get_task_schema("ti_gen_batch")
        params = s["parameters"]
        assert params["type"] == "object"
        assert "properties" in params
        assert "required" in params

    def test_ti_gen_batch_includes_expected_fields(self):
        s = get_task_schema("ti_gen_batch")
        props = s["parameters"]["properties"]
        for key in (
            "inp_path", "xyz_path", "restart_path", "output_dir",
            "targets_au", "time_range", "steps", "script_path",
        ):
            assert key in props, f"missing field in schema: {key}"

    def test_time_range_is_object_not_tuple_array(self):
        """time_range must be an object schema (not tuple-form array)."""
        s = get_task_schema("ti_gen_batch")
        tr = s["parameters"]["properties"]["time_range"]
        assert tr.get("type") == "object"
        assert "properties" in tr
        assert set(tr["properties"].keys()) == {
            "time_initial_fs", "time_final_fs", "n_points",
        }

    def test_required_fields_match_contract(self):
        s = get_task_schema("ti_gen_batch")
        required = set(s["parameters"]["required"])
        # targets_au / time_range / steps / script_path are optional
        assert required == {"inp_path", "xyz_path", "restart_path", "output_dir"}


class TestSchemaNonContractTasksUnchanged:
    """Legacy (non-contract) tasks must retain the {name, description,
    parameters} public shape."""

    def test_water_three_panel_still_works(self):
        s = get_task_schema("water_three_panel")
        assert set(s.keys()) == {"name", "description", "parameters"}
        assert s["parameters"]["type"] == "object"


class TestListTasksIncludesTiGenBatch:
    def test_ti_gen_batch_registered(self):
        names = [t["name"] for t in list_tasks()]
        assert "ti_gen_batch" in names


class TestHandlerReloadRegression:
    """Batch 0 regression: after splitting _handlers.py into task modules,
    ``_reset_registry()`` + ``reload(_handlers)`` must still produce the
    full registry (currently 11 after the Phase 3 charge legacy cleanup
    removed charge_surface / charge_tracked / charge_counterion)."""

    def test_reset_then_reload_restores_full_registry(self):
        from md_analysis.agent import _handlers
        from md_analysis.agent._core import _reset_registry

        baseline = len(list_tasks())
        assert baseline == 11
        _reset_registry()
        assert len(list_tasks()) == 0
        importlib.reload(_handlers)
        assert len(list_tasks()) == baseline

    def test_reload_preserves_task_order(self):
        from md_analysis.agent import _handlers
        from md_analysis.agent._core import _reset_registry

        expected_order = [t["name"] for t in list_tasks()]
        _reset_registry()
        importlib.reload(_handlers)
        assert [t["name"] for t in list_tasks()] == expected_order


# ---------------------------------------------------------------------------
# ti_full_analysis contract surface
# ---------------------------------------------------------------------------


class TestTIFullAnalysisSchema:
    """Contract schema for ti_full_analysis must expose only the wrapper
    signature, not analyze_ti's internal parameters."""

    def test_schema_keys(self):
        s = get_task_schema("ti_full_analysis")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_includes_wrapper_params(self):
        s = get_task_schema("ti_full_analysis")
        props = s["parameters"]["properties"]
        expected = {
            "root_dir", "output_dir", "parser", "dir_filter", "reverse",
            "equilibration", "epsilon_tol_ev", "auto_equilibration",
            "point_slice",
        }
        assert expected <= set(props.keys())

    def test_excludes_internal_analyze_ti_params(self):
        s = get_task_schema("ti_full_analysis")
        props = s["parameters"]["properties"]
        for internal in ("xi_values", "lambda_series_list", "dt",
                         "time_starts", "engine_overrides"):
            assert internal not in props, f"leaked internal param: {internal}"

    def test_parser_field_exists(self):
        s = get_task_schema("ti_full_analysis")
        parser = s["parameters"]["properties"]["parser"]
        assert parser.get("type") == "string"

    def test_dir_filter_is_nullable_string(self):
        s = get_task_schema("ti_full_analysis")
        df = s["parameters"]["properties"]["dir_filter"]
        assert df.get("type") == ["string", "null"]

    def test_point_slice_is_nullable_string(self):
        s = get_task_schema("ti_full_analysis")
        ps = s["parameters"]["properties"]["point_slice"]
        # Either ["string", "null"] or equivalent nullable form
        assert ps.get("type") == ["string", "null"] or ps.get("type") == "string"

    def test_no_required_params(self):
        """All ti_full_analysis inputs have defaults → required list is empty."""
        s = get_task_schema("ti_full_analysis")
        assert s["parameters"]["required"] == []

    def test_registered(self):
        names = [t["name"] for t in list_tasks()]
        assert "ti_full_analysis" in names


# ---------------------------------------------------------------------------
# bader_gen_batch contract surface
# ---------------------------------------------------------------------------


class TestBaderGenBatchSchema:
    """Contract schema for bader_gen_batch — script preparation only,
    no submission / parsing."""

    def test_schema_keys(self):
        s = get_task_schema("bader_gen_batch")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_includes_all_wrapper_params(self):
        s = get_task_schema("bader_gen_batch")
        props = s["parameters"]["properties"]
        expected = {
            "xyz_path", "cell_abc", "output_dir",
            "mode", "frame_start", "frame_end", "frame_step",
            "time_start_fs", "time_end_fs", "time_step_fs",
            "script_path", "element_order",
            "generate_potcar", "direct", "verbose",
        }
        assert expected == set(props.keys())

    def test_required_params_only_three(self):
        """xyz_path / cell_abc / output_dir are required; rest have defaults."""
        s = get_task_schema("bader_gen_batch")
        assert set(s["parameters"]["required"]) == {
            "xyz_path", "cell_abc", "output_dir",
        }

    def test_mode_enum(self):
        s = get_task_schema("bader_gen_batch")
        m = s["parameters"]["properties"]["mode"]
        assert set(m.get("enum", [])) == {"index", "time"}

    def test_cell_abc_is_length_3_array(self):
        s = get_task_schema("bader_gen_batch")
        ca = s["parameters"]["properties"]["cell_abc"]
        assert ca["type"] == "array"
        assert ca["minItems"] == 3 and ca["maxItems"] == 3
        assert ca["items"]["type"] == "number"

    def test_registered_task_count(self):
        names = [t["name"] for t in list_tasks()]
        assert "bader_gen_batch" in names
        # 14 -> 11 after Phase 3 charge legacy cleanup removed
        # charge_surface / charge_tracked / charge_counterion.
        assert len(names) == 11


# ---------------------------------------------------------------------------
# TaskContract.to_agent_schema / to_mcp_tool_schema
# ---------------------------------------------------------------------------


def _mini_contract() -> TaskContract:
    return TaskContract(
        inputs={
            "a": FieldSpec(
                description="int field",
                json_schema={"type": "integer"},
            ),
            "b": FieldSpec(
                description="optional string",
                json_schema={"type": "string"},
                required=False,
                default="hello",
            ),
        },
    )


class TestTaskContractExports:
    def test_agent_schema_shape(self):
        c = _mini_contract()
        s = c.to_agent_schema(name="demo", description="desc")
        assert s["name"] == "demo"
        assert s["description"] == "desc"
        assert s["parameters"]["type"] == "object"
        assert s["parameters"]["properties"]["a"]["type"] == "integer"
        assert s["parameters"]["required"] == ["a"]
        # default surfaced
        assert s["parameters"]["properties"]["b"]["default"] == "hello"

    def test_mcp_tool_schema_shape(self):
        c = _mini_contract()
        s = c.to_mcp_tool_schema(name="demo", description="desc")
        assert s["name"] == "demo"
        assert s["description"] == "desc"
        assert "inputSchema" in s
        assert s["inputSchema"]["type"] == "object"
        # Agent schema's "parameters" becomes MCP's "inputSchema" with same body
        assert "parameters" not in s


# ---------------------------------------------------------------------------
# _coerce_params_from_contract rules
# ---------------------------------------------------------------------------


class TestCoercionRules:
    def test_integer_accepts_int(self):
        inputs = {"n": FieldSpec("n", json_schema={"type": "integer"})}
        out = _coerce_params_from_contract(inputs, {"n": 5})
        assert out["n"] == 5 and isinstance(out["n"], int)

    def test_integer_accepts_int_string(self):
        inputs = {"n": FieldSpec("n", json_schema={"type": "integer"})}
        out = _coerce_params_from_contract(inputs, {"n": "7"})
        assert out["n"] == 7 and isinstance(out["n"], int)

    def test_integer_accepts_integer_float(self):
        inputs = {"n": FieldSpec("n", json_schema={"type": "integer"})}
        out = _coerce_params_from_contract(inputs, {"n": 5.0})
        assert out["n"] == 5 and isinstance(out["n"], int)

    def test_integer_rejects_non_integer_float(self):
        inputs = {"n": FieldSpec("n", json_schema={"type": "integer"})}
        with pytest.raises(ValueError):
            _coerce_params_from_contract(inputs, {"n": 5.5})

    def test_integer_rejects_bool(self):
        """bool is a subclass of int in Python — we explicitly reject it."""
        inputs = {"n": FieldSpec("n", json_schema={"type": "integer"})}
        with pytest.raises(ValueError):
            _coerce_params_from_contract(inputs, {"n": True})

    def test_number_coerces_to_float(self):
        inputs = {"x": FieldSpec("x", json_schema={"type": "number"})}
        out = _coerce_params_from_contract(inputs, {"x": "3.14"})
        assert out["x"] == 3.14 and isinstance(out["x"], float)

    def test_boolean_keeps_bool(self):
        inputs = {"b": FieldSpec("b", json_schema={"type": "boolean"})}
        out = _coerce_params_from_contract(inputs, {"b": True})
        assert out["b"] is True

    def test_boolean_rejects_string(self):
        inputs = {"b": FieldSpec("b", json_schema={"type": "boolean"})}
        with pytest.raises(ValueError):
            _coerce_params_from_contract(inputs, {"b": "true"})

    def test_array_keeps_list(self):
        inputs = {"xs": FieldSpec(
            "xs", json_schema={"type": "array", "items": {"type": "number"}},
        )}
        out = _coerce_params_from_contract(inputs, {"xs": [1.0, 2.0]})
        assert out["xs"] == [1.0, 2.0] and isinstance(out["xs"], list)

    def test_object_keeps_dict(self):
        inputs = {"d": FieldSpec("d", json_schema={"type": "object"})}
        val = {"k": 1}
        out = _coerce_params_from_contract(inputs, {"d": val})
        assert out["d"] is val

    def test_path_kind_resolves_to_path(self, tmp_path):
        inputs = {"p": FieldSpec(
            "p", json_schema={"type": "string"}, path_kind="file",
        )}
        # Use a string that exists on disk so .resolve() is stable
        f = tmp_path / "foo.txt"
        f.write_text("x")
        out = _coerce_params_from_contract(inputs, {"p": str(f)})
        from pathlib import Path as _P
        assert isinstance(out["p"], _P)
        assert out["p"] == f.resolve()

    def test_path_kind_none_value_passthrough(self):
        inputs = {"p": FieldSpec(
            "p", json_schema={"type": ["string", "null"]}, path_kind="file",
        )}
        out = _coerce_params_from_contract(inputs, {"p": None})
        assert out["p"] is None

    def test_unknown_key_passthrough(self):
        inputs = {"a": FieldSpec("a", json_schema={"type": "integer"})}
        out = _coerce_params_from_contract(inputs, {"z": "unchanged"})
        assert out["z"] == "unchanged"


# ---------------------------------------------------------------------------
# ExceptionMapping — FQN resolution
# ---------------------------------------------------------------------------


class TestCoercionErrorsWrapped:
    """Contract coercion errors must surface as classified TaskResults,
    not leak out of dispatch()."""

    def test_integer_field_with_non_integer_float_returns_validation(self):
        from md_analysis.agent import dispatch

        # ti_gen_batch: steps is integer → 5.5 should fail coercion,
        # and the failure must be wrapped (no raise).
        result = dispatch("ti_gen_batch", {
            "inp_path": "missing.inp",
            "xyz_path": "missing.xyz",
            "restart_path": "missing.restart",
            "output_dir": "out",
            "targets_au": [1.0],
            "steps": 5.5,
        })
        assert not result.success
        assert result.error_type == "validation"

    def test_boolean_field_with_string_returns_validation(self):
        from md_analysis.agent import dispatch

        # ti_full_analysis: auto_equilibration is boolean → "true" (str)
        # should fail coercion, wrapped as validation.
        result = dispatch("ti_full_analysis", {
            "root_dir": "missing",
            "output_dir": "out",
            "auto_equilibration": "true",
        })
        assert not result.success
        assert result.error_type == "validation"

    def test_integer_field_with_bool_returns_validation(self):
        """bool-is-int Python trap: contract must explicitly reject bools
        for integer fields."""
        from md_analysis.agent import dispatch

        result = dispatch("ti_gen_batch", {
            "inp_path": "missing.inp",
            "xyz_path": "missing.xyz",
            "restart_path": "missing.restart",
            "output_dir": "out",
            "targets_au": [1.0],
            "steps": True,
        })
        assert not result.success
        assert result.error_type == "validation"


class TestConfigShowContract:
    """Batch 1: config_show is contract-backed, read-only, honours
    ``config_path`` (previously ignored by the legacy handler)."""

    def test_schema_keys(self):
        s = get_task_schema("config_show")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_no_required_params(self):
        s = get_task_schema("config_show")
        assert s["parameters"]["required"] == []

    def test_config_path_is_nullable_string(self):
        s = get_task_schema("config_show")
        cp = s["parameters"]["properties"]["config_path"]
        assert cp["type"] == ["string", "null"]

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["config_show"].contract is not None

    def test_dispatch_custom_config_path(self, tmp_path):
        import json
        from md_analysis.agent import dispatch

        cfg = tmp_path / "cfg.json"
        cfg.write_text(json.dumps({"my_key": "my_val"}), encoding="utf-8")

        r = dispatch("config_show", {"config_path": str(cfg)})
        assert r.success
        assert r.summary["config"] == {"my_key": "my_val"}

    def test_dispatch_malformed_json_is_validation(self, tmp_path):
        from md_analysis.agent import dispatch

        cfg = tmp_path / "bad.json"
        cfg.write_text("{not valid json", encoding="utf-8")

        r = dispatch("config_show", {"config_path": str(cfg)})
        assert not r.success
        assert r.error_type == "validation"

    def test_dispatch_default_path_returns_dict(self):
        from md_analysis.agent import dispatch
        r = dispatch("config_show", {})
        assert r.success
        assert isinstance(r.summary.get("config"), dict)


class TestExceptionFQNResolution:
    def test_contract_exception_fqns_importable(self):
        """Every contract's exception_fqn must resolve to a real class."""
        from md_analysis.agent._core import _TASK_REGISTRY
        from md_analysis.agent._dispatch import _resolve_exception_class

        for name, task in _TASK_REGISTRY.items():
            if task.contract is None:
                continue
            for em in task.contract.exceptions:
                cls = _resolve_exception_class(em.exception_fqn)
                assert isinstance(cls, type)
                assert issubclass(cls, BaseException), (
                    f"{name}: {em.exception_fqn} is not an exception class"
                )

    def test_error_types_valid(self):
        """All exception mappings use one of the four valid error_types."""
        from md_analysis.agent._core import _TASK_REGISTRY

        valid = {"validation", "file_not_found", "analysis", "internal"}
        for name, task in _TASK_REGISTRY.items():
            if task.contract is None:
                continue
            for em in task.contract.exceptions:
                assert em.error_type in valid, (
                    f"{name}: invalid error_type {em.error_type!r}"
                )

    def test_exception_ordering_specific_first(self):
        """Within a contract, specific subclasses must precede their parents."""
        from md_analysis.agent._core import _TASK_REGISTRY
        from md_analysis.agent._dispatch import _resolve_exception_class

        for name, task in _TASK_REGISTRY.items():
            if task.contract is None:
                continue
            classes = [
                _resolve_exception_class(em.exception_fqn)
                for em in task.contract.exceptions
            ]
            for i, a in enumerate(classes):
                for b in classes[i + 1:]:
                    assert not (issubclass(b, a) and b is not a), (
                        f"{name}: {b.__name__} (subclass of {a.__name__}) "
                        f"listed after its parent"
                    )
