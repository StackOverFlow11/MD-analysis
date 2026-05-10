"""Structured Tools-layer contracts for agent tasks.

This module defines three dataclasses that together describe the **mechanical
contract** of a dispatchable task:

- ``FieldSpec``        — one input / output field (type, schema, description)
- ``ExceptionMapping`` — exception class ↔ dispatch error-type classification
- ``TaskContract``     — aggregate: inputs / outputs / preconditions /
                         side_effects / exceptions

Contracts are purely **mechanical**: they describe what a call does, returns,
writes, and raises. They do **not** encode business judgements like
"API-call-succeeded-but-result-is-bad" failure modes, which belong in a
future Resources layer.

See ``temp/agent_friendly_refactor/01_contract_design.md`` for the full
design rationale.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

ErrorType = Literal["validation", "file_not_found", "analysis", "internal"]
PathKind = Literal["file", "dir", "glob"]
OutputCategory = Literal["artifact", "metric", "raw_model"]


@dataclass(frozen=True)
class FieldSpec:
    """Describes a single input or output field.

    ``json_schema`` is the **authoritative, machine-readable** schema
    snippet (JSON Schema draft-07 compatible).  ``type`` is a human-readable
    hint only and is **not** parsed for schema generation or coercion.
    """

    description: str
    # Authoritative machine-readable schema (JSON Schema draft-07 snippet).
    # Example: {"type": "number"} / {"type": "array", "items": {"type": "number"}}
    json_schema: dict[str, Any]
    # Human-readable type hint (NOT used for schema/coercion).
    type: str = ""
    required: bool = True
    default: Any = None
    choices: tuple[str, ...] | None = None
    # Output-only: which MCP bucket this field belongs to.
    category: OutputCategory = "metric"
    # Domain annotations (purely descriptive; not consumed by dispatch).
    unit: str | None = None
    shape: str | None = None
    path_kind: PathKind | None = None


@dataclass(frozen=True)
class ExceptionMapping:
    """Maps one exception type to a dispatch error-type classification.

    ``exception_fqn`` is a fully-qualified path like
    ``"md_analysis.scripts.TIGen.TIGenError"`` or
    ``"builtins.FileNotFoundError"``, avoiding bare strings like
    ``"MDAnalysisError"`` that become ambiguous when multiple subclasses
    inherit from it.

    The ``exceptions`` tuple in :class:`TaskContract` is **ordered**:
    dispatch matches from front to back, so list more specific subclasses
    before their parents.
    """

    exception_fqn: str
    triggered_by: str
    error_type: ErrorType
    # Reserved for future MCP server layer.  Dispatch currently ignores this.
    user_visible: bool = True


@dataclass(frozen=True)
class TaskContract:
    """Tools-layer structured contract.

    One instance per task, written next to the handler and the
    ``register(TaskDef(...))`` call in ``_handlers.py``.
    """

    inputs: dict[str, FieldSpec]
    # Three-way output decomposition:
    #   artifacts  — files / dirs / globs (MCP tool returns URIs)
    #   metrics    — JSON-serializable scalars/dicts (MCP tool returns body)
    #   raw_model  — in-process Python objects (MCP does NOT return; exists
    #                so a future Resources layer can reference field paths).
    outputs_artifacts: dict[str, FieldSpec] = field(default_factory=dict)
    outputs_metrics: dict[str, FieldSpec] = field(default_factory=dict)
    outputs_raw_model: dict[str, FieldSpec] = field(default_factory=dict)
    # Input-level hard constraints (violations must raise).
    preconditions: tuple[str, ...] = ()
    # Files written / state mutated.
    side_effects: tuple[str, ...] = ()
    # Ordered: specific subclasses first, then parents.
    exceptions: tuple[ExceptionMapping, ...] = ()

    # ── Schema export (double-outlet; preserves public API shape) ────────

    def to_agent_schema(self, *, name: str, description: str) -> dict[str, Any]:
        """Produce OpenAI function-calling compatible schema.

        Shape: ``{"name", "description", "parameters": {type, properties,
        required}}``.  Used by ``get_task_schema()`` so the public API
        return shape stays identical to the pre-contract behaviour.
        """
        properties: dict[str, Any] = {}
        required: list[str] = []
        for pname, spec in self.inputs.items():
            # Copy json_schema; enrich with description / choices / default.
            prop: dict[str, Any] = dict(spec.json_schema)
            if spec.description and "description" not in prop:
                prop["description"] = spec.description
            if spec.choices and "enum" not in prop:
                prop["enum"] = list(spec.choices)
            if spec.default is not None and "default" not in prop:
                prop["default"] = spec.default
            properties[pname] = prop
            if spec.required:
                required.append(pname)
        return {
            "name": name,
            "description": description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        }

    def to_mcp_tool_schema(self, *, name: str, description: str) -> dict[str, Any]:
        """Produce MCP-style tool schema.

        Reserved for the future MCP server wrapper — **not** used by
        agent's public API today.
        """
        agent = self.to_agent_schema(name=name, description=description)
        return {
            "name": name,
            "description": description,
            "inputSchema": agent["parameters"],
        }
