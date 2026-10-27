"""Unit tests for ``md_analysis.workflows.models``.

Covers ``WorkflowResult`` shape and serialisation, plus the
``require_artifacts_exist`` validation helper.
"""

from __future__ import annotations

from dataclasses import dataclass, is_dataclass
from pathlib import Path

import pytest

from md_analysis.workflows import (
    MissingArtifactError,
    WorkflowResult,
    require_artifacts_exist,
)
from md_analysis.workflows.models import MissingArtifactError as ModelMissing


# ---------------------------------------------------------------------------
# WorkflowResult shape
# ---------------------------------------------------------------------------


def test_workflow_result_is_frozen_dataclass() -> None:
    assert is_dataclass(WorkflowResult)
    result = WorkflowResult(name="x", output_dir=Path("/tmp"))
    with pytest.raises(Exception):
        # frozen=True → assignment raises FrozenInstanceError
        result.name = "y"  # type: ignore[misc]


def test_workflow_result_defaults() -> None:
    result = WorkflowResult(name="demo", output_dir=Path("/tmp/demo"))
    assert result.name == "demo"
    assert result.output_dir == Path("/tmp/demo")
    assert result.artifacts == {}
    assert result.metadata == {}
    assert result.extra is None


def test_workflow_result_defaults_are_independent_instances() -> None:
    """Mutable defaults must not be shared between instances."""
    a = WorkflowResult(name="a", output_dir=Path("/tmp/a"))
    b = WorkflowResult(name="b", output_dir=Path("/tmp/b"))
    a.artifacts["k"] = Path("v")
    assert b.artifacts == {}


def test_facade_reexports_match_module() -> None:
    """The package facade re-exports the same symbol as the models module."""
    assert MissingArtifactError is ModelMissing


# ---------------------------------------------------------------------------
# to_dict serialisation
# ---------------------------------------------------------------------------


def test_to_dict_serialises_paths_and_metadata(tmp_path: Path) -> None:
    csv_path = tmp_path / "out.csv"
    result = WorkflowResult(
        name="water",
        output_dir=tmp_path,
        artifacts={"csv": csv_path},
        metadata={"n_frames": 42, "method": "counterion"},
    )
    dumped = result.to_dict()
    assert dumped["name"] == "water"
    assert dumped["output_dir"] == str(tmp_path)
    assert dumped["artifacts"] == {"csv": str(csv_path)}
    assert dumped["metadata"] == {"n_frames": 42, "method": "counterion"}
    assert dumped["extra"] is None


def test_to_dict_calls_extra_to_dict_when_available() -> None:
    @dataclass(frozen=True)
    class _Report:
        score: float

        def to_dict(self) -> dict[str, float]:
            return {"score": self.score}

    result = WorkflowResult(
        name="ti", output_dir=Path("/tmp/ti"), extra=_Report(score=0.5)
    )
    dumped = result.to_dict()
    assert dumped["extra"] == {"score": 0.5}


def test_to_dict_passes_through_extra_without_to_dict() -> None:
    result = WorkflowResult(
        name="composite", output_dir=Path("/tmp/c"), extra={"raw": 1}
    )
    dumped = result.to_dict()
    assert dumped["extra"] == {"raw": 1}


# ---------------------------------------------------------------------------
# require_artifacts_exist
# ---------------------------------------------------------------------------


def test_require_artifacts_exist_passes_for_existing_absolute_paths(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "out.csv"
    csv_path.write_text("a,b\n1,2\n")
    result = WorkflowResult(
        name="demo", output_dir=tmp_path, artifacts={"csv": csv_path}
    )
    require_artifacts_exist(result)  # must not raise


def test_require_artifacts_exist_resolves_relative_paths_against_output_dir(
    tmp_path: Path,
) -> None:
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "plot.png").write_bytes(b"")
    result = WorkflowResult(
        name="demo",
        output_dir=sub,
        artifacts={"png": Path("plot.png")},  # relative
    )
    require_artifacts_exist(result)


def test_require_artifacts_exist_raises_with_all_missing_keys(
    tmp_path: Path,
) -> None:
    result = WorkflowResult(
        name="demo",
        output_dir=tmp_path,
        artifacts={
            "csv": tmp_path / "missing.csv",
            "png": tmp_path / "missing.png",
        },
    )
    with pytest.raises(MissingArtifactError) as excinfo:
        require_artifacts_exist(result)
    message = str(excinfo.value)
    assert "csv" in message
    assert "png" in message
    assert "demo" in message


def test_require_artifacts_exist_empty_artifacts_passes(tmp_path: Path) -> None:
    """A workflow may legitimately produce no artifacts (e.g. a dry-run)."""
    result = WorkflowResult(name="dry", output_dir=tmp_path)
    require_artifacts_exist(result)


def test_missing_artifact_error_is_md_analysis_error() -> None:
    from md_analysis.exceptions import MDAnalysisError

    assert issubclass(MissingArtifactError, MDAnalysisError)
