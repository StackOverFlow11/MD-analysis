"""Stable return type for workflow entry points.

``WorkflowResult`` is the contract every ``md_analysis.workflows.*.run_*``
function returns. The shape is deliberately small:

* ``artifacts`` maps short keys (e.g. ``"csv"``, ``"png"``) to concrete
  on-disk paths.
* ``metadata`` carries lightweight scalars (frame counts, methods,
  units, mode flags) that callers can inspect without re-reading files.
* ``extra`` is reserved for strongly-typed diagnostic reports (e.g.
  TI convergence reports) whose schema would be lossy if flattened
  into ``metadata``.

The :func:`require_artifacts_exist` helper is the canonical validation
point: workflow tests and CLI code call it after a run to fail loudly
when a promised artifact never materialised on disk.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..exceptions import MDAnalysisError


class MissingArtifactError(MDAnalysisError):
    """Raised when a :class:`WorkflowResult` lists artifacts that are
    not present on disk.

    The error message names every missing key/path pair so failures are
    self-describing in test output and CLI logs.
    """


@dataclass(frozen=True)
class WorkflowResult:
    """Stable return type for ``md_analysis.workflows.*.run_*`` entry points.

    Parameters
    ----------
    name
        Stable identifier of the workflow (e.g. ``"water_three_panel"``).
        Use the same value across re-runs so downstream consumers can
        match results by name.
    output_dir
        Directory under which all listed artifacts live. Paths in
        ``artifacts`` may be absolute or relative to ``output_dir``.
    artifacts
        Mapping of short keys to concrete file paths. Keys are workflow
        defined and documented per ``run_*`` function.
    metadata
        Lightweight scalars describing the run (frame counts, modes,
        physical units, ...). Must be JSON-serialisable in spirit.
    extra
        Optional strongly-typed report. If supplied, the object should
        expose ``to_dict()`` so :meth:`to_dict` can serialise it.
    """

    name: str
    output_dir: Path
    artifacts: dict[str, Path] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    extra: Any | None = None

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-friendly representation.

        ``extra`` is serialised via its own ``to_dict()`` when present;
        otherwise it is returned as-is and the caller is responsible
        for ensuring it is serialisable.
        """
        if self.extra is None:
            extra_serialised: Any = None
        elif hasattr(self.extra, "to_dict"):
            extra_serialised = self.extra.to_dict()
        else:
            extra_serialised = self.extra
        return {
            "name": self.name,
            "output_dir": str(self.output_dir),
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "metadata": dict(self.metadata),
            "extra": extra_serialised,
        }


def require_artifacts_exist(result: WorkflowResult) -> None:
    """Assert that every path in ``result.artifacts`` exists on disk.

    Resolves relative artifact paths against ``result.output_dir`` before
    checking. Raises :class:`MissingArtifactError` listing all missing
    ``(key, resolved_path)`` pairs in one message so callers can fix
    them in a single pass.
    """
    missing: list[tuple[str, Path]] = []
    for key, raw_path in result.artifacts.items():
        candidate = Path(raw_path)
        resolved = (
            candidate if candidate.is_absolute() else result.output_dir / candidate
        )
        if not resolved.exists():
            missing.append((key, resolved))
    if missing:
        listing = "\n".join(f"  - {k}: {p}" for k, p in missing)
        raise MissingArtifactError(
            f"WorkflowResult '{result.name}' is missing artifacts:\n{listing}"
        )
