"""Engine-agnostic protocols and parser registry.

Defines the :class:`ConstraintMDParser` Protocol abstracting how a
constraint-MD point directory is parsed, plus a global registry of
named parser factories with auto-sniffing (:func:`infer_parser`).

Downstream TI / SG analysis depends only on the Protocol surface —
adding a new engine is implementing the Protocol + calling
:func:`register_parser`, with no upper-layer change.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Protocol, runtime_checkable

from ..exceptions import MDAnalysisError
from .models import ConstraintMetadata, LambdaSeries


class ParserInferenceError(MDAnalysisError):
    """Raised when no registered parser recognises a directory."""


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------

@runtime_checkable
class ConstraintMDParser(Protocol):
    """Engine-agnostic contract for reading a constraint-MD point.

    Implementations must be cheap to construct (no I/O in ``__init__``).
    """

    name: str

    def is_constraint_directory(self, directory: Path) -> bool:
        """Return ``True`` if *directory* contains the files this parser
        recognises.  Used both as a sniffing primitive and as the default
        ``dir_filter`` in discovery.
        """
        ...

    def parse_metadata(self, directory: Path) -> ConstraintMetadata:
        """Parse engine metadata: timestep, target value, growth rate, etc.

        Cheap operation (KB-sized file). Discovery layer uses this to
        sort points by ξ before reading the heavy λ(t) series.
        """
        ...

    def parse_lambda_series(self, directory: Path) -> LambdaSeries:
        """Parse the Lagrange-multiplier (constraint force) time series.

        Heavy operation. Called only when analysis actually needs the data.
        """
        ...


# ---------------------------------------------------------------------------
# Registry + sniffer
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, Callable[[], ConstraintMDParser]] = {}


def register_parser(name: str, factory: Callable[[], ConstraintMDParser]) -> None:
    """Register a parser factory under *name* (case-insensitive).

    The factory is a zero-arg callable returning a parser instance.
    Idempotent: re-registering the same name overwrites.
    """
    _REGISTRY[name.lower()] = factory


def get_parser(name: str) -> ConstraintMDParser:
    """Look up a registered parser by name.

    Raises
    ------
    ParserInferenceError
        If *name* is not registered.
    """
    key = name.lower()
    if key not in _REGISTRY:
        raise ParserInferenceError(
            f"Unknown parser {name!r}. Registered: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[key]()


def infer_parser(directory: Path) -> ConstraintMDParser:
    """Sniff *directory* and return the first registered parser that
    recognises it.

    Parameters
    ----------
    directory : Path
        A candidate constraint-MD point directory (or its parent — the
        caller is responsible for picking a representative sample).

    Raises
    ------
    ParserInferenceError
        If no registered parser recognises the directory.
    """
    directory = Path(directory)
    for name, factory in _REGISTRY.items():
        parser = factory()
        if parser.is_constraint_directory(directory):
            return parser
    raise ParserInferenceError(
        f"No registered parser recognises {directory}. "
        f"Registered: {sorted(_REGISTRY)}"
    )


def resolve_parser(parser: ConstraintMDParser | str) -> ConstraintMDParser:
    """Resolve a parser argument from public API.

    Accepts either an instance (returned as-is) or a name string.
    The literal string ``"auto"`` is a sentinel handled by callers
    (they have a directory to sniff against); :func:`resolve_parser`
    rejects it so misuse surfaces clearly.
    """
    if isinstance(parser, str):
        if parser == "auto":
            raise ValueError(
                "'auto' must be resolved by the caller via infer_parser(); "
                "resolve_parser() handles named parsers only."
            )
        return get_parser(parser)
    return parser


__all__ = [
    "ConstraintMDParser",
    "ParserInferenceError",
    "register_parser",
    "get_parser",
    "infer_parser",
    "resolve_parser",
]
