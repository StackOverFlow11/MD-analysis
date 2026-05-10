"""Engine-agnostic parser layer for constrained MD inputs.

Provides a :class:`ConstraintMDParser` Protocol that abstracts how a
constraint-MD point directory is parsed, plus a CP2K implementation,
a sniffer (:func:`infer_parser`), and a registry for future engines
(VASP, etc.).

Downstream TI / SG analysis code only depends on the Protocol surface;
adding a new engine = registering a new parser, no analysis-side
changes required.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Protocol, runtime_checkable

from ..exceptions import MDAnalysisError
from ..utils.RestartParser.ColvarParser import (
    ColvarRestart,
    LagrangeMultLog,
    parse_colvar_restart,
    parse_lagrange_mult_log,
)


class ParserInferenceError(MDAnalysisError):
    """Raised when no registered parser recognises a directory."""


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------

@runtime_checkable
class ConstraintMDParser(Protocol):
    """Engine-agnostic contract for reading a constraint-MD point.

    Implementations must be cheap to construct (no I/O in __init__).
    """

    name: str

    def is_constraint_directory(self, directory: Path) -> bool:
        """Return ``True`` if *directory* contains the files this parser
        recognises.  Used both as a sniffing primitive and as the default
        ``dir_filter`` in discovery.
        """
        ...

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        """Parse engine metadata: timestep, target value, growth rate, etc.

        Cheap operation (KB-sized file). Discovery layer uses this to
        sort points by ξ before reading the heavy λ(t) series.
        """
        ...

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
        """Parse the Lagrange-multiplier (constraint force) time series.

        Heavy operation. Called only when analysis actually needs the data.
        """
        ...


# ---------------------------------------------------------------------------
# CP2K implementation
# ---------------------------------------------------------------------------

class CP2KParser:
    """Parser for CP2K constraint-MD point directories.

    Recognises directories containing ``*.restart`` and
    ``*.LagrangeMultLog`` files (the standard CP2K output pair).
    """

    name = "cp2k"

    def is_constraint_directory(self, directory: Path) -> bool:
        if not directory.is_dir():
            return False
        try:
            self._find_restart(directory)
            self._find_log(directory)
        except FileNotFoundError:
            return False
        return True

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        return parse_colvar_restart(self._find_restart(directory))

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
        return parse_lagrange_mult_log(self._find_log(directory))

    # ------------------------------------------------------------------
    # File discovery (private)
    # ------------------------------------------------------------------

    @staticmethod
    def _find_restart(directory: Path) -> Path:
        """Find the primary .restart file (skips .bak and .RESTART.wfn)."""
        candidates = sorted(directory.glob("*.restart"))
        candidates = [
            p for p in candidates
            if ".bak" not in p.name and "RESTART.wfn" not in p.name
        ]
        if not candidates:
            raise FileNotFoundError(f"No .restart file in {directory}")
        # Prefer the one with highest suffix number (e.g. cMD-1_1500.restart)
        return candidates[-1]

    @staticmethod
    def _find_log(directory: Path) -> Path:
        """Find the .LagrangeMultLog file."""
        candidates = list(directory.glob("*.LagrangeMultLog"))
        if not candidates:
            raise FileNotFoundError(f"No .LagrangeMultLog file in {directory}")
        return candidates[0]


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


# ---------------------------------------------------------------------------
# Default registrations
# ---------------------------------------------------------------------------

register_parser("cp2k", CP2KParser)
