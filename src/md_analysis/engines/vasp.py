"""VASP engine adapter — placeholder.

The class implements the :class:`ConstraintMDParser` Protocol shape so
type checkers can see the surface, but every method raises
``NotImplementedError`` with a clear message.

Deliberately **not** auto-registered in ``engines/__init__.py``: VASP is
not yet a valid auto-discovery target, and silently substituting a stub
parser would mask the real "no engine recognised this directory" error.
Callers wanting to exercise the placeholder must construct it
explicitly:

    >>> from md_analysis.engines.vasp import VASPParser
    >>> parser = VASPParser()  # OK, no I/O
    >>> parser.parse_metadata(some_path)  # raises NotImplementedError
"""

from __future__ import annotations

from pathlib import Path

from .models import ConstraintMetadata, LambdaSeries


_NOT_IMPLEMENTED_MSG = (
    "VASPParser is a placeholder; VASP constraint-MD support is not "
    "implemented yet. See "
    "context4agent/requirements/overall_reconstruction_plan.md for the "
    "VASP extension notes."
)


class VASPParser:
    """Placeholder VASP constraint-MD parser.

    Implements the :class:`~md_analysis.engines.protocols.ConstraintMDParser`
    Protocol surface but every method raises ``NotImplementedError``.

    NOT registered in the default parser registry.
    """

    name = "vasp"

    def is_constraint_directory(self, directory: Path) -> bool:  # noqa: ARG002
        raise NotImplementedError(_NOT_IMPLEMENTED_MSG)

    def parse_metadata(self, directory: Path) -> ConstraintMetadata:  # noqa: ARG002
        raise NotImplementedError(_NOT_IMPLEMENTED_MSG)

    def parse_lambda_series(self, directory: Path) -> LambdaSeries:  # noqa: ARG002
        raise NotImplementedError(_NOT_IMPLEMENTED_MSG)


__all__ = ["VASPParser"]
