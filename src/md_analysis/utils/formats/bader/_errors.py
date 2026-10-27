"""Shared exception class for Bader file parsing."""

from __future__ import annotations

from ....exceptions import MDAnalysisError


class BaderParseError(MDAnalysisError):
    """Raised when an ACF/POTCAR file format is invalid or atom counts mismatch."""
