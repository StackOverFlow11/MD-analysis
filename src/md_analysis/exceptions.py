"""Base exception for the md_analysis package."""


class MDAnalysisError(Exception):
    """Base class for all md_analysis domain-specific errors.

    Callers can ``except MDAnalysisError`` to catch any domain error
    raised by the library.
    """
