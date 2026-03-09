"""Reusable input-prompt helpers for the interactive CLI."""

from __future__ import annotations

import functools
import logging

from ..exceptions import MDAnalysisError

logger = logging.getLogger(__name__)


def _prompt_str(label: str, *, default: str | None = None) -> str | None:
    """Prompt for a string value. Returns *default* on empty input."""
    suffix = f" [{default}]" if default is not None else ""
    raw = input(f"  {label}{suffix}: ").strip()
    return raw if raw else default


def _prompt_str_required(label: str) -> str:
    """Prompt for a non-empty string; loops until valid."""
    while True:
        raw = input(f"  {label}: ").strip()
        if raw:
            return raw
        print("    (required, please enter a value)")


def _prompt_int(label: str, *, default: int | None = None) -> int | None:
    """Prompt for an integer. Returns *default* on empty input."""
    suffix = f" [{default}]" if default is not None else ""
    raw = input(f"  {label}{suffix}: ").strip()
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        print(f"    Invalid integer: {raw!r}, using default {default}")
        return default


def _prompt_float(label: str, *, default: float) -> float:
    """Prompt for a float. Returns *default* on empty input."""
    raw = input(f"  {label} [{default}]: ").strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"    Invalid number: {raw!r}, using default {default}")
        return default


def _prompt_choice(label: str, choices: list[str], *, default: str) -> str:
    """Prompt with a fixed set of valid choices."""
    opts = "/".join(choices)
    raw = input(f"  {label} ({opts}) [{default}]: ").strip()
    if not raw:
        return default
    if raw in choices:
        return raw
    print(f"    Invalid choice: {raw!r}, using default {default}")
    return default


def _prompt_bool(label: str, *, default: bool = True) -> bool:
    """Prompt for yes/no. Returns *default* on empty input."""
    hint = "Y/n" if default else "y/N"
    raw = input(f"  {label} ({hint}): ").strip().lower()
    if not raw:
        return default
    return raw in ("y", "yes")


def _parse_metal_elements(s: str | None) -> set[str] | None:
    """Parse comma-separated element symbols into a set, or None."""
    if s is None:
        return None
    elements = {e.strip() for e in s.split(",") if e.strip()}
    return elements or None


def _prompt_global_params() -> dict:
    """Prompt for shared global parameters (outdir, frame slicing)."""
    params: dict = {}
    params["outdir"] = _prompt_str("Output directory", default="analysis") or "analysis"
    params["frame_start"] = _prompt_int("Frame start (0-based)", default=None)
    params["frame_end"] = _prompt_int("Frame end (exclusive)", default=None)
    params["frame_step"] = _prompt_int("Frame step", default=None)
    return params


def _handle_cmd_error(func):
    """Catch analysis errors in CLI handlers, print message, return 1."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (MDAnalysisError, FileNotFoundError, ValueError, RuntimeError) as exc:
            print(f"\n  Error: {exc}")
            return 1
        except Exception as exc:
            logger.error("Unexpected error in %s: %s", func.__name__, exc, exc_info=True)
            print(f"\n  Unexpected error ({type(exc).__name__}): {exc}")
            return 1
    return wrapper
