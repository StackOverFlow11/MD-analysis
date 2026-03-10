"""Reusable input-prompt helpers for the interactive CLI."""

from __future__ import annotations

import functools
import logging

from ..exceptions import MDAnalysisError

logger = logging.getLogger(__name__)


def _get_effective_default(key: str):
    """Return user-configured value for *key*, falling back to hardcoded default."""
    from ..config import CONFIGURABLE_DEFAULTS, get_config

    entry = CONFIGURABLE_DEFAULTS[key]
    return get_config(key, entry["default"])


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


def _prompt_cell_abc() -> tuple[float, float, float]:
    """Prompt for cell parameters with restart-first priority and retry on failure.

    Interactive flow:
    1. Ask for cell source (.restart or md.inp), default .restart
    2. Parse the chosen file
    3. On failure, offer one retry with a different file
    4. On second failure, raise CellParseError (caught by @_handle_cmd_error)
    """
    from ..utils.RestartParser.CellParser import CellParseError, parse_abc_from_md_inp, parse_abc_from_restart

    for attempt in range(2):
        source = _prompt_choice("Cell source", [".restart", "md.inp"], default=".restart")

        try:
            if source == ".restart":
                path = _prompt_str_required("CP2K .restart file")
                abc = parse_abc_from_restart(path)
            else:
                path = _prompt_str_required("CP2K input file (e.g. md.inp)")
                abc = parse_abc_from_md_inp(path)
        except (CellParseError, FileNotFoundError) as exc:
            print(f"\n  Error: {exc}")
            if attempt == 0:
                if _prompt_bool("Retry with different file?", default=True):
                    continue
            raise CellParseError(str(exc)) from exc

        print(f"  Cell: a={abc[0]:.4f}, b={abc[1]:.4f}, c={abc[2]:.4f} A")
        return abc

    # Should not reach here, but satisfy type checker
    raise CellParseError("Failed to parse cell parameters.")  # pragma: no cover


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
