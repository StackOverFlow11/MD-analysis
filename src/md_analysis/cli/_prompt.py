"""Low-level prompt functions for the interactive CLI."""

from __future__ import annotations

from typing import Callable

_input_fn: Callable[[str], str] = input


def set_input_source(fn: Callable[[str], str]) -> None:
    """Replace the input function (for testing)."""
    global _input_fn
    _input_fn = fn


def _read(prompt_text: str) -> str:
    """All prompt functions call this instead of input() directly."""
    return _input_fn(prompt_text)


def prompt_str(label: str, *, default: str | None = None) -> str | None:
    """Prompt for a string value. Returns *default* on empty input."""
    suffix = f" [{default}]" if default is not None else ""
    raw = _read(f"  {label}{suffix}: ").strip()
    return raw if raw else default


def prompt_str_required(label: str) -> str:
    """Prompt for a non-empty string; loops until valid."""
    while True:
        raw = _read(f"  {label}: ").strip()
        if raw:
            return raw
        print("    (required, please enter a value)")


def prompt_int(label: str, *, default: int | None = None) -> int | None:
    """Prompt for an integer. Returns *default* on empty input."""
    suffix = f" [{default}]" if default is not None else ""
    raw = _read(f"  {label}{suffix}: ").strip()
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        print(f"    Invalid integer: {raw!r}, using default {default}")
        return default


def prompt_float(label: str, *, default: float) -> float:
    """Prompt for a float. Returns *default* on empty input."""
    raw = _read(f"  {label} [{default}]: ").strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"    Invalid number: {raw!r}, using default {default}")
        return default


def prompt_choice(label: str, choices: list[str], *, default: str) -> str:
    """Prompt with a fixed set of valid choices."""
    opts = "/".join(choices)
    raw = _read(f"  {label} ({opts}) [{default}]: ").strip()
    if not raw:
        return default
    if raw in choices:
        return raw
    print(f"    Invalid choice: {raw!r}, using default {default}")
    return default


def prompt_bool(label: str, *, default: bool = True) -> bool:
    """Prompt for yes/no. Returns *default* on empty input."""
    hint = "Y/n" if default else "y/N"
    raw = _read(f"  {label} ({hint}): ").strip().lower()
    if not raw:
        return default
    return raw in ("y", "yes")
