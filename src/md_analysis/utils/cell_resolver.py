"""Non-interactive cell parameter resolution.

Provides :func:`resolve_cell_abc` — a pure-function wrapper that resolves
cell dimensions from multiple sources (direct value, .restart file, md.inp
file, or auto-discovery in a work directory).

This module complements the interactive ``CellAbcParam.collect()`` in
``cli/_params.py``, which adds retry logic and user prompts.  Both share
the same underlying parsers from ``RestartParser.CellParser``.
"""

from __future__ import annotations

import re
from pathlib import Path


def resolve_cell_abc(
    cell_abc: tuple[float, float, float] | list[float] | None = None,
    restart_path: str | Path | None = None,
    md_inp_path: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> tuple[float, float, float]:
    """Resolve cell dimensions by priority.

    Priority: direct value > .restart > md.inp > auto-discover in work_dir.

    Parameters
    ----------
    cell_abc : tuple or list, optional
        Direct cell lengths ``(a, b, c)`` in angstrom.
    restart_path : str or Path, optional
        Path to a CP2K ``.restart`` file.
    md_inp_path : str or Path, optional
        Path to a CP2K ``md.inp`` file.
    work_dir : str or Path, optional
        Directory to search for ``.restart`` or ``md.inp`` files.

    Returns
    -------
    tuple[float, float, float]

    Raises
    ------
    ValueError
        If no source can determine cell_abc.
    FileNotFoundError
        If *work_dir* contains no parseable cell file.
    """
    if cell_abc is not None:
        abc = tuple(float(x) for x in cell_abc)
        if len(abc) != 3:
            raise ValueError(f"cell_abc must have 3 elements, got {len(abc)}")
        return abc  # type: ignore[return-value]

    from .RestartParser.CellParser import parse_abc_from_md_inp, parse_abc_from_restart

    if restart_path is not None:
        return parse_abc_from_restart(Path(restart_path))

    if md_inp_path is not None:
        return parse_abc_from_md_inp(Path(md_inp_path))

    if work_dir is not None:
        return _auto_discover(Path(work_dir))

    raise ValueError(
        "Cannot determine cell_abc. Provide one of: "
        "cell_abc, restart_path, md_inp_path, or work_dir"
    )


def _auto_discover(work_dir: Path) -> tuple[float, float, float]:
    """Search *work_dir* for .restart or md.inp files."""
    from .RestartParser.CellParser import parse_abc_from_md_inp, parse_abc_from_restart

    # Prefer .restart (exclude _N.restart checkpoint files)
    for f in sorted(work_dir.glob("*.restart")):
        if not re.search(r"_\d+\.restart$", f.name):
            return parse_abc_from_restart(f)

    # Fallback to common input file names
    for name in ("md.inp", "cp2k.inp"):
        inp = work_dir / name
        if inp.exists():
            return parse_abc_from_md_inp(inp)

    raise FileNotFoundError(f"No .restart or md.inp found in {work_dir}")
