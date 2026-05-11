"""CP2K stdout / SP-out parsing helpers.

Regex + tiny streaming parsers for ``md.out`` and ``sp.out`` files produced
by CP2K. Extracts the per-step Fermi energy (Hartree) plus the step number
and time embedded in the surrounding ``STEP NUMBER`` / ``TIME [fs]`` lines.

Phase 7a status
---------------
Plain mechanical move from ``electrochemical.potential._frame_source`` —
no behaviour change. The dict-of-records return shape of
``parse_md_out_fermi`` is preserved exactly; introducing a ``FermiRecord``
dataclass is deliberately deferred to Phase 7b together with the
``engines/cp2k.py`` facade rollout.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from .cube import _float

# Re-usable line-level regex; centralised here so cp2k_xyz / future
# engines/cp2k facade share the exact same patterns rather than redefining
# them per call site.
FERMI_RE = re.compile(
    r"Fermi energy:\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)"
)
STEP_RE = re.compile(r"STEP NUMBER\s*=\s*(\d+)")
TIME_RE = re.compile(
    r"TIME\s*\[fs\]\s*=\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)"
)


def parse_md_out_fermi(md_out_path: Path) -> list[dict]:
    """Parse ``(step, time_fs, fermi_raw)`` records from CP2K ``md.out``.

    Returns a list of dictionaries with keys ``step``, ``time_fs``,
    ``fermi_raw`` (Hartree). Records without an associated Fermi energy
    are dropped. ``time_fs`` may be ``None`` for steps where the
    ``TIME [fs]`` line was not present.

    Streaming implementation: one pass, ``O(file size)`` memory.
    """
    records: list[dict] = []
    fermi_pending_raw: Optional[float] = None
    last_rec: Optional[dict] = None

    with md_out_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = FERMI_RE.search(line)
            if m:
                fermi_pending_raw = _float(m.group(1))
                continue

            m = STEP_RE.search(line)
            if m:
                step = int(m.group(1))
                rec = {"step": step, "time_fs": None, "fermi_raw": None}
                if fermi_pending_raw is not None:
                    rec["fermi_raw"] = fermi_pending_raw
                    fermi_pending_raw = None
                records.append(rec)
                last_rec = rec
                continue

            m = TIME_RE.search(line)
            if m and last_rec is not None and last_rec["time_fs"] is None:
                last_rec["time_fs"] = _float(m.group(1))

    return [r for r in records if r["fermi_raw"] is not None]


def parse_sp_out_fermi(sp_out_path: Path) -> float | None:
    """Extract the last ``Fermi energy:`` value (Hartree) from a single-point
    ``sp.out`` file.

    Single-point CP2K runs emit no ``STEP NUMBER`` line, so the caller is
    responsible for sourcing the step from the surrounding filesystem
    layout (e.g. ``potential_t*_i*`` directory name).

    Returns ``None`` if the file does not contain any ``Fermi energy:``
    line.
    """
    fermi_raw: float | None = None
    with sp_out_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = FERMI_RE.search(line)
            if m:
                fermi_raw = _float(m.group(1))
    return fermi_raw


__all__ = [
    "FERMI_RE",
    "STEP_RE",
    "TIME_RE",
    "parse_md_out_fermi",
    "parse_sp_out_fermi",
]
