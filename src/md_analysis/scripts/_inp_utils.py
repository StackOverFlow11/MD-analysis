"""Shared CP2K input-file text manipulation helpers.

Used by PotentialGen and SpGen to prepare single-point calculation inputs
from an MD trajectory: update ``&CELL ABC`` and ensure ``&TOPOLOGY`` points
to ``init.xyz``.
"""

from __future__ import annotations

import logging
import re

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Regex definitions
# ---------------------------------------------------------------------------

_CELL_BLOCK_RE = re.compile(
    r"(&CELL\b)(.*?)(&END\s+CELL)",
    re.DOTALL | re.IGNORECASE,
)
_ABC_LINE_RE = re.compile(
    r"^(\s*ABC\s+)(\[.*?\]\s+)?\S+\s+\S+\s+\S+",
    re.MULTILINE | re.IGNORECASE,
)
_TOPOLOGY_BLOCK_RE = re.compile(
    r"(&TOPOLOGY\b)(.*?)(&END\s+TOPOLOGY)",
    re.DOTALL | re.IGNORECASE,
)
_COORD_FILE_NAME_RE = re.compile(
    r"^(\s*COORD_FILE_NAME)\s+\S+", re.MULTILINE | re.IGNORECASE,
)
_COORD_FILE_FORMAT_RE = re.compile(
    r"^(\s*COORD_FILE_FORMAT)\s+\S+", re.MULTILINE | re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def replace_cell_abc(inp_text: str, a: float, b: float, c: float) -> str:
    """Replace the ``ABC`` line inside ``&CELL`` with new values."""

    def _replace_in_cell(match: re.Match) -> str:
        header, body, footer = match.group(1), match.group(2), match.group(3)
        new_body = _ABC_LINE_RE.sub(
            rf"\g<1>[angstrom] {a:.4f}   {b:.4f}   {c:.4f}",
            body,
        )
        return header + new_body + footer

    result = _CELL_BLOCK_RE.sub(_replace_in_cell, inp_text)
    if result == inp_text:
        logger.warning("No &CELL ABC line found in template; cell not updated")
    return result


def ensure_topology_init_xyz(inp_text: str) -> str:
    """Ensure ``&TOPOLOGY`` has ``COORD_FILE_NAME init.xyz`` and ``COORD_FILE_FORMAT XYZ``."""

    def _replace_topology(match: re.Match) -> str:
        header, body, footer = match.group(1), match.group(2), match.group(3)

        if _COORD_FILE_NAME_RE.search(body):
            body = _COORD_FILE_NAME_RE.sub(r"\1 init.xyz", body)
        else:
            body += "      COORD_FILE_NAME init.xyz\n"

        if _COORD_FILE_FORMAT_RE.search(body):
            body = _COORD_FILE_FORMAT_RE.sub(r"\1 XYZ", body)
        else:
            body += "      COORD_FILE_FORMAT XYZ\n"

        return header + body + footer

    return _TOPOLOGY_BLOCK_RE.sub(_replace_topology, inp_text)


def modify_inp_for_sp(inp_text: str, cell_abc: tuple[float, float, float]) -> str:
    """Return *inp_text* modified for single-point calculation.

    Modifications:
    - ``&CELL ABC`` → updated cell parameters
    - ``&TOPOLOGY``: ensure ``COORD_FILE_NAME init.xyz`` + ``COORD_FILE_FORMAT XYZ``
    """
    text = replace_cell_abc(inp_text, *cell_abc)
    text = ensure_topology_init_xyz(text)
    return text
