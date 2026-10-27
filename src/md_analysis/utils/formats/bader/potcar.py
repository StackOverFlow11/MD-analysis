"""Parse VASP POTCAR files for ``(element, ZVAL)`` pairs."""

from __future__ import annotations

from pathlib import Path

from ._errors import BaderParseError


def _read_potcar_zval(path: Path) -> list[tuple[str, float]]:
    """Extract (element_symbol, zval) pairs from a VASP POTCAR file.

    Parameters
    ----------
    path : Path
        Path to the POTCAR file.

    Returns
    -------
    list of (str, float)
        Element symbol and valence electron count for each species block,
        in the order they appear in POTCAR (matching POSCAR element order).
    """
    lines = Path(path).read_text().splitlines()
    results: list[tuple[str, float]] = []

    element: str | None = None
    for line in lines:
        stripped = line.strip()
        # The first line of each element block starts with PAW_PBE/PAW_LDA/US:
        #   PAW_PBE Cu_pv 06Sep2000
        # Avoid matching TITEL lines (e.g. "TITEL  = PAW_PBE Cu ...").
        if stripped.startswith(("PAW_PBE", "PAW_LDA", "PAW_GGA", "US ")):
            token = stripped.split()[1]
            # Strip suffixes like _pv, _sv, _GW, etc.
            element = token.split("_")[0]

        if "ZVAL" in line and element is not None:
            # Line format: "   POMASS =   63.546; ZVAL   =   11.000    mass and valenz"
            after_zval = line.split("ZVAL")[1]
            # Extract the number after '='
            num_str = after_zval.split("=")[1].split()[0]
            zval = float(num_str)
            results.append((element, zval))
            element = None  # reset for next block

    if not results:
        raise BaderParseError(f"No ZVAL entries found in {path}")

    return results
