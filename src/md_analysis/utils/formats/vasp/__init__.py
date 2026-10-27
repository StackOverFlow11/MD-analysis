"""VASP single-file parsers (placeholders).

Modules:

- ``report`` — VASP REPORT (constraint-MD metadata + lambda series); placeholder.
- ``outcar`` — VASP OUTCAR (Fermi series); placeholder.
- ``locpot`` — VASP LOCPOT (planar-averaged electrostatic potential); placeholder.

Each placeholder follows the explicit ``NotImplementedError`` pattern;
``engines.vasp.VASPParser`` is intentionally NOT registered in the
default parser registry, so ``infer_parser`` cannot dispatch to a stub.

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.formats.vasp.report import parse_vasp_report_metadata
"""

from __future__ import annotations

__all__: list[str] = []
