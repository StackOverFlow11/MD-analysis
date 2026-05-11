"""I/O scaffolding for `md_analysis`.

Path discovery, generic file writers, and dispatching helpers that pick
which file to read. Modules here:

- Walk directories, glob patterns, and sort frame directories
- Write generic outputs (CSV from dicts / numpy arrays)
- Dispatch to format-specific parsers in ``md_analysis.utils.formats``
- Do NOT parse the contents of specific file formats themselves
  (those belong to ``md_analysis.utils.formats``)
- Do NOT carry analysis logic (those belong to the upper modules)

Consumers should import directly from the relevant submodule, e.g.::

    from md_analysis.utils.io._frame_discovery import discover_frame_dirs
    from md_analysis.utils.io._io_helpers import _write_csv_from_arrays
    from md_analysis.utils.io.cell_resolver import resolve_cell_abc
"""

from __future__ import annotations

__all__: list[str] = []
