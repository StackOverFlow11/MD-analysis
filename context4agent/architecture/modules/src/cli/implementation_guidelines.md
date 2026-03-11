# `md_analysis.cli` implementation guidelines

> Source: `src/md_analysis/cli/`

## Role

Interactive CLI package providing a VASPKIT-style numbered menu interface. Replaces the former `CLI.py` argparse module.

## Design principles

- Pure interactive: no command-line arguments, all input via `input()` prompts
- Each sub-menu module (`_water.py`, `_potential.py`, `_charge.py`, `_scripts.py`, `_settings.py`) is self-contained with its own menu display, parameter collection, and handler dispatch
- Prompt helpers in `_prompt.py` are shared across all sub-menus
- Handlers delegate to `main.py` workflow functions or directly to sub-package analysis functions
- `KeyboardInterrupt` / `EOFError` caught at top level for clean exit
- All `_cmd_*()` / `_run_*()` handlers decorated with `@_handle_cmd_error` (from `_prompt.py`): catches `MDAnalysisError`/`FileNotFoundError`/`ValueError`/`RuntimeError` with user-friendly message, and unexpected exceptions with type name; returns 1 on error

## Dependencies

- `cli` -> `main` (for integrated workflow functions like `run_water_analysis`)
- `cli` -> `water`, `potential`, `charge` (for individual analysis functions)
- `cli` -> `potential.config` (for default constants)
- `cli` -> `scripts` (for `generate_bader_workdir`, `batch_generate_bader_workdirs`)
- `cli` -> `utils.CellParser` (for `parse_abc_from_restart`, `parse_abc_from_md_inp`)
- `cli` -> `config` (for persistent user configuration, `CONFIGURABLE_DEFAULTS` registry, and `delete_config`)
- No reverse dependencies: no other module imports from `cli`

## Logging configuration

- `main()` configures a `StreamHandler` on the `"md_analysis"` logger at `INFO` level, outputting to stderr with format `"%(levelname)s: %(message)s"`
- Configuration is idempotent: only adds the handler if no non-NullHandler handlers exist
- The library itself uses `NullHandler` (set in `md_analysis/__init__.py`), so logging is silent unless the CLI (or an application) explicitly configures a handler
- `_handle_cmd_error` logs unexpected exceptions at `ERROR` level with `exc_info=True` for full traceback in logs, while printing a concise message to stdout for the user

## Handler naming convention

Most sub-menu handler functions follow a `_cmd_<code>()` naming pattern, where `<code>` is the menu code number (e.g., `_cmd_101()`, `_cmd_201()`). This applies across `_water.py`, `_potential.py`, `_scripts.py`, and `_settings.py`. Exception: `_charge.py` uses a shared `_run_charge(params)` handler (with `_collect_params()` for parameter gathering) instead of per-code `_cmd_301()/302()/303()` functions. All handlers are decorated with `@_handle_cmd_error` for unified error handling.

## Cell parameter acquisition

All sub-menus requiring cell parameters (water 101-105, scripts 401-402) use the shared `_prompt_cell_abc()` helper from `_prompt.py`:

1. Prompt cell source: `.restart` (default) or `md.inp`
2. Parse the chosen file (`parse_abc_from_restart` or `parse_abc_from_md_inp`)
3. On failure, offer one retry with a different file
4. On second failure, raise `CellParseError` (caught by `@_handle_cmd_error`)
5. Return `(a, b, c)` tuple, passed as `cell_abc` keyword argument to analysis functions

## Parameter flow

1. User selects analysis code from sub-menu
2. Required parameters prompted (if any)
3. "Modify advanced parameters? (y/N)" gate for optional parameters
4. Handler calls the appropriate analysis function
5. Results printed, program exits

## Configurable analysis defaults

Settings menu 903-907 allow users to persistently override algorithm defaults from `utils/config.py`:

- 903: `layer_tol_A` (layer clustering tolerance)
- 904: `z_bin_width_A` (z-axis bin width)
- 905: `theta_bin_deg` (theta bin width)
- 906: `water_oh_cutoff_A` (water O-H cutoff)
- 907: reset all analysis defaults

The `_get_effective_default(key)` helper in `_prompt.py` reads the user config first, falling back to the hardcoded default from `CONFIGURABLE_DEFAULTS` registry. Analysis sub-menus (`_potential.py`, `_water.py`, `_charge.py`) use this helper to populate prompt defaults. Library function signatures remain unchanged — persistence only affects CLI prompt defaults.
