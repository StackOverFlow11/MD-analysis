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
- `cli` -> `config` (for persistent user configuration)
- No reverse dependencies: no other module imports from `cli`

## Parameter flow

1. User selects analysis code from sub-menu
2. Required parameters prompted (if any)
3. "Modify advanced parameters? (y/N)" gate for optional parameters
4. Handler calls the appropriate analysis function
5. Results printed, program exits
