# CLAUDE.md

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations.

## Quick Start

```bash
pip install numpy matplotlib ase pytest
pip install .          # install package (required before running tests)
pytest test/           # run all tests
```

## Architecture

Standard src-layout: `src/md_analysis/` with sub-packages `utils`, `water`, `potential`, `charge`, plus `main.py` / `CLI.py`.

Details → `context4agent/architecture/README.md`

## context4agent — Single Source of Truth (MANDATORY)

`context4agent/` is the authoritative reference for architecture, API contracts, units, and project status. **Every code change MUST synchronize the corresponding `context4agent/` docs.**

| What changed | Update where |
|---|---|
| Module / sub-package added/removed/renamed | `architecture/README.md` + `architecture/modules/src/<module>/` mirror docs |
| Public API (function, signature, export) | corresponding `interface_exposure.md` |
| Implementation pattern / dependency | corresponding `implementation_guidelines.md` |
| Output format (CSV columns, filenames, shapes) | `architecture/modules/data_contract.md` |
| Units / terminology | `architecture/modules/glossary_units.md` |
| Top-level `__all__` | `architecture/modules/src/interface_exposure.md` + `implementation_guidelines.md` |
| Project capabilities / status | `requirements/short_term.md` |

Mirror rule: `context4agent/architecture/modules/src/` mirrors `src/md_analysis/` — each directory has `interface_exposure.md` + `implementation_guidelines.md`.

Key entry points inside `context4agent/`:
- `architecture/README.md` — module structure, data flow, system assumptions
- `architecture/modules/data_contract.md` — cross-module data shapes, units, CSV specs
- `architecture/modules/glossary_units.md` — terminology and unit conventions
- `requirements/short_term.md` — current capabilities and near-term tasks
