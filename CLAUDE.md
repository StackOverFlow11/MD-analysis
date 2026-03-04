# CLAUDE.md

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Produces water mass-density profiles, orientation-weighted profiles, adsorbed-layer detection, Hartree potential analysis, Fermi energy extraction, electrode potential computation, Bader charge analysis, and integrated visualizations.

## Quick Start

```bash
pip install numpy matplotlib ase pytest
pip install .          # install package (required before running tests)
pytest test/           # run all tests
```

Integration tests are also runnable as standalone scripts: `python test/integration/water/test_water_three_panel_plot.py`

## CLI

```bash
md-analysis water --xyz md-pos-1.xyz --md-inp md.inp
md-analysis potential --cube-pattern "md-POTENTIAL-v_hartree-1_*.cube"
md-analysis all --xyz md-pos-1.xyz --md-inp md.inp
```

Key flags: `--thickness` (default 7.5 Å), `--thickness-end` (sweep limit, default 15 Å), `--center-mode` (interface|cell), `--fermi-unit` (au|ev), `--no-compute-u`, `--no-phi-z`.

## Architecture

Standard src-layout under `src/md_analysis/` with four sub-packages:

- **`utils`** — Single-frame tools: metal layer detection, water topology, cube I/O, density/orientation profiles, Bader charge parsing, constants
- **`water`** — Multi-frame water workflows: ensemble-averaged density/orientation, adsorbed-layer detection, three-panel PNG
- **`potential`** — Multi-frame potential workflows: center slab potential, Fermi energy, U vs SHE, φ(z) overlay, thickness sensitivity
- **`charge`** — Bader charge workflows: single-frame surface charge density, per-atom charge selection, multi-frame trajectory analysis with CSV output
- **`main.py` / `CLI.py`** — Programmatic entry points and `md-analysis` console script

Physical units: density `g/cm³`, orientation-weighted density `g/cm³`, potential `eV`, electrode potential `V vs SHE`, surface charge density `μC/cm²` (internal `e/Å²`), position `Å`, angle `degrees`.

Detailed architecture, data contracts, public API specs, and implementation rules: see `context4agent/`.

## System Assumptions

- Orthogonal simulation cell (no triclinic support)
- Always a metal/water interface → always exactly **two** interface layers
- Periodic boundary conditions; analysis uses MIC (Minimum Image Convention)
- Profiles computed from one interface to the cell midpoint, then ensemble-averaged
- CP2K cube files: z is fastest-running index, reshape as `(nx, ny, nz)`

## context4agent Synchronization (MANDATORY)

**Every code change that affects architecture, public API, data contracts, or units MUST be accompanied by corresponding updates to `context4agent/`.** This is non-negotiable — `context4agent/` is the single source of truth for cross-session collaboration.

Specifically, when you:
- **Add/remove/rename a module or sub-package** → update `context4agent/architecture/README.md` (structure & data flow) + create/update mirror docs under `context4agent/architecture/modules/src/<module>/`
- **Change public API** (new function, changed signature, removed export) → update the corresponding `interface_exposure.md`
- **Change implementation patterns** (dependency direction, data flow, internal conventions) → update the corresponding `implementation_guidelines.md`
- **Change output format** (CSV columns, file names, data shapes) → update `context4agent/architecture/modules/data_contract.md`
- **Change units or terminology** → update `context4agent/architecture/modules/glossary_units.md`
- **Change top-level exports** (`__init__.py` `__all__`) → update `context4agent/architecture/modules/src/interface_exposure.md` + `implementation_guidelines.md`
- **Change project status or capabilities** → update `context4agent/requirements/short_term.md`

Mirror structure rule: `context4agent/architecture/modules/src/` must mirror `src/md_analysis/` directory structure. Each mirror directory maintains:
- `interface_exposure.md` — public symbols, import paths, compatibility
- `implementation_guidelines.md` — responsibilities, dependencies, conventions

## Example Data

- `data_example/potential/` — minimal CP2K output set (`md.inp`, `md-pos-1.xyz`, `md.out`, `md-POTENTIAL-v_hartree-1_*.cube`)
- `data_example/bader_work_dir/` — VASP Bader analysis output (`POSCAR`, `ACF.dat`, `POTCAR`)
