# CLAUDE.md

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Produces water mass-density profiles, orientation-weighted profiles, adsorbed-layer detection, Hartree potential analysis, Fermi energy extraction, electrode potential computation, and integrated visualizations.

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

Key flags: `--thickness` (default 7 Å), `--thickness-end` (sweep limit, default 15 Å), `--center-mode` (interface|cell), `--fermi-unit` (au|ev), `--no-compute-u`, `--no-phi-z`.

## Architecture

Standard src-layout under `src/md_analysis/` with three sub-packages:

- **`utils`** — Single-frame tools: metal layer detection, water topology, cube I/O, density/orientation profiles, constants
- **`water`** — Multi-frame water workflows: ensemble-averaged density/orientation, adsorbed-layer detection, three-panel PNG
- **`potential`** — Multi-frame potential workflows: center slab potential, Fermi energy, U vs SHE, φ(z) overlay, thickness sensitivity
- **`main.py` / `CLI.py`** — Programmatic entry points and `md-analysis` console script

Physical units: density `g/cm³`, orientation-weighted density `g/cm³`, potential `eV`, electrode potential `V vs SHE`, position `Å`, angle `degrees`.

Detailed architecture, data contracts, public API specs, and implementation rules: see `context4agent/`.

## System Assumptions

- Orthogonal simulation cell (no triclinic support)
- Always a metal/water interface → always exactly **two** interface layers
- Periodic boundary conditions; analysis uses MIC (Minimum Image Convention)
- Profiles computed from one interface to the cell midpoint, then ensemble-averaged
- CP2K cube files: z is fastest-running index, reshape as `(nx, ny, nz)`

## Example Data

`data_example/potential/` — minimal CP2K output set (`md.inp`, `md-pos-1.xyz`, `md.out`, `md-POTENTIAL-v_hartree-1_*.cube`).
