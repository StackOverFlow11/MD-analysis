# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Produces water mass-density profiles, orientation-weighted profiles, adsorbed-layer detection, and integrated three-panel visualizations.

## Dependencies

```bash
pip install numpy matplotlib ase pytest
```

No `pyproject.toml` or `setup.py` — the package is used via `sys.path` manipulation in tests.

## Running Tests

```bash
# All tests
pytest test/

# Single unit test file
pytest test/unit/structure/utils/test_water_parser.py

# Integration test (end-to-end, generates PNG output)
python test/integration/structure/Analysis/test_water_three_panel_plot.py
# Output: test/_tmp_preview/water_three_panel_analysis.png
```

Integration tests are also runnable as standalone scripts (they add `REPO_ROOT` to `sys.path` themselves via `conftest.py`).

## Architecture

Two-layer design under `scripts/structure/`:

**Layer 1: `scripts.structure.utils` — Single-frame, low-level**
- Input: a single ASE `Atoms` object
- `LayerParser.py`: Detects metal layers, identifies the two interface layers (metal surfaces facing water), returns `Layer` dataclasses with `is_interface=True` and a `normal_unit` vector pointing metal→water
- `WaterParser.py`: Builds water topology (O-H bonds), computes density/orientation profiles, angle PDF
- `config.py`: Global constants (`TRANSITION_METAL_SYMBOLS`, `DEFAULT_Z_BIN_WIDTH_A=0.1`, `DEFAULT_THETA_BIN_DEG=5.0`, `DEFAULT_WATER_OH_CUTOFF_A=1.25`)

**Layer 2: `scripts.structure.Analysis` — Multi-frame workflows**
- Input: xyz trajectory path + CP2K `md.inp` path
- `WaterAnalysis/_common.py` (private): Trajectory I/O, per-frame layer/water detection, resampling to uniform normalized grid, ensemble averaging
- `WaterAnalysis/WaterDensity.py`, `WaterOrientation.py`, `AdWaterOrientation.py`: Individual analysis workflows
- `Water.py`: `plot_water_three_panel_analysis()` — the primary public entry point; reads trajectory, produces 5 CSVs + 1 TXT + 1 PNG

## Key Data Contracts

Physical units (enforced in `history/architecture/modules/data_contract.md`):
- Density: `g/cm³`
- Orientation-weighted density: `1/Å³`
- Angle: degrees
- Position axis: fractional cell coordinate (`c_fraction`) or Å depending on context

Water molecule indices returned as `(n_water, 3)` array with column order `[O_idx, H1_idx, H2_idx]`.

Cell parsing from `md.inp` expects the line: `ABC [angstrom] a b c`.

## System Assumptions

- Orthogonal simulation cell (no triclinic support)
- Always a metal/water interface → always exactly **two** interface layers
- Periodic boundary conditions apply; analysis uses MIC (Minimum Image Convention)
- Profiles computed from one interface to the cell midpoint, then ensemble-averaged

## Example Data

`data_example/potential/` contains a minimal CP2K output set:
- `md.inp` — cell parameters
- `md-pos-1.xyz` — trajectory frames

## Architecture Documentation

Detailed contracts, module interface specs, and implementation rules live in `history/`:
- `history/architecture/modules/data_contract.md` — output shapes, units, CSV headers
- `history/architecture/modules/glossary_units.md` — terminology and unit definitions
- `history/architecture/modules/scripts/**/interface_exposure.md` — per-module public API
- `history/architecture/modules/scripts/**/implementation_guidelines.md` — per-module rules
- `temp.md` — known pending issues (unit label consistency, API semantic clarifications)
