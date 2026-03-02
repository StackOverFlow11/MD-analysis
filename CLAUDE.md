# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Produces water mass-density profiles, orientation-weighted profiles, adsorbed-layer detection, Hartree potential analysis, Fermi energy extraction, electrode potential computation, and integrated visualizations.

## Dependencies

```bash
pip install numpy matplotlib ase pytest
```

Install in development mode (required before running tests):

```bash
pip install -e .
```

## Running Tests

```bash
# All tests
pytest test/

# Single unit test file
pytest test/unit/utils/test_water_parser.py

# Integration test (end-to-end, generates PNG output)
python test/integration/water/test_water_three_panel_plot.py
# Output: test/_tmp_preview/water_three_panel_analysis.png

# Potential integration tests
python test/integration/potential/test_center_potential.py
python test/integration/potential/test_phi_z_profile.py
```

Integration tests are also runnable as standalone scripts (requires `pip install -e .` from repo root).

## CLI

```bash
md-analysis water --xyz md-pos-1.xyz --md-inp md.inp
md-analysis potential --cube-pattern "md-POTENTIAL-v_hartree-1_*.cube"
md-analysis all --xyz md-pos-1.xyz --md-inp md.inp
```

## Architecture

Three-package design under `src/`:

**`src.utils` — Shared low-level tools (single-frame)**
- Input: a single ASE `Atoms` object, cube file, or raw arrays
- `config.py`: Global constants (`TRANSITION_METAL_SYMBOLS`, `HA_TO_EV`, `BOHR_TO_ANG`, cSHE constants, bin widths)
- `ClusterUtils.py`: 1D periodic clustering + largest-gap detection
- `CubeParser.py`: Gaussian cube file I/O, plane-averaged φ(z), slab-averaged potential
- `LayerParser.py`: Metal layer detection, interface identification (uses ClusterUtils)
- `WaterParser.py`: Water topology (O-H bonds), density/orientation profiles, angle PDF

**`src.water` — Water analysis workflows (multi-frame)**
- Input: xyz trajectory path + CP2K `md.inp` path
- `WaterAnalysis/_common.py` (private): Trajectory I/O, per-frame layer/water detection, ensemble averaging
- `WaterAnalysis/WaterDensity.py`, `WaterOrientation.py`, `AdWaterOrientation.py`: Individual workflows
- `Water.py`: `plot_water_three_panel_analysis()` — primary water entry point; produces 5 CSVs + 1 TXT + 1 PNG

**`src.potential` — Potential analysis workflows (multi-frame)**
- Input: cube file glob pattern + CP2K `md.out` + optional xyz trajectory
- `CenterPotential.py`: `center_slab_potential_analysis()`, `fermi_energy_analysis()`, `electrode_potential_analysis()`
- `PhiZProfile.py`: `phi_z_planeavg_analysis()` — φ(z) heatmap + overlay visualization

**`src.main` / `src.CLI` — Integration entry points**
- `main.py`: `run_water_analysis()`, `run_potential_analysis()`, `run_all()`
- `CLI.py`: argparse-based CLI registered as `md-analysis` console script

## Key Data Contracts

Physical units (enforced in `history/architecture/modules/data_contract.md`):
- Density: `g/cm³`
- Orientation-weighted density: `1/Å³`
- Potential: `eV` (Hartree converted via `HA_TO_EV`)
- Position: Å or fractional cell coordinate depending on context
- Angle: degrees

Water molecule indices returned as `(n_water, 3)` array with column order `[O_idx, H1_idx, H2_idx]`.

Cell parsing from `md.inp` expects the line: `ABC [angstrom] a b c`.

## System Assumptions

- Orthogonal simulation cell (no triclinic support)
- Always a metal/water interface → always exactly **two** interface layers
- Periodic boundary conditions apply; analysis uses MIC (Minimum Image Convention)
- Profiles computed from one interface to the cell midpoint, then ensemble-averaged
- CP2K cube files: z is fastest-running index, reshape as `(nx, ny, nz)`

## Example Data

`data_example/potential/` contains a minimal CP2K output set:
- `md.inp` — cell parameters
- `md-pos-1.xyz` — trajectory frames
- `md.out` — CP2K output with Fermi energies
- `md-POTENTIAL-v_hartree-1_*.cube` — Hartree potential cube snapshots

## Architecture Documentation

Detailed contracts, module interface specs, and implementation rules live in `history/`:
- `history/architecture/modules/data_contract.md` — output shapes, units, CSV headers
- `history/architecture/modules/glossary_units.md` — terminology and unit definitions
- `history/architecture/modules/src/**/interface_exposure.md` — per-module public API
- `history/architecture/modules/src/**/implementation_guidelines.md` — per-module rules
