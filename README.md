# MD Analysis

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations.

## Features

**Water analysis** — density and orientation profiles along the interface-to-midpoint direction:

- Water mass-density profile (g/cm³)
- Orientation-weighted density profile (g/cm³)
- Adsorbed-layer auto-detection and orientation-angle distribution
- Integrated three-panel PNG visualization

**Potential analysis** — Hartree potential and electrode potential from cube files:

- Center-slab Hartree potential time series
- Fermi energy extraction from `md.out`
- Electrode potential U vs SHE (computational SHE)
- Plane-averaged φ(z) profile overlay across frames
- Thickness sensitivity sweep (dual-axis: U vs SHE + spatial std of φ(z))

**Charge analysis** — Bader surface charge density from VASP post-processing:

- Two surface charge methods:
  - `counterion` — only non-water, non-metal species (counterions, adsorbates) contribute to σ
  - `layer` — sum of interface-layer metal atom net charges / area
- Time-series CSV with cumulative average + PNG visualization
- Per-atom indexed charge extraction across trajectory

## Quick Start

### Install

```bash
pip install numpy matplotlib ase pytest
pip install .
```

### CLI (interactive menu)

Run `md-analysis` with no arguments to enter the interactive menu:

```bash
md-analysis
```

The VASPKIT-style numbered menu guides you through:

- **1) Water Analysis** — density (101), orientation (102), adsorbed-layer (103/104), full three-panel (105)
- **2) Potential Analysis** — center slab (201), Fermi (202), electrode U vs SHE (203), φ(z) (204), thickness sweep (205), full (206)
- **3) Charge Analysis** — counterion method (301), layer method (302), full with plots (303)

Each option prompts for required inputs, then offers an optional "Modify advanced parameters?" gate for output directory and frame slicing.

### Python API

```python
from md_analysis.water import plot_water_three_panel_analysis

plot_water_three_panel_analysis(
    xyz_path="data_example/potential/md-pos-1.xyz",
    md_inp_path="data_example/potential/md.inp",
    output_dir="output/",
)
```

Or use the programmatic entry points:

```python
from md_analysis.main import run_water_analysis, run_potential_analysis, run_charge_analysis, run_all
```

### Example input data

- `data_example/potential/md-pos-1.xyz` — trajectory frames
- `data_example/potential/md.inp` — cell parameters (`ABC [angstrom] a b c`)
- `data_example/potential/md.out` — CP2K output with Fermi energies
- `data_example/potential/md-POTENTIAL-v_hartree-1_*.cube` — Hartree potential cube snapshots
- `data_example/bader_work_dir/` — POSCAR, ACF.dat, POTCAR for Bader charge tests

## Project Layout

```
src/md_analysis/
├── cli/                    # VASPKIT-style interactive CLI (md-analysis console script)
│   ├── __init__.py         #   main() entry point, banner, top menu
│   ├── _prompt.py          #   reusable input prompt helpers
│   ├── _water.py           #   water sub-menu (101-105)
│   ├── _potential.py       #   potential sub-menu (201-206)
│   └── _charge.py          #   charge sub-menu (301-303)
├── main.py                 # programmatic entry points
├── utils/                  # single-frame low-level tools
│   ├── config.py           #   constants, unit conversions, cSHE parameters
│   ├── ClusterUtils.py     #   1D periodic clustering + gap detection
│   ├── CubeParser.py       #   cube file I/O, plane-averaged φ(z), slab-averaged potential
│   ├── LayerParser.py      #   metal layer detection, interface identification
│   └── WaterParser.py      #   water topology, density/orientation/angle profiles
├── water/                  # multi-frame water analysis workflows
│   ├── Water.py            #   plot_water_three_panel_analysis() — primary entry point
│   └── WaterAnalysis/      #   density, orientation, adsorbed-layer sub-workflows
│       ├── _common.py      #     trajectory I/O, per-frame computation, ensemble averaging
│       ├── WaterDensity.py
│       ├── WaterOrientation.py
│       └── AdWaterOrientation.py
├── potential/              # multi-frame potential analysis workflows
│   ├── config.py           #   output filename constants
│   ├── CenterPotential.py  #   center-slab, Fermi, electrode potential, thickness sweep
│   └── PhiZProfile.py      #   φ(z) overlay visualization
└── charge/                 # Bader charge analysis workflows
    ├── config.py           #   unit conversion constants, default filenames
    └── BaderAnalysis.py    #   surface charge (counterion/layer), indexed charges, trajectory analysis

test/
├── unit/
│   ├── utils/              # pytest unit tests (ClusterUtils, LayerParser, WaterParser)
│   └── charge/             # unit tests for BaderAnalysis
└── integration/
    ├── water/              # end-to-end water analysis scripts
    ├── potential/          # end-to-end potential analysis scripts
    └── charge/             # end-to-end charge analysis scripts

data_example/               # minimal reproducible input data
context4agent/                    # architecture contracts, decisions, requirements
```

## Running Tests

```bash
# All tests
pytest test/

# Single unit test file
pytest test/unit/utils/test_water_parser.py

# Integration tests (standalone scripts, require pip install first)
python test/integration/water/test_water_three_panel_plot.py
python test/integration/potential/test_center_potential.py
python test/integration/potential/test_phi_z_profile.py
```

## Architecture Notes

- **Three-layer design**: `utils` (single-frame primitives) → `water` / `potential` / `charge` (multi-frame workflows) → `main` / `cli` (integration entry points)
- **Interface detection**: 1D periodic clustering on metal z-coordinates, largest gap = water region, two bounding layers = interfaces
- **Water profiles**: computed from selected interface toward the cell midpoint, normalized to [0, 1], then ensemble-averaged across frames
- **Electrode potential**: U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP (computational SHE)
- **Surface charge**: two methods — `counterion` (non-water non-metal species only) and `layer` (interface-layer metal atom net charge sum)
- **Frame selection**: `frame_start` / `frame_end` / `frame_step` supported across all workflows (prompted via "advanced parameters" in the CLI)

## For Contributors

Architecture contracts and implementation rules:

- `context4agent/architecture/modules/data_contract.md` — output shapes, units, CSV headers
- `context4agent/architecture/modules/glossary_units.md` — terminology and unit definitions
- `context4agent/architecture/modules/src/**/interface_exposure.md` — public API per module
- `context4agent/architecture/modules/src/**/implementation_guidelines.md` — implementation rules
