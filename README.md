# MD Analysis

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations.

Core outputs:

- Water mass-density profile along interface-to-midpoint distance
- Orientation-weighted density profile ρ·⟨cosφ⟩ (g/cm³) along the same axis
- Adsorbed-layer auto-detection and orientation-angle distribution
- Integrated three-panel PNG visualization

## Quick Start

### 1) Set up environment

```bash
conda activate env_md_an   # numpy, matplotlib, ase, pytest included
```

Or install manually:

```bash
pip install numpy matplotlib ase pytest
```

### 2) Run end-to-end plot generation

From the repo root:

```bash
python test/integration/structure/Analysis/test_water_three_panel_plot.py
```

Outputs written to `test/_tmp_preview/`:

- `water_three_panel_analysis.png` — three-panel figure
- `water_mass_density_z_distribution_analysis.csv`
- `water_orientation_weighted_density_z_distribution_analysis.csv`
- `adsorbed_water_orientation_profile.csv`
- `adsorbed_water_layer_range.txt`
- `adsorbed_water_theta_distribution_0_180.csv`

### 3) Example input data

- `data_example/potential/md-pos-1.xyz` — trajectory frames
- `data_example/potential/md.inp` — cell parameters (`ABC [angstrom] a b c`)

## Recommended Entry Point

```python
from src.structure.Analysis import plot_water_three_panel_analysis

plot_water_three_panel_analysis(
    xyz_path="data_example/potential/md-pos-1.xyz",
    md_inp_path="data_example/potential/md.inp",
    output_dir="output/",
)
```

This single call produces:

1. Panel 1 — water mass density (g/cm³)
2. Panel 2 — orientation-weighted density ρ·⟨cosφ⟩ (g/cm³)
3. Panel 3 — adsorbed-layer θ probability distribution (degree⁻¹)

## Project Layout

```
src/structure/utils/     # single-frame, low-level (LayerParser, WaterParser, config)
src/structure/Analysis/  # multi-frame workflows and plot composition
test/unit/               # unit tests
test/integration/        # end-to-end runnable scripts
data_example/            # minimal reproducible input data
history/                 # architecture contracts, decisions, requirements
```

## Running Tests

```bash
# All tests
python -m pytest test/ -v

# Single unit test file
python -m pytest test/unit/structure/utils/test_water_parser.py -v
```

## For Contributors

Architecture contracts and implementation rules:

- `history/architecture/modules/data_contract.md` — output shapes, units, CSV headers
- `history/architecture/modules/glossary_units.md` — terminology and unit definitions
- `history/architecture/modules/scripts/**/interface_exposure.md` — public API per module
- `history/architecture/modules/scripts/**/implementation_guidelines.md` — implementation rules
