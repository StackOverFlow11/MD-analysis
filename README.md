# MD Analysis

Lightweight analysis utilities for periodic metal-water interfaces, focused on:

- water mass-density profile along interface distance
- water orientation-weighted profile along interface distance
- adsorbed-layer detection and orientation-angle distribution
- integrated three-panel visualization output

## Who This Is For

- New contributors who want a minimal, working entrypoint
- Users who just cloned the repo and want to run one end-to-end example

## Quick Start (5 minutes)

### 1) Prepare environment

Install required dependencies in your Python environment (at least `numpy`, `matplotlib`, `ase`, `pytest`).

### 2) Run one end-to-end plot generation

From repo root:

```bash
python test/integration/structure/Analysis/test_water_three_panel_plot.py
```

Generated figure:

- `test/_tmp_preview/water_three_panel_analysis.png`

### 3) Example input data

Current integration examples use:

- `data_example/potential/md-pos-1.xyz`
- `data_example/potential/md.inp`

## Recommended Entry Function

Primary integrated plotting entry:

- `scripts.structure.Analysis.plot_water_three_panel_analysis(...)`

This function integrates:

1. water density profile
2. orientation-weighted profile
3. adsorbed-layer theta distribution

## Project Layout (Practical View)

- `scripts/structure/utils/`
  - low-level geometry and water parsing/statistics
- `scripts/structure/Analysis/`
  - workflow-level analysis and plotting composition
- `test/integration/structure/Analysis/`
  - runnable integration scripts and plot-generation tests
- `history/`
  - architecture, contracts, and implementation rules

## Notes For Contributors

- Public interface and behavior contracts are documented under:
  - `history/architecture/modules/scripts/**/interface_exposure.md`
  - `history/architecture/modules/scripts/**/implementation_guidelines.md`
- Global cross-module contracts are centralized in:
  - `history/architecture/modules/data_contract.md`
  - `history/architecture/modules/glossary_units.md`
