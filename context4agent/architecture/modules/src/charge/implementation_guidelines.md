# md_analysis.charge — Implementation Guidelines

## Layer dependency

- `md_analysis.charge` depends on `md_analysis.utils` (BaderParser, LayerParser)
- `md_analysis.charge` does NOT depend on `md_analysis.water` or `md_analysis.potential`

## Module layout

Flat structure (no sub-packages):
- `config.py` — unit conversion constant + default filenames
- `ChargeAnalysis.py` — atom selectors, single-frame analysis, trajectory analysis

## Data flow

### Single frame (`compute_frame_surface_charge`)

1. Validate `"bader_net_charge" in atoms.arrays`
2. `detect_interface_layers(atoms)` → 2 interface layers
3. Surface area: `|a × b|` (orthogonal cell, normal="c")
4. Per interface layer: `σ = Σ net_charge[layer_indices] / area`
5. Convert e/Å² → μC/cm² via `E_PER_A2_TO_UC_PER_CM2`
6. Store in `atoms.info["surface_charge_density_e_A2"]` and `atoms.info["surface_charge_density_uC_cm2"]`
7. Optional: apply `AtomSelector`, store indices & charges in `atoms.info`

### Trajectory (`trajectory_charge_analysis`)

1. Discover `calc_t*_i*` subdirectories via `sorted(root.glob(dir_pattern))`
2. Per frame: `load_bader_atoms()` → `compute_frame_surface_charge()` → collect
3. Skip frames with missing ACF.dat (warning)
4. Validate `n_selected` consistent across frames
5. Aggregate → mean/std → write CSV → return `TrajectoryChargeResult`

## Unit convention

- Internal computation: e/Å² (`σ = Σ q_i / A`)
- Output: μC/cm² (`1 e/Å² = 1.602176634 × 10³ μC/cm²`)

## Atom selector contract

- `AtomSelector.select(atoms) -> np.ndarray` returns 0-based indices
- May return empty array (no match)
- `ElementSelector`: match by `atoms.get_chemical_symbols()`
- `IndexSelector`: filter to valid range `[0, len(atoms))`

## CSV output format

### surface_charge_density.csv

```
frame,sigma_bottom_uC_cm2,sigma_top_uC_cm2
calc_t000_i000,<value>,<value>
...
mean,<value>,<value>
std,<value>,<value>
```

### selected_atom_charges.csv (optional, when selector provided)

```
frame,atom_0,atom_62,...
calc_t000_i000,<value>,<value>,...
...
mean,<value>,<value>,...
```
