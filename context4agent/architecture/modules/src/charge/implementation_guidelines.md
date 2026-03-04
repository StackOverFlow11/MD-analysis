# md_analysis.charge — Implementation Guidelines

## Layer dependency

- `md_analysis.charge` depends on `md_analysis.utils` (BaderParser, LayerParser)
- `md_analysis.charge` does NOT depend on `md_analysis.water` or `md_analysis.potential`

## Module layout

Flat structure (no sub-packages):
- `config.py` — unit conversion constant + default filenames
- `ChargeAnalysis.py` — single-frame surface charge, trajectory indexed atom charges

## Data flow

### Single frame (`compute_frame_surface_charge`)

1. Validate `normal` ∈ `{"a", "b", "c"}` → `ValueError`
2. Validate `"bader_net_charge" in atoms.arrays`
3. `detect_interface_layers(atoms)` → 2 interface layers
4. Surface area: `|cell[i] × cell[j]|` where `(i, j) = _AREA_VECTORS[normal]`
   - `_AREA_VECTORS = {"a": (1, 2), "b": (0, 2), "c": (0, 1)}`
5. Per interface layer: `σ = Σ net_charge[layer_indices] / area`
6. Convert e/Å² → μC/cm² via `E_PER_A2_TO_UC_PER_CM2`
7. Store in `atoms.info["surface_charge_density_e_A2"]` and `atoms.info["surface_charge_density_uC_cm2"]`

### Trajectory (`trajectory_indexed_atom_charges`)

1. Validate `atom_index_matrix`: must be 2-D, integer, non-negative → `ValueError`
2. Discover `calc_t*_i*` subdirectories via `sorted(root.glob(dir_pattern))`
3. Validate `t == len(frame_dirs)` → `ValueError`
4. Per frame:
   - Check POSCAR/ACF.dat/POTCAR exist → `FileNotFoundError` (includes frame name)
   - `load_bader_atoms()` → read `bader_net_charge`
   - Validate indices < n_atoms → `IndexError` (includes frame name + index value)
   - Fill `result[i, :, 0] = indices`, `result[i, :, 1] = net_charge[indices]`
5. Return `(t, N, 2)` ndarray

## Unit convention

- Internal computation: e/Å² (`σ = Σ q_i / A`)
- Output: μC/cm² (`1 e/Å² = 1.602176634 × 10³ μC/cm²`)
