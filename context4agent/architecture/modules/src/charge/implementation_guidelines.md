# md_analysis.charge — Implementation Guidelines

## Layer dependency

- `md_analysis.charge` depends on `md_analysis.utils` (BaderParser, LayerParser)
- `md_analysis.charge` does NOT depend on `md_analysis.water` or `md_analysis.potential`

## Module layout

Flat structure (no sub-packages):
- `config.py` — unit conversion constant + default filenames + output file name constants
- `ChargeAnalysis.py` — single-frame surface charge, single-frame indexed atom charges, trajectory indexed atom charges, trajectory surface charge, end-to-end surface charge analysis

## Data flow

### Single frame (`compute_frame_surface_charge`)

Counterion-based surface charge computation:

1. Validate `normal` ∈ `{"a", "b", "c"}` → `ValueError`
2. Validate `"bader_net_charge" in atoms.arrays`
3. `detect_interface_layers(atoms)` → 2 interface layers (needed for surface positions)
4. Surface area: `|cell[i] × cell[j]|` where `(i, j) = _AREA_VECTORS[normal]`
5. Identify counterion indices: `symbol ∉ metal_set ∪ SOLVENT_SYMBOLS`
6. If no counterions → σ = [0, 0], return early
7. Compute fractional coords along normal; surface positions via circular mean of interface layer atoms
8. Slab midplane via circular mean of all metal fractional coords → gap midpoint = (midplane + 0.5) % 1
9. For each counterion, MIC delta to each surface; check (a) water-side direction, (b) within `cutoff_A`
10. Sum `net_charge[assigned_counterion_indices]` for each surface → σ = sum / area
11. Convert e/Å² → μC/cm² via `E_PER_A2_TO_UC_PER_CM2`
12. Store in `atoms.info`: `surface_charge_density_e_A2`, `surface_charge_density_uC_cm2`, `n_counterions_per_surface`, `counterion_charge_per_surface_e`

### Single frame indexed (`frame_indexed_atom_charges`)

1. Validate `atom_indices`: must be 1-D, integer, non-negative → `ValueError`
2. Validate `"bader_net_charge" in atoms.arrays` → `ValueError`
3. Validate indices < n_atoms → `IndexError`
4. Return `(N, 2)` ndarray: `[:,0]=index`, `[:,1]=net_charge`

### Trajectory (`trajectory_indexed_atom_charges`)

1. Validate `atom_index_matrix`: must be 2-D, integer, non-negative → `ValueError`
2. Discover `calc_t*_i*` subdirectories via `_sorted_frame_dirs(root, dir_pattern)` (numeric sort by t value)
3. Validate `t == len(frame_dirs)` → `ValueError`
4. Per frame:
   - Check POSCAR/ACF.dat/POTCAR exist → `FileNotFoundError` (includes frame name)
   - `load_bader_atoms()` → call `frame_indexed_atom_charges(atoms, arr[i])`
   - Re-raise `IndexError` with frame name context
   - Fill `result[i] = frame_result`
5. Return `(t, N, 2)` ndarray

### Trajectory surface charge (`trajectory_surface_charge`)

1. Validate `normal` ∈ `{"a", "b", "c"}` → `ValueError`
2. Validate `root_dir` exists → `FileNotFoundError`
3. Discover `calc_t*_i*` subdirectories via `_sorted_frame_dirs` (numeric sort) → `FileNotFoundError` if none
4. Per frame:
   - Check POSCAR/ACF.dat/POTCAR exist → `FileNotFoundError` (includes frame name)
   - `load_bader_atoms()` → `compute_frame_surface_charge(atoms, ...)`
   - Collect `atoms.info["surface_charge_density_uC_cm2"]`
5. Stack → return `(t, 2)` ndarray (μC/cm²)

### End-to-end analysis (`surface_charge_analysis`)

1. Validate `normal`, `root_dir`
2. `_sorted_frame_dirs(root, dir_pattern)` → numeric sort by t value
3. Apply `frame_dirs[frame_start:frame_end:frame_step]` slice
4. Per frame: `load_bader_atoms` → `compute_frame_surface_charge` → collect `(step, σ_bottom, σ_top)`
5. Compute cumulative average for both surfaces
6. Write CSV (`_write_csv`) + PNG (`_plot_surface_charge`)
7. Return CSV path

## Directory sorting

All trajectory functions use `_sorted_frame_dirs()` which sorts by the numeric value of `_t(\d+)` in the directory name, avoiding lexicographic mis-ordering of non-zero-padded names (e.g., `calc_t1000` before `calc_t50`).

## Counterion detection

- Counterion = any atom whose symbol ∉ `metal_set ∪ SOLVENT_SYMBOLS`
- `SOLVENT_SYMBOLS = frozenset({"O", "H"})`
- Supports alkali metals (Group 1), halides (Group 17), any species
- No counterions → σ = 0 (neutral electrode, naturally handled by empty sum)

## Directional cutoff (PBC-safe)

Uses minimum image convention (MIC) on fractional coordinates:
- `gap_dir = mic(gap_frac - frac_surface)` — direction from surface toward water region
- `delta = mic(frac_ion - frac_surface)` — direction from surface toward ion
- Assigned if `delta * gap_dir > 0` (water-side) AND `|delta| * cell_len < cutoff_A`
- Default `cutoff_A = 7.0` Å

## Unit convention

- Internal computation: e/Å² (`σ = Σ q_counterion / A`)
- Output: μC/cm² (`1 e/Å² = 1.602176634 × 10³ μC/cm²`)
