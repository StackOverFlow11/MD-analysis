# md_analysis.charge — Implementation Guidelines

## Layer dependency

- `md_analysis.charge` depends on `md_analysis.utils` (BaderParser, LayerParser, WaterParser, config constants, `_io_helpers`)
- `md_analysis.charge` does NOT depend on `md_analysis.water` or `md_analysis.potential`

## Module layout

Flat structure (no sub-packages):
- `config.py` — unit conversion constant + default filenames + output file name constants
- `BaderAnalysis.py` — single-frame surface charge (two methods), single-frame indexed atom charges, trajectory indexed atom charges, trajectory surface charge, end-to-end surface charge analysis

## Key imports from `utils`

- `AXIS_MAP`, `AREA_VECTOR_INDICES` — 轴索引和面积计算向量索引，取代原有的模块局部映射字典
- `CHARGE_METHOD_COUNTERION`, `CHARGE_METHOD_LAYER` — 电荷方法名称常量，取代硬编码字符串
- `_cumulative_average`, `_write_csv` — 从 `utils._io_helpers` 导入的私有共享 helper，取代原有的模块内重复实现

## Surface charge methods

`compute_frame_surface_charge(atoms, *, method="counterion", layer_tol_A=DEFAULT_LAYER_TOL_A)` dispatches to one of (validated against `_VALID_METHODS = (CHARGE_METHOD_COUNTERION, CHARGE_METHOD_LAYER)`):

The `layer_tol_A` parameter is forwarded to `detect_interface_layers()` to control the clustering tolerance (in Å) for grouping metal atoms into layers. This enables consistent layer detection across the full call chain (CLI → analysis → layer detection).

### `method="counterion"` (`_compute_surface_charge_counterion`)

Only non-water, non-metal species (counterions, adsorbates) with non-zero net charge contribute:

1. `detect_interface_layers(atoms, layer_tol_A=layer_tol_A)` → `interface_normal_aligned()` / `interface_normal_opposed()`
2. `detect_water_molecule_indices(atoms)` → water atom set
3. Exclude = water ∪ metal → remaining charged atoms only
4. MIC-based directional assignment to nearest surface within half-gap
5. σ = Σq_assigned / area

### `method="layer"` (`_compute_surface_charge_layer`)

Sum net charges of surface-layer metal atoms (configurable depth via `n_surface_layers`):

1. `detect_interface_layers(atoms, layer_tol_A=layer_tol_A)` → `metal_layers_sorted` ordered `[aligned, interior…, opposed]`
2. Validate `1 <= n_surface_layers <= len(metal_layers_sorted)`
3. aligned side = `layers[:n_surface_layers]`，opposed side = `layers[-n_surface_layers:]`
4. For each side: σ = Σ(net_charge[all_layer_atoms]) / area
5. When `n_surface_layers=1` (default), behavior matches the single-interface-layer case

### Output column ordering

Both methods output `[σ_aligned, σ_opposed]` — ordered by stable `interface_label`, not by `center_frac`. This ensures PBC-shift invariance: a periodic translation does not swap the two columns.

## Data flow

### Single frame (`compute_frame_surface_charge`)

1. Validate `method` ∈ `{"counterion", "layer"}` → `ValueError`
2. Dispatch to `_compute_surface_charge_counterion` or `_compute_surface_charge_layer`
3. Both store results in `atoms.info`: `surface_charge_density_e_A2`, `surface_charge_density_uC_cm2`, `n_charged_atoms_per_surface`, `charge_per_surface_e`

### Single frame indexed (`frame_indexed_atom_charges`)

1. Validate `atom_indices`: must be 1-D, integer, non-negative → `ValueError`
2. Validate `"bader_net_charge" in atoms.arrays` → `ValueError`
3. Validate indices < n_atoms → `IndexError`
4. Return `(N, 2)` ndarray: `[:,0]=index`, `[:,1]=net_charge`

### Trajectory (`trajectory_indexed_atom_charges`)

1. Validate `atom_index_matrix`: must be 2-D, integer, non-negative → `ValueError`
2. Discover `bader_t*_i*` subdirectories via `_sorted_frame_dirs(root, dir_pattern)` (numeric sort by t value)
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
3. Discover `bader_t*_i*` subdirectories via `_sorted_frame_dirs` (numeric sort) → `FileNotFoundError` if none
4. Per frame:
   - Check POSCAR/ACF.dat/POTCAR exist → `FileNotFoundError` (includes frame name)
   - `load_bader_atoms()` → `compute_frame_surface_charge(atoms, ..., method=method)`
   - Collect `atoms.info["surface_charge_density_uC_cm2"]`
5. Stack → return `(t, 2)` ndarray (μC/cm²)

### End-to-end analysis (`surface_charge_analysis`)

1. Validate `normal`, `root_dir`
2. `_sorted_frame_dirs(root, dir_pattern)` → numeric sort by t value
3. Apply `frame_dirs[frame_start:frame_end:frame_step]` slice
4. Per frame: `load_bader_atoms` → `compute_frame_surface_charge(..., method=method)` → collect `(step, σ_aligned, σ_opposed)`
5. Compute cumulative average for both surfaces
6. Write CSV (`_write_csv`) + PNG (`_plot_surface_charge`)
7. Return CSV path

## CLI

`md-analysis charge --charge-method counterion|layer`

Output directory: `<outdir>/charge/<method>/`

## Directory sorting

All trajectory functions use `_sorted_frame_dirs()` which sorts by the numeric value of `_t(\d+)` in the directory name, avoiding lexicographic mis-ordering of non-zero-padded names (e.g., `bader_t1000` before `bader_t50`).

## Unit convention

- Internal computation: e/Å² (`σ = Σ q / A`)
- Output: μC/cm² (`1 e/Å² = 1.602176634 × 10³ μC/cm²`)
