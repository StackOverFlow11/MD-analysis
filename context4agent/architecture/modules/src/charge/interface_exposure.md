# md_analysis.charge — Interface Exposure

## Public API

| Symbol                              | Module              | Description                                                     |
|-------------------------------------|---------------------|-----------------------------------------------------------------|
| `compute_frame_surface_charge`      | ChargeAnalysis.py   | Single-frame surface charge density → `atoms.info`              |
| `frame_indexed_atom_charges`        | ChargeAnalysis.py   | Single-frame net charges for caller-specified atom indices → `(N, 2)` |
| `trajectory_indexed_atom_charges`   | ChargeAnalysis.py   | Per-frame net charges for caller-specified atom indices → `(t, N, 2)` |
| `trajectory_surface_charge`         | ChargeAnalysis.py   | Multi-frame surface charge density time series → `(t, 2)` μC/cm² |
| `E_PER_A2_TO_UC_PER_CM2`           | config.py           | Unit conversion: 1 e/Å² = 1.602176634×10³ μC/cm²               |

## Config Constants

All default filenames and constants are defined in `config.py`:

| Constant                     | Value          |
|------------------------------|----------------|
| `E_PER_A2_TO_UC_PER_CM2`    | `1.602176634e3`|
| `DEFAULT_DIR_PATTERN`        | `calc_t*_i*`   |
| `DEFAULT_STRUCTURE_FILENAME` | `POSCAR`       |
| `DEFAULT_ACF_FILENAME`       | `ACF.dat`      |
| `DEFAULT_POTCAR_FILENAME`    | `POTCAR`       |

## `compute_frame_surface_charge` Signature

```python
def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",   # "a" | "b" | "c"
) -> Atoms
```

Results in `atoms.info`:
- `surface_charge_density_e_A2`: `[σ_bottom, σ_top]` (e/Å²)
- `surface_charge_density_uC_cm2`: `[σ_bottom, σ_top]` (μC/cm²)

## `frame_indexed_atom_charges` Signature

```python
def frame_indexed_atom_charges(
    atoms: Atoms,
    atom_indices: np.ndarray,   # (N,), 0-based int
) -> np.ndarray   # (N, 2): [:,0]=index, [:,1]=net_charge
```

## `trajectory_indexed_atom_charges` Signature

```python
def trajectory_indexed_atom_charges(
    root_dir: str | Path,
    atom_index_matrix: np.ndarray,   # (t, N), 0-based int
    *,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray   # (t, N, 2): [:,:,0]=index, [:,:,1]=net_charge
```

## `trajectory_surface_charge` Signature

```python
def trajectory_surface_charge(
    root_dir: str | Path,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray   # (t, 2): [:,0]=σ_bottom, [:,1]=σ_top, μC/cm²
```
