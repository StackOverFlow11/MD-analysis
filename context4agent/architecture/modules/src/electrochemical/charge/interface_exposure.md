# md_analysis.charge — Interface Exposure

## Public API

| Symbol                              | Module                  | Description                                                     |
|-------------------------------------|-------------------------|-----------------------------------------------------------------|
| `compute_frame_surface_charge`      | Bader/SurfaceCharge.py  | Single-frame surface charge density → `atoms.info` (`method` selects algorithm) |
| `trajectory_surface_charge`         | Bader/SurfaceCharge.py  | Multi-frame surface charge density time series → `(t, 2)` μC/cm² |
| `surface_charge_analysis`           | Bader/SurfaceCharge.py  | End-to-end surface charge analysis → CSV + PNG output           |
| `frame_indexed_atom_charges`        | Bader/AtomCharges.py    | Single-frame net charges for caller-specified atom indices → `(N, 2)` (POSCAR order) |
| `trajectory_indexed_atom_charges`   | Bader/AtomCharges.py    | Per-frame net charges for caller-specified atom indices → `(t, N, 2)` (POSCAR order) |
| `tracked_atom_charge_analysis`      | Bader/AtomCharges.py    | Track XYZ-indexed atoms across trajectory → CSV + PNG (XYZ order) |
| `counterion_charge_analysis`        | Bader/AtomCharges.py    | Per-frame counterion detection + charge tracking → CSV + PNG (XYZ order) |
| `BaderTrajectoryData`               | Bader/BaderData.py      | Frozen dataclass: steps, times, atom_indices_xyz, net_charges (XYZ order) |
| `load_bader_trajectory`             | Bader/BaderData.py      | Load all frames → BaderTrajectoryData with IndexMap remap to XYZ order |
| `E_PER_A2_TO_UC_PER_CM2`           | config.py               | Unit conversion: 1 e/Å² = 1.602176634×10³ μC/cm²               |
| `DEFAULT_SURFACE_CHARGE_CSV_NAME`   | config.py               | Default CSV output filename                                     |
| `DEFAULT_N_SURFACE_LAYERS`          | config.py               | Default number of surface layers per interface (1)              |
| `DEFAULT_SURFACE_CHARGE_PNG_NAME`   | config.py               | Default PNG output filename                                     |

## Config Constants

All default filenames and constants are defined in `config.py`:

| Constant                          | Value                |
|-----------------------------------|----------------------|
| `E_PER_A2_TO_UC_PER_CM2`         | `1.602176634e3`      |
| `DEFAULT_DIR_PATTERN`             | `bader_t*_i*`         |
| `DEFAULT_STRUCTURE_FILENAME`      | `POSCAR`             |
| `DEFAULT_ACF_FILENAME`            | `ACF.dat`            |
| `DEFAULT_POTCAR_FILENAME`         | `POTCAR`             |
| `DEFAULT_SURFACE_CHARGE_CSV_NAME` | `surface_charge.csv` |
| `DEFAULT_N_SURFACE_LAYERS`        | `1`                  |
| `DEFAULT_SURFACE_CHARGE_PNG_NAME` | `surface_charge.png` |

## `compute_frame_surface_charge` Signature

```python
def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",   # "a" | "b" | "c"
    method: str = "counterion",  # "counterion" | "layer"
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
) -> Atoms
```

`method` options (validated against `CHARGE_METHOD_COUNTERION` / `CHARGE_METHOD_LAYER` constants from `utils.config`):
- `"counterion"` — excludes water and metal; only counterion/solute species contribute to σ. `n_surface_layers` is ignored.
- `"layer"` — sums net charges of the N outermost metal layers per interface / area. `n_surface_layers` controls how many layers inward from each interface are included (default 1 = outermost layer only).

Results in `atoms.info`:
- `surface_charge_density_e_A2`: `[σ_aligned, σ_opposed]` (e/Å²)
- `surface_charge_density_uC_cm2`: `[σ_aligned, σ_opposed]` (μC/cm²)
- `n_charged_atoms_per_surface`: `[n_aligned, n_opposed]`
- `charge_per_surface_e`: `[Σq_aligned, Σq_opposed]` (e)

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
    method: str = "counterion",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray   # (t, 2): [:,0]=σ_aligned, [:,1]=σ_opposed, μC/cm²
```

## `surface_charge_analysis` Signature

```python
def surface_charge_analysis(
    root_dir: str | Path = ".",
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = "counterion",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = DEFAULT_N_SURFACE_LAYERS,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    output_dir: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path   # path to written CSV
```

> **Calibration integration**: 函数末尾自动尝试加载 `~/.config/md_analysis/calibration.json`。
> 若存在有效标定，CSV 追加 `phi_aligned_V_vs_SHE`、`phi_opposed_V_vs_SHE` 及其累积平均列；
> PNG 添加右轴显示外推电势。无标定时仅输出 σ 列。

## `BaderTrajectoryData` Dataclass

```python
@dataclass(frozen=True)
class BaderTrajectoryData:
    steps: np.ndarray            # (n_frames,) int — MD step
    times: np.ndarray            # (n_frames,) int — time in fs
    atom_indices_xyz: np.ndarray # (n_atoms,) int — [0, 1, ..., N-1] XYZ order
    net_charges: np.ndarray      # (n_frames, n_atoms) float — XYZ order net charges
```

## `load_bader_trajectory` Signature

```python
def load_bader_trajectory(
    root_dir: str | Path = ".",
    *,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> BaderTrajectoryData
```

## `tracked_atom_charge_analysis` Signature

```python
def tracked_atom_charge_analysis(
    root_dir: str | Path = ".",
    *,
    atom_indices_xyz: Iterable[int],   # XYZ 0-based atom indices
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    output_dir: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path   # path to written CSV
```

## `counterion_charge_analysis` Signature

```python
def counterion_charge_analysis(
    root_dir: str | Path = ".",
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    output_dir: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path   # path to per-frame CSV
```
