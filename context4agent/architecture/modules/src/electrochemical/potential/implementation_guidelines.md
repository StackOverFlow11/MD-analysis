# md_analysis.potential — Implementation Guidelines

## Layer dependency

- `md_analysis.potential` depends on `md_analysis.utils` (CubeParser, ClusterUtils, config, `_io_helpers`)
- `md_analysis.potential` does NOT depend on `md_analysis.water`

## Module layout

Flat structure (no sub-packages):
- `config.py` — default output filenames + distributed mode constants
- `_frame_source.py` — frame data abstraction: `PotentialFrame` dataclass + `discover_continuous_frames()` / `discover_distributed_frames()`
- `CenterPotential.py` — center slab potential + Fermi + electrode potential + thickness sensitivity (delegates to frame source)
- `PhiZProfile.py` — φ(z) plane-averaged profile analysis (delegates to frame source for distributed mode)

## Input mode abstraction

Two input modes, unified through `PotentialFrame`:

- **Continuous (mode A)**: single directory with `*-HARTREE-*.cube` files + `md.out` → `discover_continuous_frames()`
- **Distributed (mode B)**: `potential_t{time}_i{step}/` subdirectories, each containing one cube file + `sp.out` → `discover_distributed_frames()`

All 5 analysis functions accept `input_mode` parameter (default `"continuous"`) and route to the appropriate discovery function via `_resolve_frames()`. Downstream analysis logic operates on `list[PotentialFrame]` uniformly.

In distributed mode, atoms for interface detection come from the cube file itself (`read_cube_atoms()` in CubeParser), eliminating the need for a separate xyz trajectory.

## Key imports from `utils`

- `discover_cube_files` — 从 `CubeParser` 导入，取代 `CenterPotential.py` 中原有的内联 glob+检查逻辑
- `_float` — 从 `CubeParser` 导入的 Fortran 浮点解析 helper（私有但跨模块使用），用于解析 md.out 中的数值
- `_cumulative_average`, `_write_csv` — 从 `utils._io_helpers` 导入的私有共享 helper，取代原有的模块内重复实现

## Cube file conventions

- CP2K V_HARTREE_CUBE: z is the fastest-running index → reshape as `(nx, ny, nz)`
- Units in cube files: Bohr (positions) and Hartree (values)
- Conversion constants are in `md_analysis.utils.config`

## Interface detection for slab centering

Delegates to `LayerParser.detect_interface_layers` for metal layer clustering
and interface labeling. `_extract_interface_geometry` extracts interface layer
Cartesian coordinates and water gap midpoint from the `SurfaceDetectionResult`.
`gap_midpoint_periodic` is still imported directly from `ClusterUtils`.

## Electrode potential formula (cSHE)

$$
U_{\mathrm{SHE}} = -E_{\mathrm{Fermi}} + \varphi_{\mathrm{center}} + \Delta\Psi_{a}(\mathrm{H_3O^+/w}) - \mu(\mathrm{H^+, g^0}) - \Delta E_{\mathrm{ZP}}
$$

Constants from `md_analysis.utils.config`:
- `DP_A_H3O_W_EV = 15.35`
- `MU_HPLUS_G0_EV = 15.81`
- `DELTA_E_ZP_EV = 0.35`

## slab_average_potential_ev (CubeParser)

Core single-frame computation used by `center_slab_potential_analysis` and `thickness_sensitivity_analysis`:

1. Reshape flat cube values → `(nx, ny, nz)`, xy-plane average → `phi_z_ha` array `(nz,)`
2. Convert `z_center_ang` to bohr relative coordinate `z_center_rel`
3. Compute periodic distance from each z grid point to `z_center_rel`
4. `mask = dist <= half_thickness` selects slab region
5. `phi_center_ev = mean(phi_z_ha[mask]) * HA_TO_EV`
6. `phi_z_std_ev = std(phi_z_ha[mask]) * HA_TO_EV` — spatial std within slab

Returns `(phi_center_ev, info_dict)` where `info_dict` includes `phi_z_std_ev`.

## thickness_sensitivity_analysis workflow

Sweeps slab-averaging thickness to evaluate sensitivity of U vs SHE:

1. **Read & cache**: For each cube frame, read `(header, values)`, detect interface center, look up Fermi energy from `md.out` → cache `(header, values, z_center_ang, fermi_ev)`
2. **Sweep**: For each thickness in `[3.5, thickness_end]` (step 0.5 Å):
   - Per frame: `slab_average_potential_ev()` → `phi_ev`, then `U = -fermi_ev + phi_ev + cSHE_offset`
   - Collect all frames → `mean(U)` (left axis) and `mean(phi_z_std_ev)` (right axis)
3. **Output**: CSV (`thickness_ang, mean_U_vs_SHE_V, mean_phi_z_spatial_std_eV, n_frames`) + dual-axis PNG

Right axis = ensemble-averaged spatial std of φ(z) within the slab window (measures Hartree potential flatness in the averaging region).

## PhiZProfile slab centering

`phi_z_planeavg_analysis` 使用**单次移位**策略对 slab 进行居中：

- `shift_n: int | None = None`，在帧循环前初始化
- 仅在**第一帧**通过 `_slab_center_roll()` 检测 slab 中心并计算 roll 偏移量
- 后续所有帧复用同一 `shift_n` 值
- 目的：避免 MD 原子热涨落导致逐帧偏移量不同，从而造成 min/max envelope 在 overlay 图中错位
- 若未检测到金属元素，`shift_n = 0`（不做 roll）
