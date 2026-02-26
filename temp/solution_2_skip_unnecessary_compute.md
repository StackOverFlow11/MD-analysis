# 方案 2：为合并计算增加按需跳过标志

## 对应问题

`problems.md` → 问题 2：`WaterDensity.py` 和 `WaterOrientation.py` 单独调用时浪费计算

## 现状

`_compute_density_orientation_ensemble()` 无条件地同时计算密度和取向：

```python
# _common.py
def _single_frame_density_and_orientation(...) -> tuple[centers_A, rho, orient, path_length]:
    ...
    cos_theta = _compute_bisector_cos_theta_vec(atoms, selected_oxygen, o_to_h, c_unit)  # 始终计算
    ...
```

当 `WaterDensity.py` 单独调用时，取向计算（含 `find_mic` 批量运算）完全浪费。

## 方案：增加 `compute_orientation` 标志

### 步骤 1：修改 `_single_frame_density_and_orientation`

```python
def _single_frame_density_and_orientation(
    atoms: Atoms,
    *,
    start_interface: StartInterface,
    dz_A: float,
    metal_symbols: Iterable[str] | None,
    compute_orientation: bool = True,      # 新增
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    ...
    # --- Density ---（始终计算）
    rho_g_cm3 = counts * mass_per_water_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)

    # --- Orientation ---（按需计算）
    if compute_orientation and selected_oxygen.size > 0:
        cos_theta = _compute_bisector_cos_theta_vec(atoms, selected_oxygen, o_to_h, c_unit)
        cos_sum = np.histogram(selected_delta_A, bins=edges_A, weights=cos_theta)[0].astype(float)
        orient_g_cm3 = cos_sum * mass_per_water_g / (bin_volumes_A3 * ANGSTROM3_TO_CM3)
    else:
        orient_g_cm3 = np.zeros_like(centers_A, dtype=float)

    return centers_A, rho_g_cm3, orient_g_cm3, path_length_A
```

### 步骤 2：修改 `_compute_density_orientation_ensemble`

```python
def _compute_density_orientation_ensemble(
    xyz_path: Path,
    md_inp_path: Path,
    *,
    start_interface: StartInterface = "low_c",
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
    compute_orientation: bool = True,      # 新增
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    ...
    for atoms in _iter_trajectory(xyz_path, a_A, b_A, c_A):
        per_frame.append(
            _single_frame_density_and_orientation(
                atoms,
                start_interface=start_interface,
                dz_A=dz_A,
                metal_symbols=metal_symbols,
                compute_orientation=compute_orientation,    # 透传
            )
        )
    ...
```

### 步骤 3：`WaterDensity.py` 传入 `compute_orientation=False`

```python
# WaterDensity.py
common_centers_u, mean_path_A, rho_ensemble, _ = _compute_density_orientation_ensemble(
    xyz_path, md_inp_path,
    start_interface=start_interface,
    dz_A=dz_A,
    metal_symbols=metal_symbols,
    compute_orientation=False,    # 跳过取向计算
)
```

`WaterOrientation.py` 保持 `compute_orientation=True`（默认值）。

## 收益

- `water_mass_density_z_distribution_analysis()` 单独调用时跳过所有 `find_mic` + bisector 计算
- `plot_water_three_panel_analysis()` 需要两者，传默认值即可，行为不变
- 向后兼容：默认 `True`，现有调用方无需修改

## 验证

```bash
python -m pytest test/ -v
```
