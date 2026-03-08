# `md_analysis.utils` 接口暴露约定（当前实现）

> 对应代码：`src/md_analysis/utils/__init__.py`
>
> 本文档定义 `md_analysis.utils` 的符号级公开接口与暴露边界。

## 1. 接口角色定义

- `md_analysis.utils` 是"底层实现能力的公开导出层"。
- 本层暴露的符号允许被上层（`md_analysis.water`、`md_analysis.potential`）和外部调用方直接导入。
- 本层不自动暴露模块内私有 helper；公开范围由 `__all__` 严格定义。

## 2. 当前公开接口清单（按模块）

### 2.1 `config.py` 常量（Stable）

- `TRANSITION_METAL_SYMBOLS`
- `DEFAULT_METAL_SYMBOLS`
- `DEFAULT_Z_BIN_WIDTH_A`
- `DEFAULT_THETA_BIN_DEG`
- `DEFAULT_WATER_OH_CUTOFF_A`
- `WATER_MOLAR_MASS_G_PER_MOL`
- `HA_TO_EV`
- `BOHR_TO_ANG`
- `DP_A_H3O_W_EV`
- `MU_HPLUS_G0_EV`
- `DELTA_E_ZP_EV`

说明：

- 以上常量为默认行为来源；默认值变化属于接口行为变化。
- 当前默认值（与 `src/md_analysis/utils/config.py` 对齐）：
  - `DEFAULT_Z_BIN_WIDTH_A = 0.1` Angstrom
  - `DEFAULT_THETA_BIN_DEG = 5.0` degree
  - `DEFAULT_WATER_OH_CUTOFF_A = 1.25` Angstrom

### 2.2 `ClusterUtils.py` 导出（Stable）

函数：

- `cluster_1d_periodic(...)`
  - 1D 周期性聚类
- `find_largest_gap_periodic(...)`
  - 在周期排列中找到最大间隙
- `gap_midpoint_periodic(...)`
  - 计算间隙中点

### 2.3 `CubeParser.py` 导出（Stable）

数据结构：

- `CubeHeader`

函数：

- `read_cube_header_and_values(...)`
  - 读取 Gaussian cube 文件
- `slab_average_potential_ev(...)`
  - slab 区间平均电势，返回 `(phi_center_ev, info_dict)`
- `plane_avg_phi_z_ev(...)`
  - xy 平面平均的 φ(z) profile
- `z_coords_ang(...)`
  - z 轴网格坐标（Angstrom）
- `extract_step_from_cube_filename(...)`
  - 从 cube 文件名提取 step 编号

### 2.4 `LayerParser.py` 导出（Stable）

数据结构与异常：

- `Layer`
- `SurfaceDetectionResult`
- `SurfaceGeometryError`

函数：

- `detect_interface_layers(...)`
  - 输入：单帧 `ase.Atoms` + 金属符号/法向/聚类参数
  - 输出：`SurfaceDetectionResult`
  - 语义：仅标记直接面向非金属环境的两层界面层（每侧固定 1 层，不可配置）
- `format_detection_summary(...)`
  - 输入：`SurfaceDetectionResult`
  - 输出：可读文本摘要 `str`

### 2.5 `WaterParser.py` 导出（Stable）

异常：

- `WaterTopologyError`

函数：

- `detect_water_molecule_indices(...)`
  - 输出：`(n_water, 3)`，每行为 `[O_idx, H1_idx, H2_idx]`
- `get_water_oxygen_indices_array(...)`
  - 输出：`(n_water, 1)` 氧索引数组

### 2.6 `BaderParser.py` 导出（Stable）

异常：

- `BaderParseError`
  - ACF 文件格式错误或原子数不匹配时抛出

函数：

- `load_bader_atoms(structure_path, acf_path, potcar_path) -> Atoms`
  - 输入：POSCAR/CONTCAR 路径、ACF.dat 路径、POTCAR 路径
  - 输出：增强的 `ase.Atoms`，附加 `atoms.arrays["bader_charge"]`（原始电子数）和 `atoms.arrays["bader_net_charge"]`（ZVAL - bader_charge，正值 = 失去电子）
  - 语义：解析 VASP Bader 电荷分析结果，与结构信息合并

### 2.7 `CellParser.py` 导出（Stable）

异常：

- `CellParseError`
  - CP2K 文件格式错误（缺少 &CELL 块、缺少 A/B/C 向量、非正交 cell、缺少 ABC 行）时抛出

函数：

- `parse_abc_from_restart(restart_path) -> (float, float, float)`
  - 输入：CP2K `.restart` 文件路径
  - 输出：正交 cell 的 `(a, b, c)` 长度（Angstrom）
  - 语义：解析 `&CELL ... &END CELL` 块中 A/B/C 向量的对角元素

- `parse_abc_from_md_inp(md_inp_path) -> (float, float, float)`
  - 输入：CP2K 输入文件路径（如 `md.inp`）
  - 输出：正交 cell 的 `(a, b, c)` 长度（Angstrom）
  - 语义：匹配 `ABC [angstrom] a b c` 行
  - 从 `water/WaterAnalysis/_common.py` 迁移而来，升级为公开 API

> **注**：`_compute_water_mass_density_z_distribution`、`_compute_water_orientation_weighted_density_z_distribution`、
> `_compute_water_orientation_theta_pdf_in_c_fraction_window` 三个函数已降级为内部（`_` 前缀），不再属于公开 API。
> 它们针对全 cell z 轴分箱，与 `water` 层的界面-到-中点分析语义不同，不适合作为公开接口暴露。

### 2.8 `RestartParser/SlowgrowthParser.py` 导出（Stable）

异常：

- `SlowGrowthParseError`
  - 解析 slow-growth restart 或 LagrangeMultLog 文件失败时抛出

数据结构（frozen dataclass）：

- `ColvarDef`
  - 集合变量定义：`cv_type`（`"DISTANCE"` / `"ANGLE"` / `"COMBINE_COLVAR"`）、`atoms`、`function`、`variables`、`components`
- `ConstraintInfo`
  - COLLECTIVE 约束参数：`colvar_id`、`target_au`、`target_growth_au`、`intermolecular`
- `SlowGrowthRestart`
  - restart 文件元数据：`project_name`、`step_start`、`time_start_fs`、`timestep_fs`、`total_steps`、`constraint`、`lagrange_filename`、`colvar`、`cell_abc_ang`、`fixed_atom_indices`
- `LagrangeMultLog`
  - 拉格朗日乘子时序：`shake`、`rattle`、`n_steps`、`n_constraints`；属性 `collective_shake`/`collective_rattle` 提取 CV 乘子

函数：

- `parse_slowgrowth_restart(restart_path) -> SlowGrowthRestart`
  - 解析 CP2K slow-growth restart 文件，内部复用 `CellParser.parse_abc_from_restart()`
- `parse_lagrange_mult_log(log_path) -> LagrangeMultLog`
  - 解析 LagrangeMultLog 文件，自动检测单约束/多约束格式
- `compute_target_series(restart, n_steps) -> np.ndarray`
  - 重建 ξ(t) 序列（原子单位），`xi(k) = target_au + (k - step_start) * target_growth_au`

## 3. 推荐导入方式

- `from md_analysis.utils import detect_interface_layers`
- `from md_analysis.utils import detect_water_molecule_indices`
- `from md_analysis.utils import get_water_oxygen_indices_array`
- `from md_analysis.utils import DEFAULT_Z_BIN_WIDTH_A, DEFAULT_THETA_BIN_DEG`
- `from md_analysis.utils import read_cube_header_and_values, slab_average_potential_ev`
- `from md_analysis.utils import HA_TO_EV, BOHR_TO_ANG`

## 4. 非公开边界（必须遵守）

- 未出现在 `src/md_analysis/utils/__init__.py` 的 `__all__` 中的符号，不属于公开接口。
- 以下类别默认非公开：
  - `_` 前缀 helper
  - 模块内部中间计算函数
  - 仅用于局部复用的工具函数

调用方不得依赖非公开符号路径（包括"可导入但未导出"的实现细节）。

## 5. 稳定性与兼容策略

### 5.1 稳定性等级

- 当前 `__all__` 内符号按 Stable 管理。

### 5.2 兼容策略

- 默认不破坏现有导入路径：
  - `from md_analysis.utils import <symbol>`
- 若需要替换接口：
  - 先新增新符号并保留旧符号兼容期
  - 再进行分阶段迁移

## 6. 接口变更触发条件

以下情况属于 `utils` 公开接口变更：

- `__all__` 新增/删除/重命名
- 公开函数签名变化
- 公开函数输出 shape 或单位变化
- 默认常量值变化（影响默认行为）
- 异常类型或抛错条件变化

## 7. 接口变更流程（必须执行）

1. 更新 `src/md_analysis/utils/__init__.py` 的导出与 `__all__`
2. 更新本文档公开接口清单
3. 若涉及契约变化，同步更新：
   - `context4agent/architecture/modules/data_contract.md`
   - `context4agent/architecture/modules/glossary_units.md`（若涉及术语/单位）
4. 运行导入烟雾测试与回归测试

## 8. 反模式（禁止）

- 将模块私有 helper 暴露到 `__all__`
- 公开函数改签名但不更新文档
- 更改默认常量而不声明行为变化
- 直接在外部依赖内部 helper 路径

## 9. 提交前检查清单

- [ ] `__all__` 与公开清单是否一致
- [ ] 输出 shape/单位变更是否已同步文档
- [ ] 默认常量变更是否已声明
- [ ] 导入测试与回归测试是否通过
