# `scripts/structure/utils/` 接口暴露约定（当前实现）

> 对应代码：`scripts/structure/utils/__init__.py`
>
> 本文档定义 `scripts.structure.utils` 的符号级公开接口与暴露边界。

## 1. 接口角色定义

- `scripts.structure.utils` 是“底层实现能力的公开导出层”。
- 本层暴露的符号允许被上层（`scripts.structure`）和外部调用方直接导入。
- 本层不自动暴露模块内私有 helper；公开范围由 `__all__` 严格定义。

## 2. 当前公开接口清单（按模块）

### 2.1 `config.py` 常量（Stable）

- `TRANSITION_METAL_SYMBOLS`
- `DEFAULT_METAL_SYMBOLS`
- `DEFAULT_Z_BIN_WIDTH_A`
- `DEFAULT_THETA_BIN_DEG`
- `DEFAULT_WATER_OH_CUTOFF_A`
- `WATER_MOLAR_MASS_G_PER_MOL`

说明：

- 以上常量为默认行为来源；默认值变化属于接口行为变化。
- 当前默认值（与 `scripts/structure/utils/config.py` 对齐）：
  - `DEFAULT_Z_BIN_WIDTH_A = 0.1` Angstrom
  - `DEFAULT_THETA_BIN_DEG = 5.0` degree
  - `DEFAULT_WATER_OH_CUTOFF_A = 1.25` Angstrom

### 2.2 `LayerParser.py` 导出（Stable）

数据结构与异常：

- `Layer`
- `SurfaceDetectionResult`
- `SurfaceGeometryError`

函数：

- `detect_interface_layers(...)`
  - 输入：单帧 `ase.Atoms` + 金属符号/法向/聚类参数
  - 输出：`SurfaceDetectionResult`
  - 语义：仅标记直接面向非金属环境的两层界面层（每侧一层）
- `format_detection_summary(...)`
  - 输入：`SurfaceDetectionResult`
  - 输出：可读文本摘要 `str`

### 2.3 `WaterParser.py` 导出（Stable）

异常：

- `WaterTopologyError`

函数：

- `detect_water_molecule_indices(...)`
  - 输出：`(n_water, 3)`，每行为 `[O_idx, H1_idx, H2_idx]`
- `get_water_oxygen_indices_array(...)`
  - 输出：`(n_water, 1)` 氧索引数组
- `compute_water_mass_density_z_distribution(...)`
  - 输出：`(nbins, 1)`，单位 `g/cm^3`
- `compute_water_orientation_weighted_density_z_distribution(...)`
  - 输出：`(nbins, 1)`，单位 `1/Angstrom^3`
- `compute_water_orientation_theta_pdf_in_c_fraction_window(...)`
  - 输出：`(180 / ndeg,)`，单位 `degree^-1`

## 3. 推荐导入方式

- `from scripts.structure.utils import detect_interface_layers`
- `from scripts.structure.utils import detect_water_molecule_indices`
- `from scripts.structure.utils import compute_water_orientation_theta_pdf_in_c_fraction_window`
- `from scripts.structure.utils import DEFAULT_Z_BIN_WIDTH_A, DEFAULT_THETA_BIN_DEG`

## 4. 非公开边界（必须遵守）

- 未出现在 `scripts/structure/utils/__init__.py` 的 `__all__` 中的符号，不属于公开接口。
- 以下类别默认非公开：
  - `_` 前缀 helper
  - 模块内部中间计算函数
  - 仅用于局部复用的工具函数

调用方不得依赖非公开符号路径（包括“可导入但未导出”的实现细节）。

## 5. 稳定性与兼容策略

### 5.1 稳定性等级

- 当前 `__all__` 内符号按 Stable 管理。

### 5.2 兼容策略

- 默认不破坏现有导入路径：
  - `from scripts.structure.utils import <symbol>`
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

1. 更新 `scripts/structure/utils/__init__.py` 的导出与 `__all__`
2. 更新本文档公开接口清单
3. 若涉及契约变化，同步更新：
   - `history/architecture/modules/data_contract.md`
   - `history/architecture/modules/glossary_units.md`（若涉及术语/单位）
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
