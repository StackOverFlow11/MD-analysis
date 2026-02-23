# `scripts/structure/` 接口暴露约定（当前实现）

> 对应代码：`scripts/structure/__init__.py`
>
> 本文档定义 `scripts.structure` 作为领域入口层的公开接口边界。

## 1. 接口角色定义

- `scripts.structure` 是结构分析能力的统一入口（facade）。
- 本层负责"对外暴露稳定符号"，不负责实现算法本体。
- 算法实现位于 `scripts.structure.utils`；本层做聚合与重导出。
- 面向分析工作流的高层 I/O 能力位于子包 `scripts.structure.Analysis`（不要求全部通过 facade 重导出）。

## 2. 当前公开接口清单

### 2.1 数据结构与异常（Stable）

- `Layer`
- `SurfaceDetectionResult`
- `SurfaceGeometryError`
- `WaterTopologyError`

### 2.2 Layer 相关函数（Stable）

- `detect_interface_layers(...)`
- `format_detection_summary(...)`

### 2.3 Water 相关函数（Stable）

- `detect_water_molecule_indices(...)`
- `get_water_oxygen_indices_array(...)`
- `compute_water_mass_density_z_distribution(...)`
- `compute_water_orientation_weighted_density_z_distribution(...)`
- `compute_water_orientation_theta_pdf_in_c_fraction_window(...)`

### 2.4 默认配置常量（Stable）

- `DEFAULT_METAL_SYMBOLS`
- `DEFAULT_Z_BIN_WIDTH_A`
- `DEFAULT_THETA_BIN_DEG`
- `DEFAULT_WATER_OH_CUTOFF_A`
- `WATER_MOLAR_MASS_G_PER_MOL`

### 2.5 Analysis 子包（Stable，子包路径）

- `scripts.structure.Analysis.plot_water_three_panel_analysis(...)`
  - 当前实现文件：`scripts/structure/Analysis/Water.py`
  - 输出：三联图 PNG + 相关 CSV/TXT（见 `data_contract.md`）
- `scripts.structure.Analysis.water_mass_density_z_distribution_analysis(...)`
  - 当前实现文件：`scripts/structure/Analysis/WaterAnalysis/WaterDensity.py`
  - 统计口径：A 口径（逐帧等权系综平均）
  - 默认参数来源：`scripts/structure/Analysis/config.py`
- `scripts.structure.Analysis.water_orientation_weighted_density_z_distribution_analysis(...)`
  - 当前实现文件：`scripts/structure/Analysis/WaterAnalysis/WaterOrientation.py`
  - 统计口径：A 口径（逐帧等权系综平均）
- `scripts.structure.Analysis.ad_water_orientation_analysis(...)`
  - 当前实现文件：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 输出：吸附层 profile CSV + range TXT（见 `data_contract.md`）
- `scripts.structure.Analysis.compute_adsorbed_water_theta_distribution(...)`
  - 当前实现文件：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 输出：吸附层 $\theta$ 分布 CSV（见 `data_contract.md`）
- `scripts.structure.Analysis.detect_adsorbed_layer_range_from_density_profile(...)`
  - 当前实现文件：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 说明：从密度剖面中推断吸附层距离范围

接口语义说明：

- 具体输入输出 shape、单位与统计口径，以
  `history/architecture/modules/data_contract.md` 为准。

## 3. 推荐导入方式

- `from scripts.structure import detect_interface_layers`
- `from scripts.structure import detect_water_molecule_indices`
- `from scripts.structure import compute_water_orientation_theta_pdf_in_c_fraction_window`
- `from scripts.structure import DEFAULT_THETA_BIN_DEG`
- `from scripts.structure.Analysis import water_mass_density_z_distribution_analysis`
- `from scripts.structure.Analysis import plot_water_three_panel_analysis`

## 4. 非公开接口边界

- 未在 `scripts/structure/__init__.py` 的 `__all__` 中声明的符号，不视为公开接口。
- `_` 前缀 helper、内部中间函数、内部常量，不得作为对外依赖目标。
- 对 Analysis 子包，未在对应 `__init__.py` 的 `__all__` 中声明的符号，同样不视为公开接口。

## 5. 稳定性与兼容承诺

### 5.1 稳定性等级

- 当前 `__all__` 内所有导出按 Stable 管理。

### 5.2 兼容策略

- 默认保持导入路径稳定（`from scripts.structure import ...`）。
- 若需要重命名/拆分接口：
  - 优先新增新接口并保留旧接口兼容期
  - 完成迁移后再移除旧接口

## 6. 接口变更触发条件

以下任一情况视为本层接口变更：

- `__all__` 内容变化
- 导出符号新增/删除/重命名
- 导出函数语义变化（即使签名不变）
- 常量默认值变化（影响外部默认行为）

## 7. 接口变更流程（必须执行）

1. 更新 `scripts/structure/__init__.py` 导出与 `__all__`
2. 更新本文档导出清单与稳定性条款
3. 若涉及契约变化，同步更新：
   - `history/architecture/modules/data_contract.md`
   - `history/architecture/modules/glossary_units.md`（若涉及术语/单位）
4. 运行导入烟雾测试：
   - `import scripts.structure`
   - 关键导出存在性检查
5. 执行回归测试并确认通过

## 8. 反模式（禁止）

- 将 `utils` 私有 helper 直接暴露到 `scripts.structure`
- 在未更新 `__all__` 的情况下改变导出行为
- 修改默认常量而不更新契约文档
- 接口变更后不做导入回归

## 9. 提交前检查清单

- [ ] 公开符号是否全部在 `__all__` 中
- [ ] 非公开符号是否未意外暴露
- [ ] 导出层是否仍保持"薄层聚合"定位
- [ ] 文档与测试是否同步完成
