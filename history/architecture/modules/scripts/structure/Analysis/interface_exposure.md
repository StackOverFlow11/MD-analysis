# `scripts/structure/Analysis/` 接口暴露约定（当前实现）

> 对应代码：`scripts/structure/Analysis/__init__.py` 与
> `scripts/structure/Analysis/WaterAnalysis/__init__.py`
>
> 本文档定义 `scripts.structure.Analysis` 子包的公开接口边界。

## 1. 接口角色定义

- `scripts.structure.Analysis` 是高层分析工作流子包。
- 本层负责多帧统计组织与结果导出（如 CSV）。
- 本层可组合 `scripts.structure.utils` 的底层能力，但不重写其物理定义。
- 子目录文档必须镜像代码目录：`Analysis/WaterAnalysis/` 由独立接口/实现文档维护。

## 2. 当前公开接口清单

### 2.1 配置常量（Stable）

- `DEFAULT_OUTPUT_DIR`
- `DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME`
- `DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME`
- `DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME`
- `DEFAULT_START_INTERFACE`
- `DEFAULT_WATER_MASS_DENSITY_CSV_NAME`
- `DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME`
- `dz_A` 默认继承 `scripts.structure.utils.config.DEFAULT_Z_BIN_WIDTH_A`（当前为 `0.1` Angstrom）

### 2.2 分析函数（Stable）

- `water_mass_density_z_distribution_analysis(...)`
  - 文件位置：`scripts/structure/Analysis/WaterAnalysis/WaterDensity.py`
  - 统计口径：A 口径（逐帧等权系综平均）
- `water_orientation_weighted_density_z_distribution_analysis(...)`
  - 文件位置：`scripts/structure/Analysis/WaterAnalysis/WaterOrientation.py`
  - 统计口径：A 口径（逐帧等权系综平均）
- `detect_adsorbed_layer_range_from_density_profile(...)`
  - 文件位置：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 主峰定义：直接取水密度分布最高 bin 位置
- `ad_water_orientation_analysis(...)`
  - 文件位置：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 功能：自动识别吸附层区间并输出取向分析结果
- `compute_adsorbed_water_theta_distribution(...)`
  - 文件位置：`scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 功能：统计吸附层内 `0-180` 度取向分布（PDF）

## 3. 推荐导入方式

- `from scripts.structure.Analysis import water_mass_density_z_distribution_analysis`
- `from scripts.structure.Analysis import DEFAULT_START_INTERFACE`

## 4. 非公开边界

- 未在对应 `__init__.py` 的 `__all__` 中声明的符号，不属于公开接口。
- 内部 helper（如 `_parse_*`、`_single_frame_*`）仅供模块内复用。

## 5. 接口变更流程（必须执行）

1. 更新 `scripts/structure/Analysis/**/__init__.py` 导出与 `__all__`
2. 更新本文档公开接口清单
3. 若涉及统计口径变化，同步更新：
   - `history/architecture/modules/data_contract.md`
4. 运行导入烟雾测试与相关回归测试
