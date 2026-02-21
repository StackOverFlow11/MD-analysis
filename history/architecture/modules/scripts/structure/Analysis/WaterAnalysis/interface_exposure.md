# `scripts/structure/Analysis/WaterAnalysis/` 接口暴露约定（当前实现）

> 对应代码：
> - `scripts/structure/Analysis/WaterAnalysis/__init__.py`
> - `scripts/structure/Analysis/WaterAnalysis/WaterDensity.py`
> - `scripts/structure/Analysis/WaterAnalysis/WaterOrientation.py`
> - `scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
> - `scripts/structure/Analysis/WaterAnalysis/_common.py`（私有实现，不属于对外契约）

## 1. 接口角色定义

- `WaterAnalysis` 是 `Analysis` 子层中“水相关分析”的聚合目录。
- 本层负责暴露水分析入口函数，不暴露内部 helper。
- 多帧遍历/界面检测/系综平均等公共实现集中在私有模块 `_common.py` 中（外部不得依赖）。

## 2. 当前公开接口清单

### 2.1 分析函数（Stable）

- `water_mass_density_z_distribution_analysis(...)`
  - 定义文件：`WaterDensity.py`
  - 统计口径：A 口径（逐帧等权系综平均，界面取水侧 interface 对）
- `water_orientation_weighted_density_z_distribution_analysis(...)`
  - 定义文件：`WaterOrientation.py`
  - 统计口径：A 口径（逐帧等权系综平均，界面取水侧 interface 对）
- `detect_adsorbed_layer_range_from_density_profile(...)`
  - 定义文件：`AdWaterOrientation.py`
  - 主峰定义：直接取水密度分布最高 bin 位置
- `ad_water_orientation_analysis(...)`
  - 定义文件：`AdWaterOrientation.py`
  - 功能：自动识别吸附层范围并导出吸附层取向分析文件
- `compute_adsorbed_water_theta_distribution(...)`
  - 定义文件：`AdWaterOrientation.py`
  - 功能：统计吸附层内 `0-180` 度取向分布并导出 CSV

### 2.2 非公开接口（Non-public）

- `WaterAnalysis/_common.py` **整个模块**均为非公开（包含轨迹读取、cell 解析、界面检测、单帧计算与系综平均等细节）。
- 任意 `_` 前缀符号（函数/常量/模块）均为非公开。
- 未在 `WaterAnalysis/__init__.py` 的 `__all__` 中声明的符号，均不属于公开接口。

## 3. 推荐导入方式

- `from scripts.structure.Analysis.WaterAnalysis import water_mass_density_z_distribution_analysis`
- `from scripts.structure.Analysis.WaterAnalysis import water_orientation_weighted_density_z_distribution_analysis`
- `from scripts.structure.Analysis.WaterAnalysis import detect_adsorbed_layer_range_from_density_profile`
- `from scripts.structure.Analysis.WaterAnalysis import ad_water_orientation_analysis`
- `from scripts.structure.Analysis import water_mass_density_z_distribution_analysis`
- `from scripts.structure.Analysis import water_orientation_weighted_density_z_distribution_analysis`
- `from scripts.structure.Analysis import detect_adsorbed_layer_range_from_density_profile`
- `from scripts.structure.Analysis import ad_water_orientation_analysis`

## 4. 兼容承诺与变更流程

- 公开接口以 `__all__` 为准。
- 变更流程：
  1. 更新 `WaterAnalysis/__init__.py` 导出与 `__all__`
  2. 更新本文档公开清单
  3. 若统计口径变化，同步 `history/architecture/modules/data_contract.md`
