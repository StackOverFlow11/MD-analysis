# `scripts/structure/Analysis/WaterAnalysis/` 内部实现准则（当前实现口径）

> 适用范围：`scripts/structure/Analysis/WaterAnalysis/`。

## 1. 职责边界

- 本层负责“水分析工作流函数”的实现与组织。
- 本层不得重新定义 `utils` 的底层物理口径（如水分子识别、密度单位定义）。

## 2. 文件职责

- `WaterDensity.py`
  - 实现水密度分析主流程
  - 实现单帧到系综平均的组织逻辑
  - 实现结果 CSV 导出
- `WaterOrientation.py`
  - 实现取向加权密度分析主流程
  - 实现单帧到系综平均的组织逻辑
  - 实现结果 CSV 导出
- `AdWaterOrientation.py`
  - 基于密度分布自动识别吸附层区间
  - 基于吸附层区间导出吸附层取向分析结果
  - 提供吸附层内 `0-180` 度取向分布统计（PDF）
- `__init__.py`
  - 只做稳定导出，不做业务计算

## 3. 统计与数值口径

- `water_mass_density_z_distribution_analysis(...)` 固定采用 A 口径：
- `water_orientation_weighted_density_z_distribution_analysis(...)` 固定采用 A 口径：
  - 每帧先独立计算 profile
  - 界面起点使用直接面向非金属环境的 interface 对（每侧一层）
  - 对 `nbins` 不一致帧先重采样到统一归一化网格
  - 再按帧等权平均
- 吸附层主峰位置判定固定为：水密度分布最高 bin 的距离坐标。
- 吸附层上界为主峰后的第一个局部极小值（平滑后），下界为主峰前最后一个近零点。
- 吸附层取向分布角域固定为 `0-180` 度，输出 PDF 单位为 `degree^-1`。
- 水质量密度输出单位必须保持 `g/cm^3`。
- 取向加权密度输出单位必须保持 `1/Angstrom^3`。

## 4. 依赖方向约束

- 允许：`WaterAnalysis` -> `scripts.structure.utils`
- 禁止：`scripts.structure.utils` -> `WaterAnalysis`
- 禁止引入与分析无关的上层耦合。

## 5. 变更同步要求

- 以下变更必须同步到 `history`：
  - 新增/删除公开函数
  - 统计口径变化（包括平均权重与重采样规则）
  - 输出列与单位变化
