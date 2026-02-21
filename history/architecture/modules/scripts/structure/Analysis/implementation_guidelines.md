# `scripts/structure/Analysis/` 内部实现准则（当前实现口径）

> 适用范围：`scripts/structure/Analysis/`。
>
> 目标：承载“多帧统计 + 结果导出”流程，同时保持与 `utils` 契约一致。

## 1. 职责边界

- `Analysis` 负责：
  - 多帧遍历
  - 系综统计组织
  - CSV/文本等结果导出
- `Analysis` 不负责：
  - 重新定义底层物理量
  - 覆盖 `utils` 的核心算法口径

## 2. 目录组织约定

- `Analysis/config.py`：分析层默认参数（输出目录、默认界面选择、默认文件名）
- `Analysis/Water.py`：水分析的集成绘图入口（组合多个分析结果并输出 PNG/CSV/TXT）
- `Analysis/WaterAnalysis/WaterDensity.py`：水密度相关分析实现
- `Analysis/WaterAnalysis/WaterOrientation.py`：水取向相关分析实现
- `Analysis/WaterAnalysis/AdWaterOrientation.py`：吸附层区间识别与吸附层取向分析实现
- `Analysis/WaterAnalysis/_common.py`：WaterAnalysis 私有公共实现（读轨迹、解析 cell、界面检测、系综平均等）
- `Analysis/WaterAnalysis/__init__.py`：子模块公开导出
- `Analysis/__init__.py`：子包聚合导出
- `Analysis/WaterAnalysis/` 对应文档必须位于：
  - `history/architecture/modules/scripts/structure/Analysis/WaterAnalysis/interface_exposure.md`
  - `history/architecture/modules/scripts/structure/Analysis/WaterAnalysis/implementation_guidelines.md`

## 3. 系综平均口径约定

- `water_mass_density_z_distribution_analysis(...)` 采用 A 口径：
  - 每帧先独立统计 profile
  - 再映射到统一归一化路径坐标
  - 最后按帧等权平均（equal-weight over frames）
- 当每帧 `nbins` 不一致时，必须先重采样到公共网格再平均。
- 若调用方未显式传入 `dz_A`，默认使用 `DEFAULT_Z_BIN_WIDTH_A`（当前 `0.1` Angstrom）。

补充（当前实现特性，需保持一致）：

- 密度与取向加权密度的系综平均在一次轨迹遍历中同时完成（避免重复读帧）。
- `plot_water_three_panel_analysis(...)` 固定两次遍历轨迹：
  - 第一次：密度 + 取向加权密度（共同通道）
  - 第二次：吸附层内 `$\\theta$` 分布（需按距离窗口筛选水分子）

## 4. 依赖方向约束

- 允许：`Analysis` -> `utils`
- 禁止：`utils` -> `Analysis`
- 禁止在 `Analysis` 中反向依赖 `scripts.structure` facade

## 5. 契约同步要求

- 以下变化必须同步更新 `data_contract.md`：
  - 输出列定义
  - 系综平均口径
  - 默认参数行为（影响外部可复现性）
