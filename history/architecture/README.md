# 代码架构（当前仓库结构与数据流）

用途：记录"仓库现在是什么样、模块怎么分、数据怎么流动"，用于协作对齐与避免口径漂移。

## 适用体系前提（当前实现依赖的假设）

- 三基矢**正交**的周期性体系（`md.inp` 中通过 `ABC [angstrom] a b c` 提供正交盒长）
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中存在**两个界面**（两侧表面）

## 仓库根目录结构（当前）

- `scripts/`：可导入的 Python 包（当前主实现都在这里）
- `test/`：pytest 单元/集成测试与可运行示例脚本
- `data_example/`：最小可复现实例数据（当前主要用 `data_example/potential/`）
- `history/`：架构/契约/决策/需求上下文（本目录）
- `README.md`：快速入口与最小跑通方式

## `scripts/` 包内分层（当前）

### 1) `scripts.structure.utils`（单帧/底层）

输入：单帧 `ase.Atoms`

- `LayerParser.py`
  - 金属原子筛选 → 沿法向投影做 1D 聚类成"层"
  - 在周期水–金属–水体系下标记**两侧**直接面向环境的界面层（每侧固定 1 层，共 2 层）
  - 给界面层附带 `normal_unit`（从金属指向环境一侧）
- `WaterParser.py`
  - 以 MIC 距离 + O–H 截断推断水分子拓扑，输出 `(n_water, 3)` 的 `[O, H1, H2]`
  - 基于氧索引计算：
    - 水质量密度 $\rho(z)$（`g/cm^3`）
    - 取向加权密度（`g/cm^3`，按 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）
    - 指定 c 分数窗口内的 $\theta$ PDF（`degree^-1`）

> 低层 shape/单位/窗口规则以 `history/architecture/modules/data_contract.md` 为准。

### 2) `scripts.structure.Analysis`（多帧/工作流）

输入：`md-pos-*.xyz` 轨迹 + `md.inp`（用于解析正交盒长）

- `Analysis/WaterAnalysis/_common.py`（私有）
  - 从 `md.inp` 解析 `ABC [angstrom]`
  - 用 ASE `iread` 逐帧读取 xyz，并为每帧设置 cell + PBC
  - 每帧自动检测两侧界面位置（沿 c 分数坐标）
  - 计算"选定界面 → 两界面中点"的**半路径**上的密度/取向剖面
  - 将不同帧的剖面重采样到统一归一化网格后做等权系综平均（A 口径）
- `Analysis/WaterAnalysis/WaterDensity.py`
  - 调用 `_common` 计算系综平均密度，并导出 CSV
- `Analysis/WaterAnalysis/WaterOrientation.py`
  - 调用 `_common` 计算系综平均取向加权密度，并导出 CSV
- `Analysis/WaterAnalysis/AdWaterOrientation.py`
  - 从密度剖面自动识别吸附层区间
  - 导出吸附层 profile（密度/取向/掩码）与区间文本
  - 二次遍历轨迹统计吸附层内的 $\theta$ 分布并导出 CSV
- `Analysis/Water.py`
  - 集成三联图入口 `plot_water_three_panel_analysis(...)`
  - 将密度、取向、吸附层与 $\theta$ PDF 组合输出 PNG

## 核心数据流（当前端到端）

1. 解析 `md.inp` 获取正交盒长 $a,b,c$
2. 逐帧读取 xyz → 设置 cell + PBC
3. 逐帧检测两侧界面 c 分数坐标（来自金属界面层）
4. 逐帧识别水分子 → 选定一侧界面起点 → 在半路径上分箱统计
5. 将每帧剖面映射到统一归一化坐标并等权平均
6. 导出 CSV；若走三联图入口，再额外计算吸附层与 $\theta$ 分布并输出 PNG

## 文档镜像与落位

- `history/architecture/modules/scripts/` 必须镜像 `scripts/` 的目录结构，并在每个镜像目录维护：
  - `interface_exposure.md`：对外公开符号/导入方式/兼容承诺
  - `implementation_guidelines.md`：职责边界/依赖方向/契约同步要求
