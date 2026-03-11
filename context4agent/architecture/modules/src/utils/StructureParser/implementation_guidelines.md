# `md_analysis.utils.StructureParser` 内部实现准则

> 适用范围：`src/md_analysis/utils/StructureParser/`（`ClusterUtils.py`、`LayerParser.py`、`WaterParser.py`）。

## 1. 子包职责

本子包将原来平铺在 `utils/` 下的三个结构解析模块组织到一起，按"结构解析"职责分组：

- **ClusterUtils**：纯数值的 1D 周期性聚类算法，不依赖 ASE
- **LayerParser**：基于 ClusterUtils 的金属层识别与界面标记，依赖 ASE Atoms
- **WaterParser**：水分子拓扑推断与 z 轴分布统计，依赖 ASE Atoms

## 2. 导入约定

- 子包内模块间使用同级相对导入（`from .ClusterUtils import ...`）
- 引用 `utils/config.py` 使用上跳一级（`from ..config import ...`）
- 子包不得直接导入 `water/`、`potential/`、`charge/` 等上层模块

## 3. 实现准则

各模块的具体实现准则（坐标口径、分箱规则、异常分类等）继续遵循父级 `utils/implementation_guidelines.md` 中的对应章节。

### 3.1 `ClusterUtils.py` 补充

- `_circular_mean(values, period)` 对空输入做显式校验，抛出 `ValueError("cannot compute circular mean of empty array")`。

### 3.2 `LayerParser.py` 补充

- `circular_mean_fractional()` 已提升为公开 API（去除 `_` 前缀），委托给 `ClusterUtils._circular_mean(values, period=1.0)`，避免重复实现。
- 轴映射使用 `config.py` 中的集中式常量 `AXIS_MAP`（取代原有的模块局部 `_AXIS_MAP`）。
- 界面标签使用 `config.py` 中的 `INTERFACE_NORMAL_ALIGNED` / `INTERFACE_NORMAL_OPPOSED` 常量（取代硬编码字符串）。
