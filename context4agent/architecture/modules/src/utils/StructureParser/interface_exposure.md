# `md_analysis.utils.StructureParser` 接口暴露约定

> 对应代码：`src/md_analysis/utils/StructureParser/__init__.py`
>
> 本子包封装结构解析相关模块：周期聚类、金属层识别、水分子拓扑。

## 1. 接口角色

- `StructureParser` 是 `utils` 的结构解析子包。
- 所有消费者（外部与内部）均通过 **模块直接路径** 导入：
  - `from md_analysis.utils.StructureParser.LayerParser import detect_interface_layers`
  - `from md_analysis.utils.StructureParser.WaterParser import detect_water_molecule_indices`
  - `from md_analysis.utils.StructureParser.ClusterUtils import cluster_1d_periodic`
- `utils/__init__.py` 不再 re-export；`StructureParser/__init__.py` 也不再集中 re-export。

## 2. 模块组成

| 模块 | 职责 |
|---|---|
| `ClusterUtils.py` | 1D 周期性聚类、最大间隙检测、间隙中点计算 |
| `LayerParser.py` | 金属层识别、界面层标记、法向符号判定 |
| `WaterParser.py` | 水分子拓扑标记、z 轴密度/取向分布统计 |

## 3. 公开符号（按模块）

**ClusterUtils.py**：`cluster_1d_periodic`、`find_largest_gap_periodic`、`gap_midpoint_periodic`

**LayerParser.py**：`Layer`、`SurfaceDetectionResult`、`SurfaceGeometryError`、`circular_mean_fractional`、`detect_interface_layers`、`format_detection_summary`、`mic_delta_fractional`

**WaterParser.py**：`WaterTopologyError`、`detect_water_molecule_indices`、`get_water_oxygen_indices_array`（`_compute_bisector_cos_theta_vec`、`_oxygen_to_hydrogen_map`、`_theta_bin_count_from_ndeg` 为 cross-layer 内部 helper，不稳定）

## 4. 内部依赖

- `LayerParser` → `ClusterUtils`（同级导入：`cluster_1d_periodic`、`find_largest_gap_periodic`、`_circular_mean`）
- `LayerParser` → `../config.py`（`DEFAULT_METAL_SYMBOLS`、`AXIS_MAP`、`INTERFACE_NORMAL_ALIGNED`、`INTERFACE_NORMAL_OPPOSED`）
- `WaterParser` → `../config.py`（`DEFAULT_THETA_BIN_DEG` 等）
