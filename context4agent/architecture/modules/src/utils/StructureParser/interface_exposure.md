# `md_analysis.utils.StructureParser` 接口暴露约定

> 对应代码：`src/md_analysis/utils/StructureParser/__init__.py`
>
> 本子包封装结构解析相关模块：周期聚类、金属层识别、水分子拓扑。

## 1. 接口角色

- `StructureParser` 是 `utils` 的内部子包，其公开符号通过 `utils/__init__.py` 门面层 re-export。
- 外部消费者应通过 `from md_analysis.utils import xxx` 导入，不应直接依赖子包路径。
- 内部消费者（如 `charge.BaderAnalysis`、`water.WaterAnalysis._common`）可直接导入子包路径。

## 2. 模块组成

| 模块 | 职责 |
|---|---|
| `ClusterUtils.py` | 1D 周期性聚类、最大间隙检测、间隙中点计算 |
| `LayerParser.py` | 金属层识别、界面层标记、法向符号判定 |
| `WaterParser.py` | 水分子拓扑标记、z 轴密度/取向分布统计 |

## 3. `__all__` 导出清单

```python
"cluster_1d_periodic", "find_largest_gap_periodic", "gap_midpoint_periodic",
"Layer", "SurfaceDetectionResult", "SurfaceGeometryError",
"detect_interface_layers", "format_detection_summary",
"WaterTopologyError", "detect_water_molecule_indices", "get_water_oxygen_indices_array",
```

## 4. 内部依赖

- `LayerParser` → `ClusterUtils`（同级导入）
- `LayerParser` → `../config.py`（`DEFAULT_METAL_SYMBOLS`）
- `WaterParser` → `../config.py`（`DEFAULT_THETA_BIN_DEG` 等）
