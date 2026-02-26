# 方案 1：将跨模块使用的私有函数提升为受控共享接口

## 对应问题

`problems.md` → 问题 1：跨模块导入私有（`_` 前缀）函数

## 现状

`_common.py` 和 `AdWaterOrientation.py` 从 `utils` 层导入多个 `_` 前缀函数：

| 被导入的私有函数 | 定义位置 | 导入方 |
|---|---|---|
| `_circular_mean_fractional` | `LayerParser.py` | `_common.py` |
| `_oxygen_to_hydrogen_map` | `WaterParser.py` | `_common.py`、`AdWaterOrientation.py` |
| `_compute_bisector_cos_theta_vec` | `WaterParser.py` | `_common.py`、`AdWaterOrientation.py` |
| `_theta_bin_count_from_ndeg` | `WaterParser.py` | `AdWaterOrientation.py` |

另外 `Water.py` 导入 `_common.py` 的 `_compute_density_orientation_ensemble`。

## 方案：在 `utils/__init__.py` 中新增受控内部导出区

不改变函数名（仍保留 `_` 前缀），但在 `utils/__init__.py` 中**显式导入并声明**这些函数为"内部共享接口"，让依赖关系从"偷偷跨模块导入私有函数"变为"通过包接口获取受控共享的内部工具"。

### 步骤 1：在 `utils/__init__.py` 增加内部共享导出

```python
# scripts/structure/utils/__init__.py

# --- Internal shared helpers (used by Analysis layer, NOT part of public API) ---
# These are intentionally excluded from __all__ but explicitly imported here
# to centralize cross-layer internal dependencies.
from .LayerParser import _circular_mean_fractional
from .WaterParser import _compute_bisector_cos_theta_vec
from .WaterParser import _oxygen_to_hydrogen_map
from .WaterParser import _theta_bin_count_from_ndeg
```

`__all__` 不变——这些函数对外部用户仍然不可见，但包内上层模块可以通过 `from ...utils import _circular_mean_fractional` 引用，依赖路径变得清晰可控。

### 步骤 2：修改 `_common.py` 的导入来源

```python
# 修改前
from ...utils.LayerParser import SurfaceGeometryError, _circular_mean_fractional, detect_interface_layers
from ...utils.WaterParser import (
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
    ...
)

# 修改后
from ...utils import (
    SurfaceGeometryError,
    detect_interface_layers,
    _circular_mean_fractional,
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
)
from ...utils.WaterParser import (   # 公开符号 + 常量
    AVOGADRO_NUMBER,
    ANGSTROM3_TO_CM3,
    WaterTopologyError,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)
```

### 步骤 3：修改 `AdWaterOrientation.py` 同理

```python
# 修改前
from ...utils.WaterParser import (
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
    _theta_bin_count_from_ndeg,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)

# 修改后
from ...utils import (
    _compute_bisector_cos_theta_vec,
    _oxygen_to_hydrogen_map,
    _theta_bin_count_from_ndeg,
    detect_water_molecule_indices,
    get_water_oxygen_indices_array,
)
```

### 步骤 4：修改 `Water.py` 的 `_common` 导入

`Water.py` 当前直接穿透到 `WaterAnalysis._common`：

```python
# 修改前
from .WaterAnalysis._common import StartInterface, _compute_density_orientation_ensemble
```

改为通过 `WaterAnalysis/__init__.py` 导出：

```python
# WaterAnalysis/__init__.py 新增
from ._common import StartInterface
from ._common import _compute_density_orientation_ensemble

# Water.py 修改后
from .WaterAnalysis import StartInterface, _compute_density_orientation_ensemble
```

## 收益

- 所有跨层内部依赖都经过 `__init__.py` 中转，**一目了然**
- 如果底层重构，只需改 `utils/__init__.py` 的导入映射，上层无感
- 不影响公开 API（`__all__` 不变）

## 验证

```bash
python -m pytest test/ -v
```
