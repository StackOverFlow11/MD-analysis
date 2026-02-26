# 架构评审问题清单（development 分支）

> 基于 `development` 分支代码重新评审。  
> 相比 `main` 分支，新版已引入 `_common.py` 统一轨迹 I/O，修复了 `DEFAULT_OUTPUT_DIR`，消除了大部分重复代码。

---

## 🔴 高优先级

### 1. 跨模块导入私有（`_` 前缀）函数

`_common.py` 和 `AdWaterOrientation.py` 大量导入 `utils` 层的 `_` 前缀内部函数：

```python
# _common.py
from ...utils.LayerParser import _circular_mean_fractional
from ...utils.WaterParser import _compute_bisector_cos_theta_vec, _oxygen_to_hydrogen_map

# AdWaterOrientation.py
from ...utils.WaterParser import _compute_bisector_cos_theta_vec, _oxygen_to_hydrogen_map, _theta_bin_count_from_ndeg
```

此外 `Water.py` 直接导入 `_common` 模块的私有函数：

```python
# Water.py
from .WaterAnalysis._common import StartInterface, _compute_density_orientation_ensemble
```

**问题**：`_` 前缀函数按 Python 惯例属于"模块内部实现"，跨模块依赖私有符号会破坏封装性，导致上层代码与底层实现细节强耦合。如果 `utils` 层重构内部函数名或签名，上层会直接报错。

**建议**：将被跨模块使用的函数提升为受控的内部共享接口（去掉 `_` 前缀并收入统一的共享模块），或通过 `__init__.py` 显式导出。

### 2. `WaterDensity.py` 和 `WaterOrientation.py` 单独调用时浪费计算

两者都调用 `_compute_density_orientation_ensemble()`，该函数**同时**计算密度和取向，但各自只用了一半结果：

```python
# WaterDensity.py — 只用 rho_ensemble，丢弃 orient_ensemble
common_centers_u, mean_path_A, rho_ensemble, _ = _compute_density_orientation_ensemble(...)

# WaterOrientation.py — 只用 orient_ensemble，丢弃 rho_ensemble
common_centers_u, mean_path_A, _, orient_ensemble = _compute_density_orientation_ensemble(...)
```

如果用户只需要密度分析，仍然会执行全部取向计算（包括 `find_mic` + 向量化 bisector），反之亦然。

**建议**：为 `_single_frame_density_and_orientation` 增加 `compute_orientation: bool` 标志，或拆分为独立的单帧函数。

---

## 🟡 中优先级

### 3. 包名 `scripts/` 不规范

架构文档建议使用 `src/`，但实际仍为 `scripts/`：

- 语义误导（`scripts` 通常暗示一次性脚本，不是可复用库）
- 测试中仍需要 `sys.path` hack：

```python
# test/conftest.py
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
```

**建议**：重命名为 `src/` + 添加 `pyproject.toml`，消除 `sys.path` hack。

### 4. `conftest.py` 中 `parse_abc_from_md_inp` 仍然重复

`_common.py` 已有 `_parse_abc_from_md_inp`，但 `test/conftest.py` 仍独立维护了一份相同实现。

**建议**：测试中直接从 `_common` 导入，或将其提升为 `utils` 层的公开/受控共享函数。

---

## 🟢 低优先级

### 5. `utils` 层的单帧统计函数未被 `Analysis` 层复用

`WaterParser.py` 公开导出了 3 个单帧统计函数（`compute_water_mass_density_z_distribution` 等），但 `Analysis` 层完全使用 `_common.py` 的界面映射逻辑，从未调用这些函数。

它们在 `__all__` 中导出，增加维护成本，也容易让用户混淆该调哪个。

**建议**：明确其适用场景（全 cell z 轴分箱），或降级为内部函数。

### 6. `try/except Exception` 仍用于 ASE 可选依赖

`LayerParser.py`、`WaterParser.py`、`_common.py` 中：

```python
try:
    from ase import Atoms
except Exception:  # pragma: no cover
    Atoms = object  # type: ignore
```

`Exception` 过于宽泛，应缩窄为 `ImportError`。

---

## ✅ 已在 development 分支修复的问题（存档）

| 原问题 | 修复方式 |
|---|---|
| 轨迹被重复读取 5 次 | `_common.py` 将密度+取向合并为单次遍历；`plot_water_three_panel_analysis` 从 5 次降到 2 次 |
| `DEFAULT_OUTPUT_DIR = Path.cwd()` 导入时求值 | 改为 `None`，各函数内 `Path.cwd()` 运行时 resolve |
| `_circular_mean_fractional` 在 `WaterDensity.py` 中重复 | `_common.py` 从 `LayerParser` 导入 |
| `_oxygen_to_hydrogen_map` 在 3 处重复 | `_common.py` 和 `AdWaterOrientation.py` 从 `WaterParser` 导入 |
| `_parse_abc_from_md_inp` 在 `WaterDensity.py` 中重复 | 统一到 `_common.py` |
| `detect_interface_layers` 的 `n_interface_layers` 废弃参数 | 已移除 |
| `LayerParser.py` config 导入用 `try/except` 兜底 | 改为直接 `from .config import` |

---

## 📋 建议实施顺序

```
3 → 4 → 1 → 6 → 2 → 5
```

| 步骤 | 问题编号 | 内容 | 理由 |
|---|---|---|---|
| 1st | **3** | 包名 `scripts/` → `src/` + `pyproject.toml` | 基础设施先行，后续改动在正确包结构上进行 |
| 2nd | **4** | 消除 `conftest.py` 中的重复 `parse_abc_from_md_inp` | 配合方案 3 一起做，消除 `sys.path` hack |
| 3rd | **1** | 跨模块私有函数导入 → 提升为受控共享接口 | 消除封装性问题，为未来重构打基础 |
| 4th | **6** | ASE 导入 `Exception` → `ImportError` | 独立小改动，随时可做 |
| 5th | **2** | 为合并计算增加跳过标志 | 优化性能，影响范围小 |
| 6th | **5** | 清理 `utils` 层公开 API 定位 | 影响最小，最后收尾 |

> **原则**：先整包结构（3→4），再修封装（1→6），最后优化和清理（2→5）。
