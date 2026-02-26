# 方案 5：清理 `utils` 层公开 API 定位

## 对应问题

`problems.md` → 问题 5：`utils` 层的单帧统计函数未被 `Analysis` 层复用

## 现状（development 分支）

`WaterParser.py` 公开导出 3 个单帧统计函数：

| 函数 | `__all__` 导出 | 被 Analysis 层调用 |
|---|---|---|
| `compute_water_mass_density_z_distribution()` | ✅ | ❌ |
| `compute_water_orientation_weighted_density_z_distribution()` | ✅ | ❌ |
| `compute_water_orientation_theta_pdf_in_c_fraction_window()` | ✅ | ❌ |

`Analysis` 层使用 `_common.py` 的合并计算函数，从未调用上述函数。

## 推荐方案：降级为内部函数

与 main 分支方案一致——将这 3 个函数重命名为 `_` 前缀，从 `__all__` 中移除。

### 变更清单

1. `WaterParser.py`：函数名加 `_` 前缀
2. `utils/__init__.py`：从 `__all__` 和导入中移除
3. `structure/__init__.py`：同步移除
4. 更新 `interface_exposure.md` 文档

## 验证

```bash
python -m pytest test/ -v
grep -rn "compute_water_mass_density_z_distribution" --include="*.py"
# 确认仅在 WaterParser.py 内部使用
```
