# md_analysis.electrochemical — Interface Exposure

> 对应代码：`src/md_analysis/electrochemical/__init__.py`

## 角色

分组命名空间，将电化学相关分析（电势、电荷）聚合在同一父包下。自身不包含业务逻辑。

## Public API

| Symbol | Source | Description |
|---|---|---|
| `potential` | sub-package | 电势分析工作流（re-export） |
| `charge` | sub-package | Bader 电荷分析工作流（re-export） |
| `calibration` | sub-package | σ→φ 标定映射（re-export） |

```python
__all__ = ["potential", "charge", "calibration"]
```

## 推荐导入方式

```python
# 直接访问子包（推荐）
from md_analysis.electrochemical import potential
from md_analysis.electrochemical import charge

# 顶层便捷访问（向后兼容）
from md_analysis import potential
from md_analysis import charge
```

## 稳定性

- **Stable** — `potential` 和 `charge` 的 re-export 为稳定契约
- 顶层 `md_analysis` 通过 `from .electrochemical import potential, charge` 保持向后兼容
