# 方案 6：将 ASE 导入的 `except Exception` 缩窄为 `except ImportError`

## 对应问题

`problems.md` → 问题 6：`try/except Exception` 仍用于 ASE 可选依赖

## 现状

3 个文件中存在：

```python
try:
    from ase import Atoms
except Exception:  # pragma: no cover
    Atoms = object  # type: ignore
```

- `LayerParser.py`
- `WaterParser.py`
- `_common.py`（同时还有 `iread`）

## 方案：缩窄为 `ImportError`

```python
# 修改后
try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]
```

```python
# _common.py 修改后
try:
    from ase import Atoms
    from ase.io import iread
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]
    iread = None  # type: ignore[assignment]
```

## 涉及文件

| 文件 | 修改 |
|---|---|
| `scripts/structure/utils/LayerParser.py` L27-29 | `Exception` → `ImportError` |
| `scripts/structure/utils/WaterParser.py` L21-23 | `Exception` → `ImportError` |
| `scripts/structure/Analysis/WaterAnalysis/_common.py` L23-28 | `Exception` → `ImportError` |

## 注意

`WaterParser.py` 的 config 导入（L16-19）在 development 分支中**已经修复为直接 import**，无需处理。`LayerParser.py` 的 config 导入（L32）同样已修复。

## 验证

```bash
python -m pytest test/ -v
# 故意在 ase 不可用的环境中测试 import 是否正常 fallback
```
