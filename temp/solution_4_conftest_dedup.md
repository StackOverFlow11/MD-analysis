# 方案 4：消除 `conftest.py` 中的重复 `parse_abc_from_md_inp`

## 对应问题

`problems.md` → 问题 4：`conftest.py` 中 `parse_abc_from_md_inp` 仍然重复

## 现状

`_common.py` 已有 `_parse_abc_from_md_inp`，但 `test/conftest.py` 仍独立维护一份相同实现（L15-22）。

## 方案

### 选项 A：测试直接导入 `_common` 的版本（推荐）

```python
# test/conftest.py（修改后）
"""Shared pytest fixtures and helpers for MD Analysis tests."""

from __future__ import annotations

# parse_abc_from_md_inp 不再在此定义，直接从源码导入
from src.structure.Analysis.WaterAnalysis._common import _parse_abc_from_md_inp as parse_abc_from_md_inp
```

**前提**：方案 3（包重命名 + `pip install -e .`）已完成。

### 选项 B：将函数提升为 `utils` 层的公开工具

如果认为"从 `_common` 导入"仍然是依赖私有模块，可以将 `parse_abc_from_md_inp` 提升到 `utils/` 层并加入 `__all__`，作为通用的 CP2K I/O 工具：

```python
# scripts/structure/utils/io.py（新建）
def parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    ...
```

然后 `_common.py` 和 `conftest.py` 都从 `utils.io` 导入。

## 建议

如果近期没有更多 CP2K I/O 需求，选项 A 最简单。如果计划扩展 I/O 功能（如 ener 文件解析），建议直接选项 B 建立 `utils/io.py`。

## 验证

```bash
python -m pytest test/ -v
```
