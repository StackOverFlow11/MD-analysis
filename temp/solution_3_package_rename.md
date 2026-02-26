# 方案 3：包名从 `scripts/` 重命名为 `src/`

## 对应问题

`problems.md` → 问题 3：包名 `scripts/` 不规范

## 方案：重命名为 `src` + 添加 `pyproject.toml`

### 步骤 1：目录重命名

```bash
git mv scripts/ src/
```

重命名后结构：

```
MD-analysis/
├── src/                      # 原 scripts/
│   ├── __init__.py
│   └── structure/
│       ├── __init__.py
│       ├── utils/
│       └── Analysis/
├── test/
├── data_example/
├── history/
├── pyproject.toml            # 新增
└── README.md
```

### 步骤 2：创建 `pyproject.toml`

```toml
[build-system]
requires = ["setuptools>=68.0"]
build-backend = "setuptools.build_meta"

[project]
name = "md-analysis"
version = "0.1.0"
description = "Lightweight analysis utilities for periodic metal-water interfaces"
requires-python = ">=3.10"
dependencies = [
    "numpy",
    "ase",
    "matplotlib",
]

[project.optional-dependencies]
test = ["pytest"]

[tool.setuptools.packages.find]
include = ["src*"]
```

### 步骤 3：全局替换测试文件中的导入路径

包内代码全部使用相对导入（`.` / `..`），**无需修改**。仅测试文件需要批量替换：

| 旧路径 | 新路径 |
|---|---|
| `from scripts.` | `from src.` |
| `import scripts` | `import src` |

### 步骤 4：消除 `sys.path` hack

```python
# test/conftest.py — 删除以下代码
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
```

各测试文件中的 `sys.path` hack 也一并删除。

### 步骤 5：以开发模式安装

```bash
pip install -e .
```

### 步骤 6：更新文档

- `README.md`、`CLAUDE.md` 中所有 `scripts.` 改为 `src.`
- `history/architecture/` 中相关路径更新

## 验证

```bash
pip install -e .
python -m pytest test/ -v
python -c "from src.structure.Analysis import plot_water_three_panel_analysis; print('OK')"
```
