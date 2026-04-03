# CLAUDE.md

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Standard src-layout, Python 3.10+, distributed as `md-analysis` (pip), importable as `md_analysis`.

## Quick Start

```bash
pip install numpy matplotlib ase pytest tqdm
pip install .          # install package (required before running tests)
pytest test/           # run all tests
```

Entry point: `md-analysis` console script → `md_analysis.cli:main` (VASPKIT-style interactive menu).

## Package Map

| 目录 | 用途 | 详情 |
|---|---|---|
| `src/md_analysis/` | 包根：re-export, config, main.py 入口 | → `src/md_analysis/CLAUDE.md` |
| `src/md_analysis/cli/` | 交互式菜单 CLI | → `cli/CLAUDE.md` |
| `src/md_analysis/utils/` | 底层解析器、常量、共享工具 | → `utils/CLAUDE.md` |
| `src/md_analysis/water/` | 水分析工作流 | → `water/CLAUDE.md` |
| `src/md_analysis/electrochemical/` | 电化学分组包（potential + charge） | → `electrochemical/CLAUDE.md` |
| `src/md_analysis/enhanced_sampling/` | 增强采样（SG + TI 准备） | → `enhanced_sampling/CLAUDE.md` |
| `src/md_analysis/scripts/` | 自动化脚本（BaderGen, TIGen） | → `scripts/CLAUDE.md` |
| `context4agent/` | 详细 API 文档（单一真相源） | 见下方 Sync Rules |

## Key Conventions

- **Interface labels:** `"normal_aligned"` (+axis facing) / `"normal_opposed"` (-axis facing)
- **Layer ordering:** `[normal_aligned, slab_interior..., normal_opposed]`
- **Surface normal:** only `"a"/"b"/"c"` cell axes supported (no custom vectors)
- **Frame dirs:** `bader_t*_i*` (Bader) / `potential_t*_i*` (distributed SP) patterns, numerically sorted by `_t(\d+)` regex
- **Output dir structure:** mirrors CLI menu tree via `output_name` on each `MenuNode`:
  - `<outdir>/water/`
  - `<outdir>/electrochemical/potential/<sub>/`
  - `<outdir>/electrochemical/charge/<method>/`
  - `<outdir>/electrochemical/calibration/{fit,predict}/`
  - `<outdir>/enhanced_sampling/slowgrowth/`
  - `<outdir>/enhanced_sampling/constrained_ti/`
- **Relative imports** inside `md_analysis/`; **absolute imports** in tests

## Logging Conventions

- 库代码: `getLogger(__name__)` + NullHandler (PEP 282)
- CLI: `StreamHandler(stderr)` at INFO, 仅在无 handler 时添加
- `verbose` 仅控制 tqdm; log 消息始终发射

## context4agent — Sync Rules (MANDATORY)

`context4agent/` is the authoritative reference for architecture, API contracts, units, and project status. **Every code change MUST synchronize the corresponding `context4agent/` docs.**

| What changed | Update where |
|---|---|
| Module / sub-package added/removed/renamed | `architecture/README.md` + `architecture/modules/src/<module>/` mirror docs |
| Public API (function, signature, export) | corresponding `interface_exposure.md` |
| Implementation pattern / dependency | corresponding `implementation_guidelines.md` |
| Output format (CSV columns, filenames, shapes) | `architecture/modules/data_contract.md` |
| Units / terminology | `architecture/modules/glossary_units.md` |
| Top-level `__all__` | `architecture/modules/src/interface_exposure.md` + `implementation_guidelines.md` |
| Project capabilities / status | `requirements/short_term.md` |

## Workflow Rules

- After completing any module's code changes, always run /commit.
- 始终使用中文回复用户。
