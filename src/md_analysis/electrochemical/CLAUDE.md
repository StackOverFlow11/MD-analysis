# electrochemical — 开发备忘

## 定位

纯命名空间聚合包，不含任何分析逻辑。`__init__.py` 仅 re-export `potential` 和 `charge` 子包。

## 子目录

| 目录 | 用途 |
|---|---|
| `potential/` | Hartree 势 / 电极电位分析 → `potential/CLAUDE.md` |
| `charge/` | Bader 电荷表面电荷密度 → `charge/CLAUDE.md` |
