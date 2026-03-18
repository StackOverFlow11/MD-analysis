# md_analysis.electrochemical — Implementation Guidelines

## 角色

纯分组包，将电化学分析（`potential`、`charge`、`calibration`）组织在同一命名空间下。

## `__init__.py` 准则

- 仅允许：`from . import calibration, charge, potential` + `__all__` 声明
- 禁止：数值计算、文件 I/O、配置读取、三方库导入

## 依赖方向

- `electrochemical.potential` → `utils`（CubeParser, LayerParser, ClusterUtils, config, _io_helpers）
- `electrochemical.charge` → `utils`（BaderParser, LayerParser, WaterParser, config, _io_helpers）
- `electrochemical.calibration` → numpy, matplotlib（scipy 仅 spline 拟合时 lazy import）
- `potential`、`charge`、`calibration` **互不依赖**
- 其他模块不应直接依赖 `electrochemical/__init__.py`（应通过子包访问）

## 变更流程

新增/删除/重命名 `electrochemical/` 下的子包时，必须同步：
1. `src/md_analysis/electrochemical/__init__.py` 的导入与 `__all__`
2. `src/md_analysis/__init__.py` 的 re-export（若需要顶层便捷访问）
3. 本目录及对应子目录的 `interface_exposure.md` / `implementation_guidelines.md`
