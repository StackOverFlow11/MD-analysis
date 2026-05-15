# 目标模块职责边界

> 用途：给后续 agent 提供稳定的架构边界说明，避免重构计划执行时职责漂移。
> 本文件不描述具体迁移步骤；迁移顺序见 `overall_reconstruction_plan.md`。

## 总原则

- 先清理底层 `utils`，再设计 `engines` 数据接口，再迁移业务层和 CLI。
- 业务层不应直接解析 CP2K / VASP / Bader / cube / stdout / restart 等原始文件。
- 文件格式解析先进入 `utils.formats`，多文件组合与 engine-neutral 数据结构由 `engines` 提供。
- `structure` 可以处理解析后的 `ase.Atoms` / arrays / dataclass 字段，但不自己 import CP2K/VASP parser。

## `utils`

底层工具层，只提供可复用基础能力。

- `utils.formats`：单文件 / 单格式 parser。
  - `formats.cp2k`：CP2K 文件解析，如 cell、COLVAR restart、LagrangeMultLog、stdout、xyz。
  - `formats.vasp`：VASP 文件解析或 placeholder，如 REPORT、OUTCAR、LOCPOT。
  - `formats.common`：真正通用格式，如纯 cube reader。
  - `formats.bader`：Bader/ACF/POTCAR 后处理格式。
- `utils.io`：真正通用 I/O helper，如 CSV 写出、通用文件系统辅助；不承载 CP2K/VASP 目录语义。
- `utils.structure`：engine-neutral 结构/几何/拓扑工具，如 layer、water、cluster；只吃结构化对象，不读原始文件。
- `utils.constants`：物理常量、单位换算、默认参数。

示范性目录结构（仅供后续计划参考，不是硬性要求）：

```text
src/md_analysis/utils/
  __init__.py
  constants.py
  formats/
    __init__.py
    common/
      __init__.py
      cube.py
    cp2k/
      __init__.py
      cell.py
      colvar.py
      stdout.py
      xyz.py
    vasp/
      __init__.py
      report.py
      outcar.py
      locpot.py
    bader/
      __init__.py
      acf.py
      potcar.py
  io/
    __init__.py
    csv.py
    discovery.py
  structure/
    __init__.py
    layer.py
    water.py
    cluster.py
```

## `engines`

计算引擎适配层，负责把 engine-specific 文件组合成上层可用的 engine-neutral 数据结构。

- `models`：稳定数据结构，例如 `ConstraintMetadata`、`LambdaSeries`、`PotentialFrame`、`FermiRecord`，未来可扩展 `ChargeFrame` / `ChargeTrajectory`。
- `protocols`：按数据契约抽象接口，不按 workflow 功能抽象 parser。
  - 例如 `ConstraintMDParser`、未来可能的 `PotentialFrameReader`、`ChargeTrajectoryReader`。
- `registry`：未来如出现多个数据接口，可考虑泛型 `ParserRegistry[T]`。
- `engines.cp2k`：CP2K facade / adapter，可拆成 package，但外部仍应能从 `md_analysis.engines.cp2k` 导入稳定入口。
- `engines.vasp`：VASP adapter；placeholder 未实现前不得自动注册为可用 parser。

示范性目录结构（仅供后续计划参考，不是硬性要求）：

```text
src/md_analysis/engines/
  __init__.py
  models.py
  protocols.py
  registry.py
  cp2k/
    __init__.py
    constraint.py
    potential.py
    charge.py
    common.py
  vasp/
    __init__.py
    constraint.py
    potential.py
    charge.py
```

## 业务模块

业务模块只做物理/化学分析，不直接理解原始文件格式。

- `water`：水密度、取向、吸附层等分析；依赖 `utils.structure` 和结构化轨迹数据。
- `electrochemical`：电势、电荷、标定等分析；依赖 `engines` 提供的电势帧、电荷轨迹或其他结构化数据。
- `enhanced_sampling`：slow-growth / constrained TI 等分析；依赖 `engines` 提供的约束 MD metadata 和 lambda series。

## 入口层

- `workflows`：程序化入口，组织业务模块调用，返回 `WorkflowResult`。
- `agent`：非交互式 agent/MCP 入口；不放业务逻辑。
- `cli`：交互式菜单入口；最后迁移，尽量只做参数收集和 workflow 调用。

## 依赖方向

```text
cli / agent
  -> workflows
  -> business modules
  -> engines
  -> utils.formats / utils.io / utils.structure / utils.constants
```

允许业务模块直接使用：

- `utils.structure`
- `utils.constants`
- 真正通用的 `utils.io`

不鼓励业务模块直接使用：

- `utils.formats.cp2k`
- `utils.formats.vasp`
- `utils.formats.bader`
- 带 CP2K/VASP 文件命名或目录语义的 parser/helper
