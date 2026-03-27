# Enhanced Sampling CLI 重构 + Constrained TI 入口设计

> 状态：审阅后修订（2026-03-24）

## 背景

`constrained_ti/` 模块已实现完整的收敛诊断和自由能积分功能，但 CLI 中仅有 SG 绘图（301-302）和 TI 工作目录生成（421-422），**缺少 TI 数据分析入口**。同时 SG 的 301/302 目前直接挂在 Enhanced Sampling 组下，缺少子组层级。

## 菜单编号方案

将 Enhanced Sampling（菜单 3）拆为两个子组：`30x` = Slow-Growth，`31x` = Constrained TI。现有 301/302 移入 `30)` 子组，新增 311/312。

### 变更前

```
3) Enhanced Sampling
   301  Slow-Growth Quick Plot
   302  Slow-Growth Publication Plot
```

### 变更后

```
3) Enhanced Sampling
   30)  Slow-Growth
        301  Quick Plot                      ← 已有，移入子组
        302  Publication Plot                ← 已有，移入子组
   31)  Constrained TI Analysis             ← 新增子组
        311  Single-Point Diagnostics        ← 新增
        312  Full TI Analysis                ← 新增
```

### build_menu_tree() 变更

```python
# --- Enhanced Sampling ---
enhanced = MenuGroup("3", "Enhanced Sampling",
                     output_name="enhanced_sampling")

slowgrowth = MenuGroup("30", "Slow-Growth",
                       output_name="slowgrowth")
slowgrowth.add(
    SGQuickPlotCmd("301", "Quick Plot"),
    SGPublicationPlotCmd("302", "Publication Plot"),
)

constrained_ti = MenuGroup("31", "Constrained TI Analysis",
                           output_name="constrained_ti")
constrained_ti.add(
    TISingleDiagCmd("311", "Single-Point Diagnostics"),
    TIFullAnalysisCmd("312", "Full TI Analysis"),
)

enhanced.add(slowgrowth, constrained_ti)
```

### 用户交互示例

```
---------- MD-Analysis ----------
 3) Enhanced Sampling
 Input: 3

---------- Enhanced Sampling ----------
 30) Slow-Growth
 31) Constrained TI Analysis
  0) Back / Exit
 Input: 31

---------- Constrained TI Analysis ----------
 311  Single-Point Diagnostics
 312  Full TI Analysis
  0) Back / Exit
 Input: 312
```

flat index 机制保证用户也可从根直接输入 `311` 或 `312` 跳转。

### 编号分配理由

- `30x` / `31x` 平级，语义清晰
- SG 的 301/302 编号不变，仅多了一层父节点，**用户习惯不受影响**（flat index 仍可直达）
- 留出 `303-309`、`313-319` 供未来扩展

## 输出目录

遵循 `output_name` 机制，父链自动拼接：

| 命令 | 父链 output_name | 完整路径 |
|---|---|---|
| 301/302 | enhanced_sampling → slowgrowth | `<outdir>/enhanced_sampling/slowgrowth/` ← 不变 |
| 311 | enhanced_sampling → constrained_ti | `<outdir>/enhanced_sampling/constrained_ti/` |
| 312 | enhanced_sampling → constrained_ti | `<outdir>/enhanced_sampling/constrained_ti/` |

301/302 的实际输出路径不变（原来 `MenuCommand.output_name = "slowgrowth"`，现在改由父 `MenuGroup("30", output_name="slowgrowth")` 提供，命令自身不再需要 `output_name`）。

> **⚠️ 实施警告**：**必须删除** `_SlowgrowthPlotCmd` 类的 `output_name = "slowgrowth"` 属性（`_enhanced_sampling.py:84`）。否则父链变为 `self("slowgrowth") → group("slowgrowth") → enhanced("enhanced_sampling")`，路径错误拼接为 `enhanced_sampling/slowgrowth/slowgrowth`。

311/312 共享 `output_name`，因为输出文件名不冲突（单点 `ti_diag_xi*.png` + `ti_single_point.csv`，多点 `ti_convergence_report.csv` + `ti_free_energy.*`）。

---

## 311 — Single-Point Diagnostics

### 用途

对单个约束 MD 目录做四步收敛诊断（ACF、F&P block average、running average、Geweke），输出诊断图 + CSV。

### 典型目录结构（即 `data_example/ti/ti_target_-0.781213/` 中的内容）

```
ti_target_-0.781213/
  cMD-1.restart                        # 主 restart（连字符 -1，非检查点）
  cMD-1_1000.restart                   # 中间检查点（下划线 _1000，被发现逻辑排除）
  cMD-constraint_force.dat-1.LagrangeMultLog
  cMD.inp
  init.xyz
  ...
```

### 参数采集方式

覆写 `_collect_all_params()`（与 SG 命令模式一致），不使用 `params`/`advanced_params` 元组。原因：需先自动发现文件再提供默认值，且 SEM target 为可空 float（现有 `FloatParam` 不支持 `None`）。

### 参数采集流程

```
Input: 311

  [自动发现 .restart 和 .LagrangeMultLog]
  Found restart: cMD-1.restart
  Found log:     cMD-constraint_force.dat-1.LagrangeMultLog

  COLVAR restart file [/path/to/cMD-1.restart]:
  LagrangeMultLog file [/path/to/cMD-...LagrangeMultLog]:
  Equilibration frames to discard [0]:
  SEM target (a.u., empty=none):
  Colvar ID (empty=primary):
  Output directory [analysis]:
```

### 参数详情

| 参数 | 键常量 | 类型 | 默认值 | 说明 |
|---|---|---|---|---|
| restart file | `K.RESTART_PATH` | str | 自动发现 | 复用 SG 的 `_discover_restart_file` / `_discover_log_file` |
| log file | `K.LOG_PATH` | str | 自动发现 | 同上 |
| equilibration | `K.EQUILIBRATION` | int | 0 | 丢弃的初始帧数 |
| SEM target | `K.SEM_TARGET` | float\|None | None | 精度目标（a.u.），无则不判 pass/fail。用 `prompt_str` + 手动 `float()` 转换实现可空 |
| colvar_id | `K.COLVAR_ID` | int\|None | None | 默认用 primary CV |
| output_dir | `K.OUTDIR` | str | "analysis" | |

### 调用的底层 API

```python
from md_analysis.enhanced_sampling.constrained_ti.workflow import standalone_diagnostics

result = standalone_diagnostics(
    restart_path=ctx[K.RESTART_PATH],
    log_path=ctx[K.LOG_PATH],
    equilibration=ctx[K.EQUILIBRATION],
    sem_target=ctx[K.SEM_TARGET],
    colvar_id=ctx[K.COLVAR_ID],
    output_dir=ctx[K.OUTDIR_RESOLVED],
)
```

### 输出文件

| 文件 | 说明 |
|---|---|
| `ti_diag_xi{value}.png` | 2×2 诊断图（running avg / ACF / F&P / summary） |
| `ti_single_point.csv` | 单行诊断报告 |

> **注意**：PNG 文件名由 `plot_point_diagnostics` 生成，格式为 `ti_diag_xi{xi:.4f}.png`（如 `ti_diag_xi-0.7812.png`），不是 `ti_single_point_diag.png`。`config.py` 中的 `DEFAULT_STANDALONE_PNG_NAME` 常量当前未被 `standalone_diagnostics` 使用。

### 终端摘要（execute 结束时 print）

```
  ξ = -0.781213 a.u.
  N = 1000 frames (after discarding 0)
  ⟨λ⟩ = -0.003245 a.u.
  τ_corr = 12.3 frames, N_eff = 40.7
  SEM_final = 0.000512 a.u. (method: plateau)
  Geweke |z| = 0.84  (PASS)
  Overall: PASS / FAIL / N/A

  diagnostics_png: /path/to/ti_diag_xi-0.7812.png
  csv:             /path/to/ti_single_point.csv
```

---

## 312 — Full TI Analysis

### 用途

扫描包含多个 `ti_target_*` 子目录的根目录，对每个约束点做收敛诊断，然后梯形积分得到自由能曲线 ΔA(ξ)。

### 典型目录结构（即 `data_example/ti/` 整体）

```
ti_work/                          # ← 用户指定的根目录
  ti_target_-0.781213/
    cMD-1.restart
    cMD-constraint_force.dat-1.LagrangeMultLog
    ...
  ti_target_-0.510226/
    ...
  ti_target_-0.239239/
    ...
  ...（更多约束点）
```

### 参数采集方式

同 311，覆写 `_collect_all_params()`。`TI_DIR_PATTERN` 使用 `prompt_str` + 手动验证（仅允许 `"ti_target"/"xi"/"auto"`），无效输入时提示重试。

### 参数采集流程

```
Input: 312

  TI root directory [.]:
  Directory pattern (ti_target/xi/auto) [auto]:
  Equilibration frames to discard (all points) [0]:
  Free-energy tolerance ε (eV) [0.05]:
  Output directory [analysis]:
```

### 参数详情

| 参数 | 键常量 | 类型 | 默认值 | 说明 |
|---|---|---|---|---|
| root_dir | `K.TI_ROOT_DIR` | str | "." | 包含 ti_target_* 子目录的根路径 |
| pattern | `K.TI_DIR_PATTERN` | str | "auto" | 目录发现模式，仅 `"ti_target"/"xi"/"auto"` 有效 |
| equilibration | `K.EQUILIBRATION` | int | 0 | 默认丢弃帧数（作为逐点的默认值） |
| epsilon_tol_ev | `K.EPSILON_TOL_EV` | float | 0.05 | 自由能精度容差（eV） |
| output_dir | `K.OUTDIR` | str | "analysis" | |

> **equilibration 按点设置**：发现约束点后，CLI 会询问 `"Set per-point equilibration frames? (y/N)"`。选 Yes 则逐点提示输入（回车使用默认值）；选 No 则所有点使用统一默认值。底层 `analyze_ti()` 接收 `int | list[int]`。

### execute 内部流程

```python
from md_analysis.enhanced_sampling.constrained_ti.io import discover_ti_points, load_ti_series
from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    analyze_ti, write_convergence_csv, write_free_energy_csv,
)
from md_analysis.enhanced_sampling.constrained_ti.plot import (
    plot_point_diagnostics, plot_free_energy_profile,
)

# 1. 发现约束点
point_defs = discover_ti_points(Path(root_dir), pattern=pattern)
print(f"  Found {len(point_defs)} constraint points")
for p in point_defs:
    print(f"    ξ = {p.xi:.6f}")

# 1.5 逐点 equilibration（可选）
default_equil = ctx[K.EQUILIBRATION]
if len(point_defs) > 1 and prompt_bool("Set per-point equilibration frames?", default=False):
    equilibration = [prompt_int(f"  ξ={p.xi:.6f} ...", default=default_equil) for p in point_defs]
else:
    equilibration = default_equil

# 2. 加载序列
series_data = load_ti_series(point_defs)
xi_values = np.array([x for x, _, _ in series_data])
lambda_list = [s for _, s, _ in series_data]

# 3. dt 一致性校验
dts = [dt for _, _, dt in series_data]
if len(set(f"{d:.10f}" for d in dts)) > 1:
    print(f"  WARNING: inconsistent timesteps detected: {set(dts)}")
    print(f"  Using dt = {dts[0]:.6f} fs from first point")
dt = dts[0]

# 4. 分析
ti_report = analyze_ti(
    xi_values, lambda_list, dt,
    epsilon_tol_ev=epsilon_tol_ev,
    equilibration=equilibration,
)

# 5. 输出
write_convergence_csv(ti_report, output_dir=outdir)
write_free_energy_csv(ti_report, output_dir=outdir)
plot_free_energy_profile(ti_report, output_dir=outdir)

# 6. 所有点均输出诊断图
for r in ti_report.point_reports:
    plot_point_diagnostics(r, output_dir=outdir)
```

> **诊断图策略**：默认对所有点生成 `ti_diag_xi*.png`，方便直接对比各点收敛质量，无需再单独调用 311。

### 输出文件

| 文件 | 说明 |
|---|---|
| `ti_free_energy.png` | 双轴图：左=dA/dξ 误差棒，右=A(ξ) 累积曲线 |
| `ti_free_energy.csv` | ξ, weight, dA/dξ, SEM, A(eV), σ_A(eV) |
| `ti_convergence_report.csv` | 每点一行的诊断汇总 |
| `ti_diag_xi{value}.png` | 所有点的 2×2 诊断图（K 张） |

### 终端摘要

```
  Found 8 constraint points
  Analyzing: ξ = -0.7812, -0.5102, -0.2392, 0.0314, 0.3024, 0.5733, 0.8440, 1.1149
  Timestep: 1.000000 fs

  Point  ξ          ⟨λ⟩         SEM          Status
  ─────────────────────────────────────────────────
  0      -0.781213  -0.003245   0.000512     PASS
  1      -0.510226   0.001123   0.000834     PASS
  ...
  7       1.114938   0.000021   0.000015     PASS

  ΔA = 0.1234 ± 0.0056 eV
  Status: ALL PASS / 2 FAILED (indices: 3, 5)

  [如有失败点]
  Suggested time allocation (relative):
    ξ = -0.7812: 0.12  ξ = -0.5102: 0.15  ...

  [如全部通过]
  All points converged. Use 311 for per-point diagnostics.

  Output files:
    ti_free_energy.png
    ti_free_energy.csv
    ti_convergence_report.csv
```

---

## 需新增的 K 常量（_params.py）

```python
K.EQUILIBRATION = "equilibration"
K.SEM_TARGET = "sem_target"
K.TI_ROOT_DIR = "ti_root_dir"
K.TI_DIR_PATTERN = "ti_dir_pattern"       # 注意区分已有的 K.DIR_PATTERN = "dir_pattern"（Bader 用）
K.EPSILON_TOL_EV = "epsilon_tol_ev"
```

已有可复用：`K.RESTART_PATH`、`K.LOG_PATH`、`K.COLVAR_ID`、`K.OUTDIR`。

## 代码组织

### 新增文件

- `cli/_constrained_ti.py` — `TISingleDiagCmd` 和 `TIFullAnalysisCmd` 两个命令类

### 修改文件

- `cli/__init__.py`：
  1. import `TISingleDiagCmd`, `TIFullAnalysisCmd`
  2. `build_menu_tree()` 中将 301/302 从直挂 enhanced 改为挂在 `MenuGroup("30", "Slow-Growth", output_name="slowgrowth")` 下
  3. 新增 `MenuGroup("31", "Constrained TI Analysis", output_name="constrained_ti")`，挂 311/312
  4. `enhanced.add(slowgrowth, constrained_ti)`
- `cli/_enhanced_sampling.py`：
  - **删除** `_SlowgrowthPlotCmd` 类的 `output_name = "slowgrowth"` 属性（改由父 MenuGroup 提供，防止路径重复）
  - `_discover_restart_file` / `_discover_log_file` 保留原位，311 通过 import 复用

## 文件发现的复用与差异

| 场景 | 发现逻辑 | 来源 |
|---|---|---|
| 301/302（SG） | 在 cwd 发现 .restart + .LagrangeMultLog | `_enhanced_sampling._discover_restart_file` / `_discover_log_file` |
| 311（TI 单点） | 同上（TI 目录中 `cMD-1.restart` 的 `-1` 是连字符，不被 `_\d+\.restart$` 排除） | 复用 SG 的发现函数 |
| 312（TI 多点） | 在根目录发现 `ti_target_*` / `xi_*` 子目录 | `constrained_ti.io.discover_ti_points()` |

## 未来扩展预留

| 编号 | 可能功能 |
|---|---|
| 313 | TI Compare（正/反向自由能对比） |
| 314 | TI Reweight（MBAR/BAR 重加权） |
| 315 | TI Add Points（交互式补点建议） |

---

## 审阅记录（2026-03-24）

| # | 发现 | 处置 |
|---|---|---|
| 1 | output_name 双重拼接风险：移入子组后若不删命令属性 → `slowgrowth/slowgrowth` | 已加 ⚠️ 警告，明确要求删除类属性 |
| 2 | 311 PNG 文件名：`standalone_diagnostics` 实际输出 `ti_diag_xi*.png` 而非 `ti_single_point_diag.png` | 已修正文档描述 |
| 3 | SEM_TARGET 可空 float：现有 `FloatParam` 不支持 None | 已明确用 `prompt_str` + 手动转换，覆写 `_collect_all_params()` |
| 4 | dt 一致性：312 取首点 dt 未校验 | 已在 execute 流程中加入一致性检查 + WARNING |
| 5 | TI_DIR_PATTERN 输入验证 | 已明确用手动验证限制为三个有效值 |
| 6 | output_dir 应使用 `ctx[K.OUTDIR_RESOLVED]` 而非原始值 | 已修正 311 的 API 调用示例 |
| 7 | K.TI_DIR_PATTERN vs K.DIR_PATTERN 命名易混 | 已加注释说明区别 |
