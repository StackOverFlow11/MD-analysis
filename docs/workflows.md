# 主线工作流

[← 回索引](README.md)

按"想做什么"组织。每条工作流给：**目标 → 输入 → 步骤 → 输出 → 常见坑**。

---

## 1. 算电极电势 U vs SHE（cSHE）

**目标**：从 CP2K AIMD 拿到 U vs SHE（绝对电极电势），含 Fermi 能量、中心势、ΔΨ 修正。

### 输入

| 类型 | 文件 |
|---|---|
| 连续轨迹模式 | `*.cube`（每帧）+ `md.out` + `md-pos-1.xyz` 在同一目录 |
| 分布式 SP 模式 | `potential_t<step>_i<frame>/sp.out + sp-v_hartree-1_0.cube` 子目录 |

CLI 第一步会问 `input_mode`：`continuous` / `distributed`。

### 步骤

走完整组合分析，直接 **216**（包含 211-215 全部）：

```
 Input: 216

 ---------- Full Potential Analysis  (includes 211-215) ----------

 Input mode  1) continuous  2) distributed [1]: 1
 Trajectory dir (contains *.cube, md.out, md-pos-1.xyz) [.]: ./md_run/
 Cell parameters [auto from md.inp]: <Enter>
 Slab thickness for averaging (A) [3.0]:
 Modify advanced parameters? (y/n) [n]: n
```

跑出来你会拿到：

```
 Computing center slab potential... done.
 Computing Fermi energy... done.
 Combining → U vs SHE...
 Computing phi(z) profile... done.
 Sweeping thickness sensitivity (1.0 → 15.0 A)... done.

 Output:
   center:                output/electrochemical/potential/center/center_slab_potential.csv
   fermi:                 output/electrochemical/potential/fermi/fermi_energy.csv
   electrode (U vs SHE):  output/electrochemical/potential/electrode/electrode_potential.csv
   phi(z) overlay:        output/electrochemical/potential/phi_z/phi_z_overlay.png
   thickness sensitivity: output/electrochemical/potential/thickness_sensitivity/...
```

### 单步骤跑法

如果你只要 Fermi 时序：用 `212`。只要 U vs SHE：用 `213`（自动调 211 + 212 合并）。

### 常见坑

- **`md.out` 找不到 Fermi**：CP2K 没开 `&PRINT &FORBIDDEN_INTERACTIONS`，或者 SCF 没收敛。看 211 的 stdout 报错。
- **cube 文件找不到**：`continuous` 模式需要 cube 文件按帧顺序排列；`distributed` 模式要 `potential_t*_i*` 命名。
- **电势参考不对**：默认 SHE。要切 RHE / vs PZC，去菜单 **931**（见 [settings.md](settings.md)）。

---

## 2. 表面电荷 σ + φ 标定

**目标**：从 VASP Bader 后处理拿到表面电荷密度 σ(t)，加上 σ→φ 标定外推电极电势。

### 输入

| 阶段 | 文件 |
|---|---|
| 帧目录 | `bader_t<step>_i<frame>/{POSCAR, ACF.dat, POTCAR}` |
| 标定数据 | CSV 表格（σ, φ 两列）或手输几个点 |

### 步骤

#### 2a. 跑 σ 时序

选方法：`counterion`（只算反离子+吸附物）/ `layer`（界面金属层电荷之和 / 面积）。

```
 Input: 221    # 或 222

 ---------- Surface Charge (Counterion) ----------

 Bader root dir (contains bader_t*_i*/) [.]: ./bader_runs/
 Surface normal axis (a/b/c) [c]: c
 Modify advanced parameters? (y/n) [n]: n
```

输出：

```
 output/electrochemical/charge/counterion/
 ├── surface_charge.csv     # 列：t_fs, sigma_aligned_uC_cm2, sigma_opposed_uC_cm2, ...
 └── surface_charge.png     # σ(t) 双线图
```

> 想同时画 φ(t) 右轴？需要先标定，菜单 224 走"σ + φ 外推"流程。

#### 2b. 标定 σ→φ

从已经跑过若干 (σ, φ) 点的数据 CSV 标定：

```
 Input: 231

 ---------- Calibrate from CSV File ----------

 CSV path: ./calibration_data.csv
 sigma column name [sigma]:
 phi column name [phi]:
 Fit type:
   1) linear   2) polynomial   3) spline   4) differential capacitance
 Choice [1]: 1
```

会输出 `output/electrochemical/calibration/fit/calibration.json`（可被 224 / 313 复用）。

#### 2c. 把 φ 外轴叠加到 σ 时序图

跑 **224**（单侧 σ + 自动加载 calibration.json）：

```
 Input: 224
 Bader root dir [.]: ./bader_runs/
 Side (aligned/opposed) [aligned]:
 Calibration JSON [./output/electrochemical/calibration/fit/calibration.json]:
```

输出 `surface_charge.csv` 多了 φ 列 + PNG 加了右轴 φ(t)。

### 常见坑

- **POTCAR 元素带 `_pv`/`_sv`**：自动剥后缀，不用手改。
- **`bader_t*_i*` 排序乱**：用 `_t(\d+)` 数值排序，不是字典序。如果你目录命名漏了 `_t`/`_i` 段会报错。
- **calibration.json 路径错**：默认在 `output/electrochemical/calibration/fit/` 下。如果你跑 231/232 用了别的 outdir，要在 224 / 313 里手工指定。

---

## 3. Slow-Growth 粗扫 → TI 精算

**目标**：CP2K 反应坐标 ξ 上的自由能 ΔA(ξ)。SG 是粗扫看趋势，TI 是定点精采样。

### 整体流程

```
 1. SG 跑一段长 traj（CP2K 端，TARGET_GROWTH > 0）
 2. 菜单 301/302 看 SG 自由能曲线，定 7-13 个 TI 点的 TARGET 值
 3. 菜单 422 在每个点生成 TI 工作目录（init.xyz + cMD.inp）
 4. 提交 CP2K 集群作业，每个目录跑约束 MD（TARGET_GROWTH = 0）
 5. 菜单 312 全 TI 收敛诊断 + ΔA 积分
 6. (可选) 菜单 313 加恒电势修正（Nørskov）
```

### 步骤

#### 3a. SG 看曲线（菜单 301）

```
 Input: 301

 ---------- Quick Plot ----------

 Restart file (e.g. slowgrowth-1_5000.restart): ./sg_run/slowgrowth-1_5000.restart
 LagrangeMultLog file [auto-discover]: <Enter>
   Auto-discovered: ./sg_run/slowgrowth-constraint_force.dat-1.LagrangeMultLog

 Plotting... → output/enhanced_sampling/slowgrowth/slowgrowth_quick.png
```

PNG 是双轴：左 ΔA(ξ) (eV)，右 λ(t) (a.u.)，下方 x 轴 = CV(a.u.)，上方 x 轴 = MD step。

#### 3b. 批量生成 TI 工作目录（菜单 422）

按时间区间 linspace 选若干点（snap 到最近 SG 帧）：

```
 Input: 422

 ---------- Batch Generate TI Work Directories ----------

 SG restart (parent of TI runs): ./sg_run/slowgrowth-1_5000.restart
 SG xyz: ./sg_run/slowgrowth-pos-1.xyz
 SG inp template (modified for TI): ./sg_run/sg.inp
 Generation mode: 1) time range  2) numeric values [1]: 1
 Initial time (fs): 1000
 Final time (fs): 5000
 Number of points: 9
 Number of MD steps per TI run [10000]:
 Output dir [./ti_runs/]:
 CP2K submit script [auto from config]:
```

会生成 9 个 `ti_runs/ti_target_<value>/`，每个含 `init.xyz` + `cMD.inp`（自动改 `PROJECT cMD` / `TARGET_GROWTH 0` / `STEPS 10000`）。直接 `qsub` 到集群。

> ⚠️ inp 模板里的 `TARGET [unit] val` 会被改成 `TARGET <au_value>`（**去掉单位标注**），统一用 a.u.。

#### 3c. TI 收敛诊断 + ΔA（菜单 312）

集群跑完 TI 之后：

```
 Input: 312

 ---------- Full TI Analysis ----------

 TI root directory [.]: ./ti_runs/
 Default equilibration frames to discard [0]: 500
 Free-energy tolerance ε (eV) [0.05]: 0.01
 Reverse integration direction (initial state = max ξ)? [n]: n
 Auto-equilibration (iteratively discard first half until converged)? [n]: n

 Found 9 constraint points:
   [0] ξ = -0.510226
   [1] ξ = -0.374922
   ...
 Select points (Python slice, e.g. 3:8, :8, 3::2, empty=all): <Enter>

 Set per-point equilibration frames? [n]: n
 Timestep: 1.000000 fs

 Point  ξ            ⟨λ⟩            SEM           N        Time range (fs)         Status
 ──────────────────────────────────────────────────────────────────────────────────
 0      -0.510226    -0.001234       0.000234     9500     1000.0 – 10500.0     PASS
 1      -0.374922    -0.000891       0.000189     9500     1000.0 – 10500.0     PASS
 ...

 Status: ALL PASS

 Output → output/enhanced_sampling/constrained_ti/
   ti_convergence_report.csv  # 每点诊断指标 + sem_final + passed
   ti_free_energy.csv         # ξ, dA/dξ, ΔA(ξ), ±σ
   ti_free_energy.png         # ΔA(ξ) 曲线 + 误差带
   ti_diagnostics_<i>.png     # 每点的 2x2 诊断（running avg / ACF / block / Geweke）
```

诊断不通过怎么办：看 `failure_reasons` 和 `time_total_fs` / `n_eff` 列；通常是采样不够、加 `equilibration` 或延长 traj。

#### 3d. 恒电势修正（可选，菜单 313）

如果 TI 是恒电荷做的（每个 ξ 点 σ 不同），想要恒电势的 ΔF(ξ)，需要 Bader 数据 + calibration：

```
 Input: 313

 ---------- Constant-Potential Correction ----------

 [先跑 Phase 1 = 走一遍 312 流程]
 ...
 [Phase 2 修正]
 Target side (aligned/opposed) [aligned]:
 Calibration JSON [.../calibration/fit/calibration.json]:

 Computing per-point sigma from bader_t*_i*/...
 Predicting Φ via mapper...
 Reference: midpoint = (σ_IS + σ_FS) / 2
 Output:
   ti_const_potential_correction.csv
   ti_const_potential_correction.png   # 原 ΔA + 修正后 ΔF 重叠
```

> 每个 `ti_target_*/` 下需要有 `bader/` 子目录（含 bader_t*_i*）。缺则 WARN + 跳过修正。

### 常见坑

- **TI 找不到点**：默认目录命名 free，但要包含 `*.restart` + `*.LagrangeMultLog`。如果你目录命名离谱（如 `pt_001/`）会被自动识别（通过 parser 内容判定）。
- **目录命名 vs 真实 ξ**：ξ 永远从 restart 的 `TARGET` 字段读，**不**从目录名 `ti_target_<x>` 解析。重命名目录不影响 ξ。
- **dt 不一致**：所有 TI 点的 timestep 必须相同，否则报错。
- **LagrangeMultLog 出现 `***`**：CP2K 数值溢出，自动转 NaN，可视化时会有 gap。诊断会记录但不崩。
- **autocorrelation N_eff < 50**：采样不够，要么延长 traj，要么调 equilibration 重新分析。

---

## 4. Bader 批量生成 VASP 工作目录

**目标**：从 CP2K MD 轨迹按帧抽出来，生成一批 VASP single-point 工作目录跑 Bader 后处理。

### 输入

| 文件 | 用途 |
|---|---|
| MD trajectory (`*.xyz`) | 帧来源 |
| cell ABC | 自动从 `.restart` 或 `md.inp` 解析 |
| INCAR / KPOINTS / 提交脚本 | 模板（在菜单 911 设置） |

### 步骤

#### 单帧（菜单 411）

```
 Input: 411

 XYZ file: ./md-pos-1.xyz
 Frame index: 0
 Cell parameters [auto]: <Enter>
 Output dir: ./bader_out/
 Generate POTCAR? (y/n) [y]: y
 VASP submit script [from config]: <Enter>
```

#### 批量（菜单 412）

```
 Input: 412

 XYZ file: ./md-pos-1.xyz
 Cell parameters [auto]: <Enter>
 Frame selection mode:
   1) Step range (start, end, stride)
   2) Time range (fs, fs, n_points)
   3) Explicit indices
 Choice [1]: 1
 Frame start: 0
 Frame end: 1000
 Frame stride: 10
 Output dir: ./bader_runs/
 Generate POTCAR? (y/n) [y]: y

 Generated 100 directories:
   bader_runs/bader_t0_i0/
   bader_runs/bader_t10_i1/
   ...
   bader_runs/bader_t990_i99/
 Each contains: POSCAR, INCAR, KPOINTS, run.sh (+ POTCAR if requested)
```

直接 `for d in bader_runs/bader_*; do (cd $d && qsub run.sh); done` 提交。

### 常见坑

- **POTCAR 不生成**：环境变量 `VASP_PP_PATH` 没设，或者元素映射缺失。
- **VASP 提交脚本路径**：第一次用要先在 **911** 设置一次，之后所有 411/412 自动用。
- **cell 解析失败**：默认从 `md.inp` 找 `ABC [angstrom]`，找不到就回退手输。

---

## 5. SP 单点为 DeePMD 训练集抽帧

**目标**：从 CP2K AIMD 轨迹按帧抽出来，生成 CP2K SP（single-point）工作目录跑高精度能量+力，作为 DeePMD 训练数据。

### 步骤

类似 Bader 批量（菜单 412），但用菜单 **442**（DeePMD SP）。区别：
- 输出 CP2K 输入文件（`sp.inp`）而非 VASP
- 用 **914**（DP SP inp template）设置的模板

```
 Input: 442

 ---------- Batch Generate SP Work Directories for DP ----------

 XYZ file: ./md-pos-1.xyz
 Cell parameters [auto]: <Enter>
 Frame selection: ...（同 412）
 Output dir: ./sp_runs/
 SP inp template [from config]: <Enter>
 CP2K submit script [from config]: <Enter>

 Generated 100 directories:
   sp_runs/sp_t0_i0/
   ...
```

### 常见坑

- **SP inp template 没设**：先去 **914** 配。
- **重名碰撞**：如果 `sp_runs/sp_t<step>_i*/` 已存在，会跳过且 warn（不覆盖）。

---

## 非交互场景出口

CLI 是给人看的。如果你要在脚本 / Jupyter / 集群作业里跑：

### 程序化入口（`main.py`）

```python
from md_analysis.main import (
    run_water_analysis, run_potential_analysis, run_charge_analysis,
    run_tracked_charge_analysis, run_counterion_charge_analysis, run_all,
)

# 跟 CLI 菜单一一对应；output_dir 是最终写入目录（不再前置 water/ 等）
run_water_analysis(
    xyz_path="md-pos-1.xyz",
    cell_abc=(10.22, 10.22, 26.42),
    output_dir="output/water/",
)
```

### Agent 入口（`agent.dispatch`）

JSON 序列化、JSON Schema 自描述，方便给 LLM / MCP / 批处理用：

```python
from md_analysis.agent import dispatch, list_tasks, get_task_schema

# 列出所有 task
for t in list_tasks():
    print(t["name"], "—", t["description"])

# 拿单个 task 的 JSON Schema
schema = get_task_schema("ti_full_analysis")

# 执行
result = dispatch("ti_full_analysis", {
    "root_dir": "./ti_runs/",
    "output_dir": "./output/",
    "epsilon_tol_ev": 0.01,
    "equilibration": 500,
})
print(result.success, result.summary["delta_A_eV"])
```

任务列表（14 个）：见 [menu_reference.md 末尾](menu_reference.md#agent-任务列表)。

---

[← 回索引](README.md)
