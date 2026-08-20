# 工作流 03 — Slow-Growth 粗扫 → TI 精算

[← 回 workflows 索引](../workflows.md) | [← 回主索引](../README.md)

**目标**：CP2K 反应坐标 ξ 上的自由能 ΔA(ξ)。SG 粗扫看趋势，TI 在选定点精采样。

---

## 1. 物理量与公式

| 量 | 含义 | 单位 |
|---|---|---|
| ξ(t) | CV（集体变量）随时间 | a.u. |
| λ(t) | SHAKE 约束力（Lagrange 乘子） | a.u. |
| dξ/dt | TARGET_GROWTH | a.u./a.u.time |
| dt | timestep | fs / a.u.time |

**SG 积分**（midpoint, CP2K 负号约定）：

```
ΔA_k = -∑_{i=0}^{k-1} (λ_i + λ_{i+1})/2 × Δξ_per_step
其中 Δξ_per_step = TARGET_GROWTH × (dt_fs / 0.02418884326585)
```

**TI 积分**（梯形法 + 取负号）：

```
dA/dξ_k = -⟨λ_k⟩
ΔA = ∑_k w_k × dA/dξ_k        # w_k 是 trapezoid 权重
σ_A² = ∑_k (w_k × SEM_k)²     # 误差传播
```

---

## 2. 整体流程

```
 1. SG 跑一段长 traj（CP2K 端，TARGET_GROWTH > 0）
                     │
                     ▼
 2. 菜单 301/302 看 SG 自由能曲线，定 7-13 个 TI 点的 TARGET 值
                     │
                     ▼
 3. 菜单 422 在每个点生成 TI 工作目录（init.xyz + cMD.inp）
                     │
                     ▼
 4. 提交 CP2K 集群作业（每个目录独立约束 MD，TARGET_GROWTH = 0）
                     │
                     ▼
 5. 菜单 312 全 TI 收敛诊断 + ΔA 积分
                     │
                     ▼
 6. (可选) 菜单 313 加恒电势修正（Nørskov）
```

---

## 3. SG 输入数据结构

```
sg_run/
├── sg.inp                                    # CP2K 输入（含 &CONSTRAINT 块）
├── slowgrowth-1.restart                      # 最新 restart（用于解析元数据）
├── slowgrowth-1_<step>.restart               # 中间 restart 备份
├── slowgrowth-pos-1.xyz                      # 轨迹（生成 TI init.xyz 用）
├── slowgrowth-1.ener                         # 能量时序（不直接读，参考用）
├── slowgrowth-constraint_force.dat-1.LagrangeMultLog   # ★ 约束力时序
└── *.bak                                     # CP2K backup files（自动跳过）
```

文件命名约定来自 CP2K 默认（`<PROJECT>-<run_id>_<step>.restart`、`<PROJECT>-<filename>-1.LagrangeMultLog`）。

### 3.1 SG inp 关键 section

```
&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 8000
    TIMESTEP 0.5
    ...
  &END MD

  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR .TRUE.
      TARGET [angstrom]  5.64                ! 起始 ξ
      TARGET_GROWTH [fs^-1*angstrom] -0.5E-3 ! 增长率（带方向，可正可负）
    &END COLLECTIVE

    &LAGRANGE_MULTIPLIERS  SILENT             ! ★ 必须开！
      COMMON_ITERATION_LEVELS 1
      FILENAME constraint_force.dat           ! 输出名后缀
    &END LAGRANGE_MULTIPLIERS
  &END CONSTRAINT

  &PRINT
    &TRAJECTORY SILENT
      &EACH MD 5 &END EACH                    ! 帧输出频率（每 5 步一帧）
    &END TRAJECTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  &SUBSYS
    &COLVAR
      ! 这里定义 CV 类型，比如 DISTANCE / ANGLE / COORDINATION
      &DISTANCE
        ATOMS 12 23
      &END DISTANCE
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL
```

> ⚠️ `&LAGRANGE_MULTIPLIERS` 是**必须**开的，没它就没有 `*.LagrangeMultLog` 文件，分析没法做。
> ⚠️ TARGET 和 TARGET_GROWTH 可以带单位标注（CP2K 解析 OK），但 restart 写出来时**统一转成 a.u.**，所以工具内部读到的永远是 a.u.

---

## 4. SG 内部数据流

```
sg_run/slowgrowth-1.restart
   │
   ▼ parse_colvar_restart() (regex on &MD, &CONSTRAINT, &CELL)
   │
ColvarRestart {
   project_name, step_start, time_start_fs, timestep_fs,
   colvars: ColvarInfo {
     constraints: (
       ConstraintInfo {colvar_id, target_au, target_growth_au, ...},
       ...
     )
   },
   cell_abc_ang, fixed_atom_indices, ...
}

sg_run/*.LagrangeMultLog
   │
   ▼ parse_lagrange_mult_log() (line-by-line)
   │
LagrangeMultLog {
   shake: np.ndarray (N_steps,) 或 (N_steps, K)
   rattle: 同形状
   n_steps, n_constraints
}

   合并 → ColvarMDInfo (.steps, .times_fs, .target_series_au, ...)
   │
   ▼ SlowgrowthFull.from_md_info(md_info, colvar_id=None)
   │   计算每步 Δξ_per_step = target_growth_au × dt_au
   │   midpoint integration: ΔA[k] = -cumsum((λ[i]+λ[i+1])/2 × Δξ)
   │
Slowgrowth dataclass {
   steps, times_fs, target_au (ξ 时序),
   lagrange_shake (λ 时序),
   free_energy_au (ΔA(t) 累积),
   timestep_fs, target_growth_au (per-step)
}

   ▼ plot_slowgrowth_quick() / publication()
   │
output/enhanced_sampling/slowgrowth/
├── slowgrowth_quick.png        # 双轴：左 ΔA(eV), 右 λ(a.u.)；底 x = CV(a.u.), 顶 x = step
├── slowgrowth_publication.png
└── slowgrowth_data.csv         # step, time_fs, target_au, lambda_shake, free_energy_eV
```

---

## 5. SG 步骤（菜单 301）

```
$ md-analysis
 Input: 301

 ---------- Quick Plot ----------

 Restart file: ./sg_run/slowgrowth-1.restart
 LagrangeMultLog [auto-discover]: <Enter>
   Auto-discovered: ./sg_run/slowgrowth-constraint_force.dat-1.LagrangeMultLog

 Loading... 8000 steps, dt=0.5 fs, ξ ∈ [5.64, 1.64] (decreasing)
 Plotting...

 Output: ./output/enhanced_sampling/slowgrowth/slowgrowth_quick.png
```

---

## 6. 从 SG 生成 TI 工作目录（菜单 422）

### 6.1 选 TI 点的两种模式

| 模式 | 用法 | 特点 |
|---|---|---|
| time range | 给 (t_start_fs, t_end_fs, n_points) | linspace 在时间上 → 映射 SG 帧 → 取该帧 ξ |
| numeric | 给一组 a.u. 值 | 直接 snap 到最近 SG 帧 |

为什么不直接给 ξ 值就好？因为 TI 需要 SG 轨迹的真实快照（init.xyz）作为初态——必须跟某个真实帧对齐。snap 到最近帧是为了保证 init.xyz + TARGET 自洽。

### 6.2 输出目录结构

```
ti_runs/
├── ti_target_-0.510226/
│   ├── init.xyz                         # 从 SG 轨迹抽该帧
│   └── cMD.inp                          # 从 sg.inp 改造
├── ti_target_-0.374922/
│   ├── init.xyz
│   └── cMD.inp
└── ti_target_<...>/
```

⚠️ 目录名里的数值是 **a.u.**（CP2K 默认），不是用户输入的 Å / deg。

### 6.3 cMD.inp 跟 sg.inp 的差异

工具自动改的 4 个地方：

| 关键字 | sg.inp（原） | cMD.inp（改后） |
|---|---|---|
| `PROJECT` | `slowgrowth` | `cMD` |
| `TARGET [unit] val` | `[angstrom] 5.64` | `<au_value>`（**去掉单位标注**，统一 a.u.）|
| `TARGET_GROWTH [unit] val` | `[fs^-1*angstrom] -0.5E-3` | `0` |
| `STEPS` | 8000 | 用户指定（默认 10000） |

**为什么去单位标注？** restart 里所有 ξ / λ 内部都是 a.u.，统一用避免量纲坑（特别是配位数 CV 这种"数量"很难写单位的）。

### 6.4 步骤

```
$ md-analysis
 Input: 422

 ---------- Batch Generate TI Work Directories ----------

 SG restart: ./sg_run/slowgrowth-1.restart
 SG xyz: ./sg_run/slowgrowth-pos-1.xyz
 SG inp template: ./sg_run/sg.inp
 Generation mode (1=time, 2=numeric) [1]: 1
 Initial time (fs): 1000
 Final time (fs): 5000
 Number of points: 9
 Number of MD steps per TI run [10000]:
 Output dir [./ti_runs/]:
 CP2K submit script [auto from config]:

 Found 9 target points (snapped to nearest SG frames):
   ti_target_-0.510226  (frame 100, t=1000.0 fs)
   ti_target_-0.374922  (frame 150, t=1500.0 fs)
   ...

 Generated 9 directories under ./ti_runs/.
```

提交：`for d in ti_runs/ti_target_*; do (cd $d && qsub run.sh); done`

---

## 7. TI 输入数据结构（提交完成后）

```
ti_runs/
├── ti_target_-0.510226/
│   ├── init.xyz
│   ├── cMD.inp                                      # 修改后的 inp
│   ├── cMD-1.restart                                # ★ 元数据来源
│   ├── cMD-1_<step>.restart                         # 中间 restart
│   ├── cMD-pos-1.xyz                                # 该点轨迹
│   └── cMD-constraint_force.dat-1.LagrangeMultLog   # ★ λ 时序
├── ti_target_-0.374922/
└── ...
```

⚠️ 目录名 `ti_target_<value>` 只是命名习惯，**不影响识别**。重构后（commit `4541369`）TI 流程是 parser-driven + 内容判定：
- 默认 `parser="auto"` 嗅探 → 找到 CP2K 的 `*.restart + *.LagrangeMultLog`
- 默认 `dir_filter=None` → 内容过滤（任何含这两个文件的子目录都算 TI 点）
- ξ 永远从 restart 的 `TARGET` 字段读，不从目录名解析

---

## 8. TI 内部数据流（菜单 312）

```
ti_runs/
└── ti_target_*/
    ├── *.restart  ─► CP2KParser.parse_metadata ─► ColvarRestart
    └── *.LagrangeMultLog ─► CP2KParser.parse_lambda_series ─► LagrangeMultLog

   discover_ti_points(root, parser="auto", dir_filter=None)
   │   每个子目录预读 metadata 拿 ξ
   │   按 ξ 升序排序（ties: 目录名）
   ▼
TIPointDefinition[K] {
   directory, parser, metadata: ColvarRestart
   xi (property → metadata.colvars.primary.target_au)
}

   load_ti_series(point_defs)  → [(xi, λ_series, dt_fs), ...]
   │
   ▼ analyze_ti(xi_values, lambda_list, dt, equilibration, ...)
   │
   for each point:
     analyze_single_point(λ', xi, dt) → ConstraintPointReport
       ┌─ 1. ACF (Sokal 1997 self-consistent cutoff) → τ_corr, N_eff, SEM_auto
       ├─ 2. Block average (Flyvbjerg-Petersen 1989, pow2 + δSEM plateau) → SEM_block
       ├─ 3. Running average drift check → D < c × SEM
       └─ 4. Geweke stationarity test → |z| < 1.96
     ⟨λ⟩, σ_λ, sem_final, passed, failure_reasons

   collective:
     trapezoidal weights w_k → forces = -⟨λ⟩ → ΔA = ∑w_k × force_k
     SEM 传播 → σ_A
   ▼
TIReport {
   point_reports: tuple[ConstraintPointReport, ...]
   xi_values, weights, forces, force_errors
   delta_A, sigma_A, all_passed, failing_indices
}

   ▼ write_convergence_csv / write_free_energy_csv / plot_*
   ▼
output/enhanced_sampling/constrained_ti/
├── ti_convergence_report.csv   # 每点 24 列诊断（sem_final, sem_inflated, passed, failure_reasons, ...）
├── ti_free_energy.csv          # ξ, dA_dxi_au, dA_dxi_eV, delta_A_eV (cumulative), sigma_eV
├── ti_free_energy.png          # ΔA(ξ) 曲线 + 误差带
└── ti_diag_xi<...>.png × K     # 每点 2x2 诊断图（running avg / ACF / block / Geweke）
```

### 4 步诊断阈值表

| 诊断 | 阈值 | 意义 |
|---|---|---|
| ACF | `N_eff ≥ 10`（地板） | 低于地板则 τ_corr 估计失效，直接 fail |
| Block average | F&P plateau detected | δSEM 趋于平缓后取值 |
| Running average | drift `D < c × SEM` | 累积均值不漂移 |
| Geweke | `|z| < 1.96` | 前后段统计无显著差异 |
| SEM 报告值 | `SEM_report ≤ SEM_max` | 平台 SEM × 卡方膨胀因子 √(ν/χ²₀.₀₅(ν))，ν = 平台块数 − 1（单侧 95% 上界） |

`sem_final` 来源优先级：F&P plateau → ACF fallback。

---

## 9. TI 步骤（菜单 312）

```
$ md-analysis
 Input: 312

 ---------- Full TI Analysis ----------

 TI root directory [.]: ./ti_runs/
 Default equilibration frames to discard [0]: 500
 Free-energy tolerance ε (eV) [0.05]: 0.01
 Reverse integration direction (max ξ = initial)? [n]: n
 Auto-equilibration (iteratively discard front half)? [n]: y

 Found 9 constraint points:
   [0] ξ = -0.5102
   [1] ξ = -0.3749
   ...
 Select points (Python slice, e.g. 3:8): <Enter>
 Set per-point equilibration? [n]: n
 Timestep: 0.5 fs

 Point  ξ           ⟨λ⟩         SEM         N      Time range (fs)    Status
 ─────────────────────────────────────────────────────────────────────────
 0     -0.5102     +0.005244    0.000434   8500   1000.0 – 9500.0   PASS
 1     -0.3749     +0.003251    0.000570   7800   ...               PASS
 ...

 Status: ALL PASS
 ΔA = -0.2136 ± 0.0050 eV

 Output → output/enhanced_sampling/constrained_ti/
```

---

## 10. 恒电势修正（菜单 313，可选）

如果 TI 是恒电荷做的（每个 ξ 点 σ 不同），想要恒电势 ΔF(ξ)，需要 **每个 TI 点下面有 bader/ 子目录** 作为 σ 数据源 + 已经标定好的 calibration.json：

```
ti_runs/
├── ti_target_-0.5102/
│   ├── *.restart
│   ├── *.LagrangeMultLog
│   └── bader/                          # ★ 恒电势修正必需
│       ├── bader_t<step>_i<frame>/{POSCAR, ACF.dat, POTCAR}
│       └── ...
└── ...
```

修正公式（Nørskov，midpoint reference）：

```
ΔF_Φ(ξ) = ΔF_q(ξ) + [σ(ξ) − σ_ref] × [Φ(ξ) − Φ_ref] × A / 2

σ_ref = (σ_IS + σ_FS) / 2     # 初末态 σ 中点
Φ_ref = (Φ_IS + Φ_FS) / 2     # 初末态 Φ 中点（外推自 calibration.json）
A     = 电极表面积（POSCAR cell 沿法向的两向量叉积）
```

---

## 11. 常见坑

- **目录命名不影响**（重构后）：但 ξ 排序是看 restart 内容，确认 restart 没问题。
- **dt 不一致**：所有 TI 点 timestep 必须相同，否则报错。
- **`LagrangeMultLog` 出现 `***`**：CP2K 数值溢出，自动转 NaN，可视化有 gap。多了说明数值不稳定。
- **N_eff < 10（地板）/ Geweke 不平稳**：采样不够 / 前期没平衡。先试 `auto_equilibration=True`，自动二分砍前半。N_eff ≥ 10 但 SEM 不达标时看 `SEM_report`（卡方 95% 上界）与 `sem_inflation_factor`——块数少会放大报告 SEM，延长 traj 即可缩小。
- **ΔA 数量级离谱**：检查 SG 的 TARGET_GROWTH 单位（per-fs vs per-au-time）。本工具内部只信 restart（已转 a.u.），但若你手改了 inp 的 TARGET_GROWTH 数值忘了同步，会出错。
- **`ti_target_<value>` 跟 restart 里 TARGET 对不上**：重构后只信 restart；目录名只是标识符。但如果对不上你应该自查（可能是 422 生成时 SG 帧 snap 到错的步）。
- **`auto_equilibration=True` 报 `InsufficientSamplingError`**：二分到 < 100 帧仍不通过，意味着采样彻底不够，必须延长 traj。

---

[← 02_surface_charge.md](02_surface_charge.md) | [→ 04_bader_batch.md](04_bader_batch.md)
