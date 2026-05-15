# 菜单参考

[← 回索引](README.md)

每个 CLI 菜单代码 → 一句话用途 + 必填参数 + 输出位置。三位数字直接在 CLI 输入即可直达。

---

## 1xx — 水分析

输出根：`<outdir>/water/`

| 代码 | 用途 | 必填参数 | 输出 |
|---|---|---|---|
| 101 | 水质量密度沿 z | `xyz_path`, `cell_abc` | `water_density.csv` |
| 102 | 水取向加权密度沿 z | 同上 | `water_orientation.csv` |
| 103 | 吸附水取向沿 z | 同上 | `ad_water_orientation.csv` |
| 104 | 吸附水 θ 分布 | 同上 | `ad_water_theta.csv` |
| 105 | 105) 水三联图（包含 101-104） | 同上 | 上述 4 个 CSV + `water_three_panel.png` |

高级参数：`layer_tol`（金属层聚类容差，Å）、`outdir`、`frame_slice`。

---

## 21x — 电极电势

输出根：`<outdir>/electrochemical/potential/<sub>/`

| 代码 | 用途 | sub | 必填 |
|---|---|---|---|
| 211 | 中心 slab Hartree 势时序 | `center/` | `input_mode`, traj dir |
| 212 | Fermi 能量时序（解析 `md.out`） | `fermi/` | traj dir |
| 213 | U vs SHE（自动调 211+212 合并） | `electrode/` | traj dir |
| 214 | φ(z) 平面平均叠加（多帧） | `phi_z/` | traj dir |
| 215 | 厚度敏感性扫描（双轴：U + 空间 std） | `thickness_sensitivity/` | traj dir + 厚度上限 |
| 216 | 全套（211-215） | 各自 sub | 同上 |

`input_mode`：`continuous`（同目录 cube 文件）/ `distributed`（`potential_t*_i*` SP 子目录）。

输出文件名：`{center,fermi,electrode}_*.csv` + 同名 PNG。

---

## 22x — 表面电荷 σ

输出根：`<outdir>/electrochemical/charge/<method>/`

| 代码 | method | sub | 用途 |
|---|---|---|---|
| 221 | counterion | `counterion/` | 只算反离子 + 吸附物的 σ |
| 222 | layer | `layer/` | 界面金属层电荷之和 / 面积 |
| 223 | 全（无 method 限制） | `full/` | 跑两种方法对比 |
| 224 | counterion (单侧 + φ) | `counterion_aligned/` 或 `_opposed/` | 加 calibration → φ(t) 右轴 |
| 225 | tracked atoms | `tracked/` | 按 XYZ 索引追踪指定原子电荷 |
| 226 | counterion tracking | `counterion_tracking/` | 自动检测每帧反离子 |

必填：bader root（含 `bader_t*_i*/`）+ surface normal axis（a/b/c）。

输出：`surface_charge.csv` / `tracked_charges.csv` / `counterion_tracking.csv` + 同名 PNG。

---

## 23x — σ↔φ 标定

输出根：`<outdir>/electrochemical/calibration/{fit,predict}/`

| 代码 | 用途 | 输入 | 输出 |
|---|---|---|---|
| 231 | 从 CSV 拟合标定（linear/poly/spline/diff-cap） | CSV with σ, φ 列 | `fit/calibration.json` + `fit_curve.png` |
| 232 | 手输 (σ, φ) 数据点拟合 | 控制台输入 | 同上 |
| 233 | 已有 calibration.json 给 σ 预测 φ | calibration.json + σ 值或 CSV | `predict/predicted.csv` |

`calibration.json` 是 224 / 313 默认查找路径。

---

## 30x — Slow-Growth

输出根：`<outdir>/enhanced_sampling/slowgrowth/`

| 代码 | 用途 | 必填 |
|---|---|---|
| 301 | Quick plot（双轴：ΔA + λ） | restart + LagrangeMultLog |
| 302 | Publication plot（增强格式化） | 同上 |

输出：`slowgrowth_{quick,publication}.png` + `slowgrowth_data.csv`。

---

## 31x — Constrained TI

输出根：`<outdir>/enhanced_sampling/constrained_ti/`

| 代码 | 用途 | 必填 |
|---|---|---|
| 311 | 单点诊断（一个约束点的 4 步收敛检验） | restart + LagrangeMultLog |
| 312 | 全 TI 分析（多点 + ΔA 积分） | TI root（含多个约束点子目录） |
| 313 | 恒电势修正（Nørskov） | 312 输入 + bader 子目录 + calibration.json |

支持参数：
- `equilibration`（弃帧数）/ `epsilon_tol_ev`（精度）/ `reverse`（max ξ 为初态）/ `auto_equilibration`（二分自动弃前半）
- `point_slice`：Python 切片选子集，如 `3:8`、`::2`、`:8`

输出：
```
ti_convergence_report.csv     # 每点诊断 (sem_final, passed, failure_reasons, ...)
ti_free_energy.csv            # ξ, dA/dξ, ΔA(ξ), ±σ
ti_free_energy.png            # ΔA(ξ) 曲线 + 误差带
ti_diagnostics_<i>.png        # 每点 2x2 图（running avg / ACF / block / Geweke）
```

313 额外：`ti_const_potential_correction.{csv,png}`。

---

## 41x — Bader 工作目录生成

`<outdir>` 是用户指定，不自动加 `electrochemical/...` 前缀。

| 代码 | 用途 |
|---|---|
| 411 | 单帧 → 单个 VASP SP 工作目录 |
| 412 | 批量帧 → 多个 `bader_t<step>_i<frame>/` |

每个目录含：`POSCAR`, `INCAR`, `KPOINTS`, `run.sh`（+ `POTCAR` if requested）。

---

## 42x — TI 工作目录生成

| 代码 | 用途 |
|---|---|
| 421 | 单 TARGET 值 → 一个 `ti_target_<value>/` |
| 422 | 批量 TARGET → 多个 |

每个含：`init.xyz`（从 SG 抽帧）+ `cMD.inp`（自动改 `PROJECT cMD` / `TARGET_GROWTH 0` / `STEPS <user>`）。

模式：
- 时间区间 + n_points → linspace 在时间上 → 映射到 SG 帧的 ξ 值
- 数值列表（直接给 a.u. 值）→ snap 到最近帧

---

## 43x — SP Potential 工作目录生成

| 代码 | 用途 |
|---|---|
| 431 | 单帧 → 一个 CP2K SP 目录 |
| 432 | 批量帧 → `potential_t<step>_i<frame>/` |

输出 SP 用 211-216 的 `distributed` 模式作为输入。

---

## 44x — DeePMD SP 工作目录生成

| 代码 | 用途 |
|---|---|
| 441 | 单帧 |
| 442 | 批量帧 → `sp_t<step>_i<frame>/` |

跟 432 类似但用 **914**（DP SP inp template）的模板，专门为 DeePMD 训练集准备。

---

## 9xx — 设置

详见 [settings.md](settings.md)。简表：

| 代码 | 用途 |
|---|---|
| 900 | 显示当前配置 |
| 909 | 重置所有默认值 |
| 911 | 设置 VASP 提交脚本路径 |
| 912 | 设置 CP2K 提交脚本路径 |
| 913 | 设置 SP inp template 路径 |
| 914 | 设置 DP SP inp template 路径 |
| 921 | 金属层聚类容差（Å） |
| 922 | Z 轴 bin 宽度（Å） |
| 923 | θ bin 宽度（deg） |
| 924 | 水 O-H 截止（Å） |
| 931 | 电势输出参考（SHE / RHE / PZC） |

---

## Agent 任务列表

非交互入口 `agent.dispatch(task, params)` 的完整任务清单（**全部 contract-backed**，
入口重构期间已移除 water/potential/charge/composite legacy task）：

| Task | 关联菜单 | 契约层级 |
|---|---|---|
| `calibration_fit_csv` | 231 | **TaskContract** |
| `calibration_predict` | 233 | **TaskContract** |
| `slowgrowth_quick` | 301 | **TaskContract** |
| `ti_full_analysis` | 312 | **TaskContract**（含 parser/dir_filter） |
| `bader_gen_batch` | 412 | **TaskContract** |
| `ti_gen_batch` | 422 | **TaskContract** |
| `sp_gen_batch` | 442 | **TaskContract** |
| `config_show` | 900 | **TaskContract**（read-only） |

> 所有任务都有完整的 `FieldSpec` + `ExceptionMapping` 描述，schema 自动生成，能映射
> MCP tool schema。Phase 3 / 5B 删除的旧 task：`water_three_panel`、`potential_full`、
> `charge_surface`、`charge_tracked`、`charge_counterion`、`run_all`。对应业务现在
> 通过 `md_analysis.workflows.run_*` 直接调用即可，无需经过 agent dispatch。

---

[← 回索引](README.md)
