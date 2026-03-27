# MD-analysis 项目架构报告

> 生成日期：2026-03-27 | 版本：0.1.0 | Python 3.10+

---

## 1. 项目概览

**MD-analysis** 是面向 CP2K 周期性金属-水界面 MD 模拟的轻量后处理工具包。采用标准 src-layout，分发名 `md-analysis`（pip），导入名 `md_analysis`。

**入口点**：`md-analysis` 控制台命令 → `md_analysis.cli:main`（VASPKIT 风格编号菜单，无 argparse）。

**核心依赖**：numpy, matplotlib, ase, tqdm（scipy 仅标定模块 spline 拟合时 lazy import）。

---

## 2. 包结构总览

```
src/md_analysis/
├── __init__.py              # __version__, re-export: utils/water/electrochemical/potential/charge
├── config.py                # 用户持久化配置 (~/.config/md_analysis/config.json)
├── exceptions.py            # MDAnalysisError 基类
├── main.py                  # 编程入口 (run_water_analysis, run_potential_analysis, ...)
│
├── cli/                     # 交互式 CLI（菜单 1xx-9xx）
│   ├── _framework.py        # MenuNode / MenuGroup / MenuCommand / lazy_import
│   ├── _params.py           # K 键常量 + ParamCollector ABC
│   ├── _prompt.py           # 输入提示 helpers
│   ├── _water.py            # 101-105
│   ├── _potential.py        # 211-216
│   ├── _charge.py           # 221-225
│   ├── _calibration.py      # 231-233
│   ├── _enhanced_sampling.py# 301-302 (Slow-Growth)
│   ├── _constrained_ti.py   # 311-312 (Constrained TI)
│   ├── _scripts.py          # 411-412, 421-422
│   └── _settings.py         # 901-909
│
├── utils/                   # 底层解析器 + 物理常量
│   ├── config.py            # HA_TO_EV, DEFAULT_LAYER_TOL_A, ... (硬编码)
│   ├── _io_helpers.py       # _cumulative_average, _write_csv
│   ├── CubeParser.py        # CP2K cube I/O + slab_average_potential_ev
│   ├── BaderParser.py       # VASP Bader 数据加载 (POSCAR + ACF.dat + POTCAR)
│   ├── StructureParser/     # ClusterUtils, LayerParser, WaterParser
│   └── RestartParser/       # CellParser, ColvarParser
│
├── water/                   # 水分析工作流 (密度/取向/吸附层)
│   ├── Water.py             # 顶层分析函数
│   └── WaterAnalysis/       # 统计核心 (z-profile, theta PDF)
│
├── electrochemical/         # 电化学分组包
│   ├── potential/           # Hartree 电势 + 电极电势
│   │   ├── CenterPotential.py
│   │   └── PhiZProfile.py
│   ├── charge/              # Bader 电荷分析
│   │   └── Bader/           # SurfaceCharge, AtomCharges, BaderData
│   └── calibration/         # sigma -> phi 标定映射
│       ├── _data.py, _mapper.py, _plot.py
│       └── CalibrationWorkflow.py
│
├── enhanced_sampling/       # 增强采样
│   ├── slowgrowth/          # 慢增长 TI
│   │   ├── SlowGrowth.py    # 数据类 + 积分
│   │   └── SlowGrowthPlot.py# 绘图 + CSV + 统一入口
│   └── constrained_ti/      # 约束 TI 收敛诊断
│       ├── analysis/        # ACF, Block Avg, Running Avg, Geweke
│       ├── integration.py   # 梯形积分 + 误差传播
│       ├── workflow.py      # 编排器
│       ├── plot.py          # 诊断图 + 自由能图
│       └── io.py            # 目录发现
│
└── scripts/                 # 自动化脚本
    ├── BaderGen.py          # VASP Bader 工作目录生成
    ├── TIGen.py             # CP2K constrained-MD 工作目录生成
    └── utils/IndexMapper.py # XYZ <-> POSCAR 索引映射
```

---

## 3. 模块依赖图

```
                          ┌─────────┐
                          │   CLI   │  (cli/)
                          └────┬────┘
                               │ lazy_import
          ┌────────────┬───────┼───────────┬──────────┐
          v            v       v           v          v
    ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌────────┐ ┌─────────┐
    │  water   │ │potential │ │  charge  │ │calibra-│ │enhanced │
    │          │ │          │ │  (Bader) │ │ tion   │ │sampling │
    └────┬─────┘ └────┬─────┘ └────┬─────┘ └────────┘ └────┬────┘
         │            │            │                        │
         │            │            │     ┌──────────────────┤
         │            │            │     v                  v
         │            │            │  slowgrowth     constrained_ti
         │            │            │     │                  │
         └────────────┴────────────┴─────┴──────────────────┘
                               │
                          ┌────v────┐
                          │  utils  │  (config, CubeParser, BaderParser,
                          │         │   LayerParser, WaterParser, ColvarParser)
                          └─────────┘

    scripts ──→ utils.RestartParser (TIGen)
    scripts ──→ utils.config (BaderGen)
    scripts ──→ scripts.utils.IndexMapper (BaderGen)
```

**关键约束**：
- 所有箭头**单向**，无反向依赖
- `potential`、`charge`、`calibration` 三者**互不依赖**
- `enhanced_sampling` 不从 `md_analysis.__init__` re-export
- `scripts` 不从 `md_analysis.__init__` re-export
- CLI 通过 `lazy_import()` 延迟加载所有分析模块（启动零 numpy/matplotlib 开销）

---

## 4. CLI 菜单结构

```
md-analysis (root)
├── 1) Water
│   ├── 101  Water Mass Density z-Distribution
│   ├── 102  Water Orientation Weighted Density
│   ├── 103  Adsorbed Water Orientation Analysis
│   ├── 104  Adsorbed Water Theta Distribution
│   └── 105  Water Three-Panel Plot
│
├── 2) Electrochemical
│   ├── 21) Potential
│   │   ├── 211  Center Slab Potential
│   │   ├── 212  Fermi Energy
│   │   ├── 213  Electrode Potential (U vs SHE)
│   │   ├── 214  phi(z) Plane-Average Profile
│   │   ├── 215  Thickness Sensitivity
│   │   └── 216  Full Potential Analysis
│   ├── 22) Charge
│   │   ├── 221  Surface Charge (Counterion)
│   │   ├── 222  Surface Charge (Layer)
│   │   ├── 223  Surface Charge (Full)
│   │   ├── 224  Tracked Atom Charges
│   │   └── 225  Counterion Charge Tracking
│   └── 23) Calibration
│       ├── 231  Calibrate from CSV
│       ├── 232  Calibrate Manually
│       └── 233  Predict Potential
│
├── 3) Enhanced Sampling
│   ├── 30) Slow-Growth
│   │   ├── 301  Quick Plot
│   │   └── 302  Publication Plot
│   └── 31) Constrained TI Analysis
│       ├── 311  Single-Point Diagnostics
│       └── 312  Full TI Analysis
│
├── 4) Scripts / Tools
│   ├── 41) Bader Charge Preparation
│   │   ├── 411  Single Frame
│   │   └── 412  Batch from Trajectory
│   └── 42) TI Preparation
│       ├── 421  Single Work Directory
│       └── 422  Batch Generation
│
├── 9) Settings
│   ├── 901  Set VASP Script Path
│   ├── 902  Set CP2K Script Path
│   ├── 903  Set Layer Tolerance
│   ├── 904  Set Z Bin Width
│   ├── 905  Set Theta Bin Width
│   ├── 906  Set Water OH Cutoff
│   ├── 907  Reset All Defaults
│   ├── 908  Show Config
│   └── 909  Set Potential Reference
│
└── 0) Exit
```

**输出目录结构**（由 `output_name` 机制自动推导）：

```
<outdir>/
├── water/
├── electrochemical/
│   ├── potential/{center,fermi,electrode,phi_z,sensitivity}/
│   ├── charge/{counterion,layer,tracked,counterion_tracking}/
│   └── calibration/{fit,predict}/
└── enhanced_sampling/
    ├── slowgrowth/
    └── constrained_ti/
```

---

## 5. 核心数据流

### 5.1 水分析流 (101-105)

```
XYZ 轨迹 + cell_abc
      │
      v
  ase.io.iread (逐帧)
      │
      v
  LayerParser.detect_interface_layers  ──→  SurfaceDetectionResult
      │                                          │
      v                                          v
  WaterParser.detect_water_molecule_indices   界面分数坐标 (两侧)
      │                                          │
      v                                          v
  WaterAnalysis 统计核心:
    ├── z-profile 分箱 (密度 / 取向加权)       每帧独立 profile
    ├── 吸附层范围检测 (峰值 + 极小值)         → 吸附层区间
    └── theta PDF 统计                          → 角度分布
      │
      v
  系综平均 (equal-weight over frames)
      │
      v
  CSV + PNG 输出
```

**关键约定**：
- 坐标基于分数坐标第三轴（c 轴）
- 取向参考方向：晶胞 `+c`
- 界面标签：`normal_aligned`（+axis facing）/ `normal_opposed`（-axis facing）

### 5.2 电势分析流 (211-216)

```
V_HARTREE cube 文件 (CP2K)           md.out (Fermi 能级)
      │                                    │
      v                                    v
  CubeParser.read_cube()              正则解析 E_Fermi
      │                                    │
      v                                    │
  LayerParser → 界面检测 → slab 中心        │
      │                                    │
      v                                    │
  slab_average_potential_ev()              │
  (xy 平面均值 → z-profile → slab 窗口均值) │
      │                                    │
      v                                    v
  phi_center_ev ─────────────────────► 合并 cSHE 公式:
                                      U = -E_Fermi + phi_center + 15.35 - 15.81 - 0.35
                                           │
                                           v
                                      CSV + PNG (逐帧 + 累积平均)
```

**thickness sensitivity (215)**：扫描 slab 窗口厚度 [3.5, 15.0] A，输出双轴图（左：U vs SHE 均值；右：phi(z) 空间标准差）。

### 5.3 电荷分析流 (221-225)

```
bader_t*_i*/ 目录 (每帧一个)
├── POSCAR
├── ACF.dat
└── POTCAR
      │
      v
  BaderParser.load_bader_atoms()  ──→  ASE Atoms (含 bader_net_charge)
      │
      ├── method="counterion" ─────┐
      │   (非水非金属物种净电荷)     │
      │                             v
      ├── method="layer" ──────► compute_frame_surface_charge()
      │   (金属表面层净电荷)         │
      │                             v
      v                        atoms.info[surface_charge_density_uC_cm2]
  IndexMapper.remap_array()         │
  (POSCAR→XYZ 索引重映射)           v
      │                        逐帧收集 → trajectory_surface_charge (t, 2)
      v                             │
  tracked_atom_charge_analysis      v
  counterion_charge_analysis   surface_charge_analysis
      │                             │
      v                             ├── 加载 calibration.json (可选)
  CSV + PNG                         │   mapper.predict(sigma) → phi
                                    v
                               CSV + PNG (sigma ± phi 双轴)
```

**sigma -> phi 标定链**：
```
用户提供 (phi, sigma) 数据对
      │
      v
  calibrate() → 拟合 (Linear/Polynomial/Spline/DiffCap)
      │
      v
  calibration.json (持久化到 ~/.config/md_analysis/)
      │
      v
  surface_charge_analysis() 自动加载 → predict(sigma) → phi 列追加到 CSV
```

**电势参考转换**：
- SHE（默认）
- RHE：phi_RHE = phi_SHE + (RT/F) * ln(10) * pH
- PZC：phi_PZC = phi_SHE - phi_pzc

### 5.4 慢增长 TI 流 (301-302)

```
.restart + .LagrangeMultLog
      │
      v
  ColvarParser.ColvarMDInfo.from_paths()
      │
      v
  SlowgrowthFull.from_md_info()
      ├── target_growth_au: per-a.u.-time → per-step (乘 dt_au)
      └── 梯形积分: dA = -sum(lambda * d_xi)
      │
      v
  切片 [initial_step, final_step)
      │
      ├── initial > final → .reversed() (翻转 + 符号反转)
      │
      v
  slowgrowth_analysis()
      │
      ├── Quick plot: 双轴 (左: dA eV, 右: lambda a.u.) + barrier 标注
      ├── Publication plot: serif 字体 + legend 标注
      └── CSV: step, time_fs, target_au, lagrange_au, free_energy_au, free_energy_ev
```

### 5.5 约束 TI 分析流 (311-312)

#### 单点诊断 (311)

```
.restart + .LagrangeMultLog
      │
      v
  ColvarParser → lambda 时间序列
      │
      v
  analyze_standalone()
      │
      ├── Step 1: Running Average 漂移检查 (D < 3.0 * SEM)
      ├── Step 2: ACF → tau_corr, N_eff, SEM_auto
      ├── Step 3: F&P Block Average → SEM_block (pow2 平台检测)
      └── Step 4: Geweke z-test (前 10% vs 后 50%)
      │
      v
  SEM 选择: plateau → ACF fallback
      │
      v
  ConstraintPointReport
      │
      v
  2x2 诊断图 + CSV
```

#### 全 TI 分析 (312)

```
TI 根目录
├── ti_target_-0.781213/
├── ti_target_-0.239239/
├── ...
└── ti_target_1.114938/
      │
      v
  discover_ti_points(reverse=False|True)  ──→  按 xi 排序的 TIPointDefinition 列表
      │
      v
  load_ti_series()  ──→  [(xi, lambda_series, dt), ...]
      │
      v
  analyze_ti(xi_values, lambda_list, dt, epsilon_tol_ev=0.05)
      │
      ├── 每点: analyze_single_point() → ConstraintPointReport
      │
      ├── 梯形权重: compute_trapezoid_weights(xi)
      │   (支持升序/降序, reverse 时权重为负)
      │
      ├── SEM 目标: SEM_max,k = epsilon / (|w_k| * sqrt(K))
      │
      ├── 力: forces = -lambda_mean (取反)
      │
      ├── 积分: delta_A = sum(weights * forces)  [Hartree]
      │         sigma_A = sqrt(sum(w^2 * sem^2)) [Hartree]
      │
      └── 转换: delta_A_eV = delta_A * HA_TO_EV
      │
      v
  TIReport
      │
      ├── ti_convergence_report.csv  (每点一行: xi, lambda_mean, SEM, pass/fail, ...)
      ├── ti_free_energy.csv         (xi, weight, dA_dxi, sem, A_integrated_eV, sigma)
      ├── ti_free_energy.png         (双轴: 左 dA/dxi a.u., 右 A(xi) eV)
      └── ti_diag_xi*.png            (每点 2x2 诊断图)
```

### 5.6 工作目录生成流 (411-422)

#### Bader 目录生成 (411-412)

```
CP2K XYZ 轨迹 + cell_abc
      │
      v
  逐帧: ase.io.iread
      │
      v
  IndexMapper.compute_index_map()  ──→  IndexMap (XYZ↔POSCAR 双射)
      │
      v
  write_poscar_with_map()  ──→  POSCAR (元素分组 + 映射编码在注释行)
      │
      v
  生成目录: bader_t{time}_i{step}/
  ├── POSCAR     (IndexMapper 生成)
  ├── INCAR      (模板复制)
  ├── KPOINTS    (模板复制)
  ├── POTCAR     (vaspkit 103 生成)
  └── script.sh  (用户脚本复制)
```

#### TI 目录生成 (421-422)

```
SG .inp + .restart + XYZ 轨迹
      │
      ├── 时间模式: linspace(t_i, t_f, n) → CV(t) → snap 到最近帧
      └── 数值模式: 直接给 a.u. 值 → snap 到最近帧
      │
      v
  每个 TARGET:
  ├── _modify_inp_for_ti():
  │   PROJECT→cMD, STEPS→user, TARGET→snapped_au, TARGET_GROWTH→0
  │   确保 TOPOLOGY 引用 init.xyz
  └── 提取帧 → init.xyz (含 cell + PBC)
      │
      v
  生成目录: ti_target_{cv:.6f}/
  ├── cMD.inp   (修改后的输入文件)
  └── init.xyz  (提取的初始结构)
```

---

## 6. 全局约定

### 6.1 界面标签与层排序

| 标签 | 含义 | 法线方向 |
|------|------|----------|
| `normal_aligned` | 界面法线沿 +axis | 面向真空/水 |
| `normal_opposed` | 界面法线沿 -axis | 面向真空/水 |

层排序固定为：`[normal_aligned, slab_interior..., normal_opposed]`

仅支持 `"a"/"b"/"c"` 晶胞轴作为表面法线（不支持自定义向量）。

### 6.2 单位体系

| 物理量 | 内部单位 | 输出单位 |
|--------|----------|----------|
| 能量 | Hartree | eV (HA_TO_EV = 27.211386) |
| 长度 | Bohr (cube) / Angstrom (其他) | Angstrom |
| 电荷密度 | e/A^2 | uC/cm^2 (E_PER_A2_TO_UC_PER_CM2 = 1602.18) |
| CV 约束力 | Hartree/xi_unit (a.u.) | a.u. (仅积分后 dA 转 eV) |
| 时间 | fs | fs |
| 电势 | V vs SHE | V vs SHE/RHE/PZC |
| 温度 | K | K |

### 6.3 坐标约定

- 0-based 原子索引（统一）
- 分数坐标 `[0, 1)` 区间（wrap 后）
- z/c 方向统计基于分数坐标第三轴与 |c|
- 取向参考方向：晶胞 `+c`

### 6.4 输出文件约定

- CSV: 固定列名 + 顺序，首行 header
- PNG: 180 DPI（诊断图）/ 160 DPI（电荷/电势图），Agg backend
- 文件名由 `config.py` 中 `DEFAULT_*` 常量定义
- 帧目录按 `_t(\d+)` 正则数值排序（非字典序）

### 6.5 导入约定

- 包内：相对导入（`.`/`..`/`...`）
- 测试：绝对导入（`from md_analysis.xxx import ...`）
- CLI execute()：通过 `lazy_import()` 延迟加载

### 6.6 日志约定

- 库代码：`getLogger(__name__)` + `NullHandler`（PEP 282）
- CLI：`StreamHandler(stderr)` at INFO
- `verbose` 参数仅控制 tqdm 进度条

---

## 7. 配置体系

### 7.1 硬编码默认值 (`utils/config.py`)

物理常量和分析参数的单一真相源：

```python
HA_TO_EV = 27.211386245988
BOHR_TO_ANG = 0.529177249
DEFAULT_LAYER_TOL_A = 0.8      # 层聚类容差
DEFAULT_Z_BIN_WIDTH_A = 0.1    # z 方向 bin 宽度
DEFAULT_THETA_BIN_DEG = 5.0    # 角度 bin 宽度
DEFAULT_WATER_OH_CUTOFF_A = 1.25
```

### 7.2 用户持久化配置 (`~/.config/md_analysis/config.json`)

通过 CLI 设置菜单 (901-909) 管理，`_get_effective_default()` 优先读取用户配置。

### 7.3 标定文件 (`~/.config/md_analysis/calibration.json`)

sigma->phi 映射的持久化拟合结果，`surface_charge_analysis()` 自动加载。

---

## 8. 异常层次

```
MDAnalysisError                    # 基类 (exceptions.py)
├── ConvergenceError               # constrained_ti 基类
│   └── InsufficientSamplingError  # N < N_MIN
├── BaderGenError                  # scripts.BaderGen
├── TIGenError                     # scripts.TIGen
├── IndexMapError                  # scripts.utils.IndexMapper
│   └── IndexMapParseError         # POSCAR 注释行解码失败
└── CellParseError                 # CLI cell 参数解析失败
```

---

## 9. 测试策略

```
test/
├── unit/                          # 纯逻辑测试 (无 I/O)
│   ├── calibration/               # 数据加载、拟合、预测
│   ├── charge/                    # 表面电荷、原子电荷
│   ├── enhanced_sampling/
│   │   └── constrained_ti/        # ACF, Block Avg, Integration, Workflow, I/O
│   ├── potential/                 # 电势计算
│   ├── scripts/                   # BaderGen, TIGen, IndexMapper
│   └── utils/                     # LayerParser, WaterParser, ColvarParser
│
├── integration/                   # 真实数据测试
│   └── enhanced_sampling/         # data_example/ti/ 回归测试
│
└── data_example/                  # 测试数据
    ├── potential/                  # cube 文件、md.out、xyz
    ├── bader/                     # POSCAR、ACF.dat、POTCAR
    └── ti/                        # 8 个 ti_target_* 约束点目录
```

**运行**：`pip install . && pytest test/`

---

## 10. 端到端工作流示例

### 电催化模拟后处理的典型流程

```
CP2K AIMD 模拟
      │
      ├── 1) 水结构分析 (101-105)
      │   XYZ 轨迹 → 密度/取向/吸附层 → 理解界面水结构
      │
      ├── 2) 电极电势计算 (211-216)
      │   cube + md.out → U vs SHE → 确定电极电势
      │
      ├── 3) Bader 电荷 (411-412 → 外部 VASP → 221-225)
      │   XYZ → BaderGen → VASP 单点 → ACF.dat → 表面电荷密度
      │
      ├── 4) sigma-phi 标定 (231-233)
      │   多个电势下的 (U, sigma) 数据 → 拟合 → calibration.json
      │   后续 sigma 分析自动外推 phi
      │
      └── 5) 自由能计算 (301-302 → 421-422 → 311-312)
          SG 粗扫 → TIGen 生成约束 MD → 跑 CP2K → TI 分析 → delta_A
```
