# Data Contract（跨模块核心契约）

> 仅记录跨模块共享、且会影响结果可比性的核心数据契约。
>
> 接口暴露清单与实现细则请查看：
> - `context4agent/architecture/modules/src/interface_exposure.md`
> - `context4agent/architecture/modules/src/utils/interface_exposure.md`
> - `context4agent/architecture/modules/src/charge/interface_exposure.md`
> - `context4agent/architecture/modules/src/utils/implementation_guidelines.md`
>
> 记录落位硬约束：本文件仅承载"全局契约"；非全局目录级记录不得写入本文件，
> 必须写入对应目录的 `interface_exposure.md` / `implementation_guidelines.md`。

## 适用范围

- 处理对象：单帧 `ase.Atoms`
- 原子索引：统一使用 **0-based**
- 分数坐标：默认按 wrap 后区间处理（`[0, 1)`）
- 界面标签：`"normal_aligned"` (outward normal = +axis) / `"normal_opposed"` (outward normal = −axis)
- 水分析起始界面：`StartInterface = Literal["normal_aligned", "normal_opposed"]`（默认 `"normal_aligned"`）

## 核心数据载体（输出形状契约）

- 水分子标记：`(n_water, 3)`，每行 `[O_idx, H1_idx, H2_idx]`
- 水氧索引：`(n_water, 1)`（输入可接受 `(n,)` 或 `(n, 1)`）
- 质量密度 z 分布：`(nbins, 1)`，单位 `g/cm^3`
- 取向加权 z 分布：`(nbins, 1)`，单位 `g/cm^3`
- c 窗口角度 PDF：`(180 / ndeg,)`，单位 `degree^-1`

## 坐标与方向契约

- z/c 方向统计基于分数坐标第三轴与 `|c|`
- 取向参考方向统一为晶胞 `+c`（`c_unit = c / |c|`）
- c 轴窗口采用周期区间（左闭右开）：
  - `start < end`：`[start, end)`
  - `start > end`：`[start, 1) U [0, end)`
  - `start == end (mod 1)`：全区间

## 统计量定义（规范口径）

### 1) 水质量密度 z 分布

分箱与体积：

$$
nbins = \left\lceil \frac{L_z}{dz} \right\rceil,\quad
V_{\mathrm{bin}} = A_{xy}\,\Delta z,\quad
A_{xy} = \lVert \mathbf{a}\times\mathbf{b}\rVert
$$

密度：

$$
\rho_{\mathrm{bin}} = \frac{N_{O,\mathrm{bin}}\cdot M_{\mathrm{H_2O}}/N_A}{V_{\mathrm{bin}}}
$$

说明：一个 O 计作一个 H2O 分子。

### 2) 取向加权 z 分布

角度定义：$\theta$ 为 H-O-H 角平分线与 `+c` 方向夹角。

$$
P_{\mathrm{orient},\mathrm{bin}} = \frac{\sum_{i\in \mathrm{bin}}\cos\theta_i \cdot m_{\mathrm{H_2O}}}{V_{\mathrm{bin,cm^3}}}
$$

说明：$m_{\mathrm{H_2O}} = M_{\mathrm{H_2O}} / N_A$，$V_{\mathrm{bin,cm^3}} = V_{\mathrm{bin,Å^3}} \times 10^{-24}$。
输出单位与质量密度一致，均为 `g/cm^3`。

### 3) 指定 c 窗口的角度 PDF

- 角域固定 `[0, 180]` 度
- 要求 `180 / ndeg` 为整数
- 窗口内无样本时，返回全零数组

## 默认参数来源

- 默认值统一来自：`src/md_analysis/utils/config.py`
- 常用项包括：`DEFAULT_Z_BIN_WIDTH_A`、`DEFAULT_THETA_BIN_DEG`、`DEFAULT_WATER_OH_CUTOFF_A`
- 当前默认值（与代码同步）：
  - `DEFAULT_Z_BIN_WIDTH_A = 0.1` Angstrom
  - `DEFAULT_THETA_BIN_DEG = 5.0` degree
  - `DEFAULT_WATER_OH_CUTOFF_A = 1.25` Angstrom

## Water 层输出契约

### `water_mass_density_z_distribution_analysis(...)` CSV

- 输出文件默认参数来自：`src/md_analysis/water/config.py`
- 界面定义：使用直接面向非金属环境的 interface 对（每侧一层）
- CSV 列定义：
  - `path_fraction_center`：界面 -> 两界面中点路径的归一化坐标（`[0, 1]`）
  - `distance_A`：对应平均物理距离（Angstrom）
  - `rho_ensemble_avg_g_cm3`：系综平均水质量密度（`g/cm^3`）

### `water_orientation_weighted_density_z_distribution_analysis(...)` CSV

- 输出文件默认参数来自：`src/md_analysis/water/config.py`
- 界面定义：使用直接面向非金属环境的 interface 对（每侧一层）
- CSV 列定义：
  - `path_fraction_center`：界面 -> 两界面中点路径的归一化坐标（`[0, 1]`）
  - `distance_A`：对应平均物理距离（Angstrom）
  - `orientation_ensemble_avg_g_cm3`：系综平均取向加权密度（`g/cm^3`）

### `ad_water_orientation_analysis(...)` 输出契约

- 输出文件默认参数来自：`src/md_analysis/water/config.py`
- 吸附层主峰定义：密度分布最高 bin 的距离位置
- 吸附层区间：
  - 下界：主峰前最后一个近零点
  - 上界：主峰后第一个局部极小值（平滑后）
- 结果文件：
  - profile CSV：`distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_g_cm3,is_adsorbed_layer_bin`
  - range TXT：`adsorbed_layer_start_A, adsorbed_layer_end_A, main_peak_distance_A, ...`

### `compute_adsorbed_water_theta_distribution(...)` 输出契约

- 统计对象：吸附层区间内水分子的取向角 `theta`
- 角域：`0-180` degree
- 输出 CSV：`theta_degree,pdf_degree_inv`

### `plot_water_three_panel_analysis(...)` 输出契约

- 入口定位：集成三联图的推荐入口（同时落盘各中间结果，便于复现）
- 输出目录：`output_dir`（默认 `Path.cwd()`）
- 输出文件（默认命名来自：`src/md_analysis/water/config.py`；PNG 名可覆盖）：
  - PNG：`DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME`（参数 `output_png_name` 可覆盖）
  - density CSV：`DEFAULT_WATER_MASS_DENSITY_CSV_NAME`
  - orientation CSV：`DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME`
  - adsorbed profile CSV：`DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME`
  - adsorbed range TXT：`DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME`
  - theta distribution CSV：`DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME`
- CSV 列（header 与顺序固定）：
  - density CSV：`path_fraction_center,distance_A,rho_ensemble_avg_g_cm3`
  - orientation CSV：`path_fraction_center,distance_A,orientation_ensemble_avg_g_cm3`
  - adsorbed profile CSV：`distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_g_cm3,is_adsorbed_layer_bin`
  - theta distribution CSV：`theta_degree,pdf_degree_inv`
- range TXT 键（key 固定，value 为数值）：
  - `adsorbed_layer_start_A=...`
  - `adsorbed_layer_end_A=...`
  - `main_peak_distance_A=...`
  - `near_zero_ratio=...`
  - `smoothing_window_bins=...`
- 轨迹读取次数：2（第一次密度+取向；第二次吸附层 `theta` 分布）

### 系综平均口径（A）

- 每帧独立计算从选定界面到中点的单帧 profile
- 若各帧 `nbins` 不一致，先重采样到统一归一化坐标网格
- 最后按帧等权平均（equal-weight over frames）

## Charge 层输出契约

### `compute_frame_surface_charge(atoms, *, method="counterion", ...)` 输出

- `method` 参数选择计算逻辑：
  - `"counterion"`：排除水分子和金属原子，仅非水非金属物种（反离子/溶质）的净电荷贡献 σ
  - `"layer"`：对界面层金属原子的净电荷求和除以面积。`n_surface_layers` 参数控制每侧取几层金属层（默认 1 = 仅最外层；设为 N 则取最外 N 层）
- 结果存入 `atoms.info`（原地修改）
- `atoms.info["surface_charge_density_e_A2"]`：`[σ_aligned, σ_opposed]`，单位 e/Å²
- `atoms.info["surface_charge_density_uC_cm2"]`：`[σ_aligned, σ_opposed]`，单位 μC/cm²
- `atoms.info["n_charged_atoms_per_surface"]`：`[n_aligned, n_opposed]`
- `atoms.info["charge_per_surface_e"]`：`[Σq_aligned, Σq_opposed]`，单位 e
- 排列顺序由稳定的 `interface_label` 决定（aligned=index 0, opposed=index 1），不受 PBC 平移影响
- `normal` 参数控制面积计算：`_AREA_VECTORS = {"a": (1,2), "b": (0,2), "c": (0,1)}`

### `frame_indexed_atom_charges(atoms, atom_indices)` 输出

- 输入 `atom_indices`：`(N,)` 0-based 整型数组
- 返回 `np.ndarray`：`(N, 2)`
  - `[:, 0]`：回显的原子索引
  - `[:, 1]`：对应的 Bader 净电荷（单位 e）

### `trajectory_indexed_atom_charges(root_dir, atom_index_matrix, ...)` 输出

- 输入 `atom_index_matrix`：`(t, N)` 0-based 整型数组
- 返回 `np.ndarray`：`(t, N, 2)`
  - `[:, :, 0]`：回显的原子索引
  - `[:, :, 1]`：对应的 Bader 净电荷（单位 e）
- 内部逐帧调用 `frame_indexed_atom_charges`

### `trajectory_surface_charge(root_dir, *, method="counterion", ...)` 输出

- 返回 `np.ndarray`：`(t, 2)`
  - `[:, 0]`：σ_aligned（μC/cm²）
  - `[:, 1]`：σ_opposed（μC/cm²）
- 逐帧调用 `load_bader_atoms` + `compute_frame_surface_charge(..., method=method)`，收集 `surface_charge_density_uC_cm2`

### `surface_charge_analysis(root_dir, *, method="counterion", ...)` 输出

- CSV 基础列：`step,sigma_aligned_uC_cm2,sigma_opposed_uC_cm2,sigma_aligned_cumavg_uC_cm2,sigma_opposed_cumavg_uC_cm2`
- CSV 可选列（当 `~/.config/md_analysis/calibration.json` 存在时自动追加）：`phi_aligned_V_vs_{REF},phi_opposed_V_vs_{REF},phi_aligned_cumavg_V_vs_{REF},phi_opposed_cumavg_V_vs_{REF}`（REF 由 `potential_reference` 参数决定，默认 SHE；可选 RHE/PZC）
- PNG：9×4.8 inch, 160 DPI，左轴 σ（aligned 蓝/opposed 橙，inst. + cum. avg）；右轴 φ（aligned 绿/opposed 红，仅有标定时显示，轴标签 `φ (V vs {REF})`）+ fit RMSE 标注
- 默认文件名：`surface_charge.csv`、`surface_charge.png`
- 返回 CSV 路径

## Potential 帧数据抽象

### `PotentialFrame` dataclass (`_frame_source.py`)

- `step: int` — MD 步数
- `time_fs: float | None` — 模拟时间（fs，分布式模式从目录名提取；连续模式为 None）
- `cube_path: Path` — cube 文件路径
- `header: CubeHeader` — 已解析的 cube 文件头
- `values: np.ndarray` — flat cube 数据，length `nx * ny * nz`
- `fermi_raw: float | None` — Fermi 能（Hartree），None 表示无数据
- `atoms: Atoms | None` — 用于 interface detection（含 cell）

### 输入模式发现契约

- **Continuous**：所有 cube 文件在单个目录，glob 匹配；Fermi 能从单个 md.out
- **Distributed**：多个 `potential_t{time}_i{step}/` 子目录，每个含一个 cube + sp.out
- 两种模式均支持 `frame_start/frame_end/frame_step` 切片
- 分布式模式下未完成计算的子目录（无 cube 文件）自动跳过

## Potential 层输出契约

### `center_slab_potential_analysis(...)` 输出

- CSV：`step,phi_center_ev,phi_center_cumavg_ev`
- 附加 CSV（`slab_center_and_interfaces.csv`）：`step,center_source,z_center_ang,z_iface_lower_ang,z_iface_upper_ang,z_iface_mid_ang,water_gap_ang,n_metal_layers,Lz_ang`
- PNG：逐帧 + 累积平均曲线

### `fermi_energy_analysis(...)` 输出

- CSV：`step,time_fs,fermi_raw,fermi_ev,fermi_cumavg_ev`
- PNG：逐帧 + 累积平均曲线

### `electrode_potential_analysis(...)` 输出

- 内部调用 `center_slab_potential_analysis` + `fermi_energy_analysis`
- 合并公式：`U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP`
- CSV：`step,U_vs_SHE_V,U_cumavg_V`
- PNG：逐帧 + 累积平均曲线

### `phi_z_planeavg_analysis(...)` 输出

- CSV：`z_ang,phi_mean_ev,phi_std_ev,phi_min_ev,phi_max_ev`
- PNG：所有帧的 φ(z) overlay + 系综平均 ± std 带

### `thickness_sensitivity_analysis(...)` 输出

- 扫描 thickness 范围：`[thickness_start, thickness_end]`（默认 3.5–15.0 Å，步长 0.5 Å）
- 每个 thickness 计算所有帧的 `U vs SHE`，取系综平均
- CSV：`thickness_ang,mean_U_vs_SHE_V,mean_phi_z_spatial_std_eV,n_frames`
- PNG 双轴图：
  - 左轴：`mean U vs SHE (V)` — 系综平均电极电势
  - 右轴：`spatial std of φ(z) in slab (eV)` — slab 区间内 Hartree 势沿 z 的空间标准差（帧系综平均）
- 连续模式需要 `md.out`（Fermi 能级）；分布式模式从各子目录 `sp.out` 获取；均无时跳过

## Enhanced Sampling 层输出契约

### Constrained TI 输出契约

#### `write_convergence_csv(ti_report)` — 多点收敛报告

- 默认文件名：`ti_convergence_report.csv`
- 每约束点一行，按 ξ 排列
- λ 相关量单位为 a.u.
- CSV 列：`xi, lambda_mean, sigma_lambda, tau_corr, n_eff, sem_auto, sem_block, delta_sem_block, plateau_B, plateau_reached, sem_final, sem_final_method, sem_max, geweke_z, geweke_reliable, drift_D, passed, failure_reasons`

#### `write_free_energy_csv(ti_report)` — 自由能曲线

- 默认文件名：`ti_free_energy.csv`
- dA/dξ 单位 a.u.；积分后 A 单位 eV
- CSV 列：`xi, weight, dA_dxi, sem, A_integrated_eV, sigma_A_cumulative_eV`

#### `write_single_point_csv(report)` — 单点诊断报告

- 默认文件名：`ti_single_point.csv`
- 单行 CSV，列同 `write_convergence_csv`

#### 诊断图（PNG）

- **单点诊断图**（`ti_diag_xi{value}.png`）：2×2 子图
  - 左上：累积均值 λ̄(n) + ±SEM 带
  - 右上：ACF C(j) + 截断点标记
  - 左下：SEM(B) vs B + 平台区域高亮 + δSEM 误差棒
  - 右下：文字摘要表（τ_corr, N_eff, SEM, Geweke z, pass/fail）
- DPI: 180

#### 自由能曲线图（PNG）

- **自由能图**（`ti_free_energy.png`）：双轴图
  - 左轴：dA/dξ with error bars (a.u.)，点按 PASS/FAIL 着色（绿/红）
  - 右轴：积分 A(ξ) (eV) + ±σ_A 不确定度带
  - ξ 降序时自动翻转 x 轴（`invert_xaxis`）
  - 标题：`ΔA = ... ± ... eV (ALL PASS / N FAILED)`
- DPI: 180

#### Constant-Potential Correction 输出契约

##### `write_corrected_free_energy_csv(result)` — 修正后自由能曲线

- 默认文件名：`ti_corrected_free_energy.csv`
- CSV 列：`xi, weight, dA_dxi, sem, A_const_q_eV, sigma_A_cumulative_eV, sigma_uC_cm2, phi_V_SHE, correction_eV, A_const_phi_eV`

##### 修正后自由能图（PNG）

- **自由能图**（`ti_corrected_free_energy.png`）：双轴图
  - 左轴：dA/dξ with error bars (a.u.)
  - 右轴：A_const_q (C1 橙, 虚线 + 误差带) + A_const_phi (C3 红, 实线)
  - 标题：`ΔA(const-q) = ... ± ... eV | ΔA(const-Φ) = ... eV`
- DPI: 180

### `slowgrowth_analysis(...)` 输出

- 统一入口：解析 restart + log 文件，切片/反转，输出 CSV + PNG
- 返回 `dict[str, Path]`，键可包含 `"csv"`、`"quick_png"`、`"publication_png"`

### CSV 输出

- 默认文件名：`slowgrowth_data.csv`（`DEFAULT_SG_CSV_NAME`）
- CSV 列（header 与顺序固定）：`step,time_fs,target_au,lagrange_au,free_energy_au,free_energy_ev`
- `free_energy_ev = free_energy_au * HA_TO_EV`

### PNG 输出

- Quick plot：`slowgrowth_quick.png`（`DEFAULT_SG_QUICK_PNG_NAME`）
  - 双轴图：左轴 Lagrange 乘子 (a.u.)，右轴自由能 (eV)
  - 顶轴：MD 步编号
  - 标注：barrier peak + 总自由能变
- Publication plot：`slowgrowth_publication.png`（`DEFAULT_SG_PUBLICATION_PNG_NAME`）
  - 双轴图：同上布局，legend 标注能量值

### 积分公式与单位

$$
\Delta A_k = -\sum_{i=0}^{k-1} \frac{\lambda_i + \lambda_{i+1}}{2} \Delta\xi
$$

- $\lambda_i$：Lagrange 乘子（a.u.）
- $\Delta\xi$：CV 步长（a.u./step）= `target_growth_au`（`Slowgrowth` 中已转换为 per-step；原始 `ConstraintInfo.target_growth_au` 为 per a.u. time，需乘 `dt_au = timestep_fs / AU_TIME_TO_FS`）
- 内部计算单位：Hartree；输出同时提供 Hartree 和 eV（`HA_TO_EV` 转换）

---

## 标定模块（`electrochemical/calibration/`）

### 输入

- CSV 文件：列 1 = 电势 φ (V vs SHE)，列 2 = 表面电荷密度 σ (μC/cm²)
  - 自动检测首行是否为表头（`_detect_header`：尝试 float 转换）
- 或手动 `data_points: list[tuple[float, float]]`（(φ, σ) 对）

### 标定 JSON（`calibration.json`）

默认位置：`~/.config/md_analysis/calibration.json`

```json
{
  "version": 1,
  "created": "ISO 8601",
  "reference": "SHE",
  "metadata": {},
  "data": {
    "potentials_V": [float, ...],
    "charge_densities_uC_cm2": [float, ...]
  },
  "fit": {
    "method": "linear|polynomial|spline",
    "params": { ... },
    "r_squared": float,
    "rmse": float,
    "equation": "str"
  }
}
```

### 输出

- CSV：`calibration_data.csv`（列 `potential_V`, `sigma_uC_cm2`）
- PNG：`calibration_fit.png`（散点 + 拟合曲线 + R²/RMSE 标注，11×4.8 inch, 160 DPI）

### Mapper 算法

- **LinearMapper**：`φ = slope · σ + intercept`，最小二乘拟合
- **PolynomialMapper**：`φ = Σ c_i · σ^i`，可配置阶数（默认 2）
- **SplineMapper**：按 σ 排序后构造三次样条插值（需 scipy）
- **DifferentialCapacitanceMapper**：
  1. 按 φ 升序排列标定点
  2. 计算相邻点间微分电容 `C_i = Δσ_i / Δφ_i` (μF/cm²)
  3. 分段线性插值；超出数据范围时使用最近端点的电容外推
  - 约束：σ 必须随 φ 单调递增（正微分电容），否则抛出 `ValueError`
  - 重复 φ 值同样抛出 `ValueError`

### 错误条件

- `ValueError`：`csv_path` 与 `data_points` 同时或均未提供
- `ValueError`：未知拟合方法或参考标度
- `ValueError`：PZC 转换缺少 `phi_pzc`
- `ValueError`：`DifferentialCapacitanceMapper` 重复 φ 值或 σ 非单调
- `RuntimeError`：mapper 未拟合即调用 `predict()`
- `ImportError`：`SplineMapper` 缺少 scipy
- `FileNotFoundError`：标定 JSON 不存在

### 电势参考转换

- SHE ↔ RHE：`φ_RHE = φ_SHE + (RT/F)·ln(10)·pH`
- SHE ↔ PZC：`φ_PZC = φ_SHE − φ_pzc`
