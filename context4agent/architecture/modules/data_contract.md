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

### `compute_frame_surface_charge(atoms, ...)` 输出

- 结果存入 `atoms.info`（原地修改）
- `atoms.info["surface_charge_density_e_A2"]`：`[σ_bottom, σ_top]`，单位 e/Å²
- `atoms.info["surface_charge_density_uC_cm2"]`：`[σ_bottom, σ_top]`，单位 μC/cm²
- `normal` 参数控制面积计算：`_AREA_VECTORS = {"a": (1,2), "b": (0,2), "c": (0,1)}`

### `trajectory_indexed_atom_charges(root_dir, atom_index_matrix, ...)` 输出

- 输入 `atom_index_matrix`：`(t, N)` 0-based 整型数组
- 返回 `np.ndarray`：`(t, N, 2)`
  - `[:, :, 0]`：回显的原子索引
  - `[:, :, 1]`：对应的 Bader 净电荷（单位 e）

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
- 需要 `md.out`（Fermi 能级）；无 `md.out` 时跳过
