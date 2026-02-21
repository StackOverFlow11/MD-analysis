# Data Contract（跨模块核心契约）

> 仅记录跨模块共享、且会影响结果可比性的核心数据契约。
>
> 接口暴露清单与实现细则请查看：
> - `history/architecture/modules/scripts/structure/interface_exposure.md`
> - `history/architecture/modules/scripts/structure/utils/interface_exposure.md`
> - `history/architecture/modules/scripts/structure/utils/implementation_guidelines.md`
>
> 记录落位硬约束：本文件仅承载“全局契约”；非全局目录级记录不得写入本文件，
> 必须写入对应目录的 `interface_exposure.md` / `implementation_guidelines.md`。

## 适用范围

- 处理对象：单帧 `ase.Atoms`
- 原子索引：统一使用 **0-based**
- 分数坐标：默认按 wrap 后区间处理（`[0, 1)`）

## 核心数据载体（输出形状契约）

- 水分子标记：`(n_water, 3)`，每行 `[O_idx, H1_idx, H2_idx]`
- 水氧索引：`(n_water, 1)`（输入可接受 `(n,)` 或 `(n, 1)`）
- 质量密度 z 分布：`(nbins, 1)`，单位 `g/cm^3`
- 取向加权 z 分布：`(nbins, 1)`，单位 `1/Angstrom^3`
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
P_{\mathrm{orient},\mathrm{bin}} = \frac{\sum_{i\in \mathrm{bin}}\cos\theta_i}{V_{\mathrm{bin}}}
$$

### 3) 指定 c 窗口的角度 PDF

- 角域固定 `[0, 180]` 度
- 要求 `180 / ndeg` 为整数
- 窗口内无样本时，返回全零数组

## 默认参数来源

- 默认值统一来自：`scripts/structure/utils/config.py`
- 常用项包括：`DEFAULT_Z_BIN_WIDTH_A`、`DEFAULT_THETA_BIN_DEG`、`DEFAULT_WATER_OH_CUTOFF_A`
- 当前默认值（与代码同步）：
  - `DEFAULT_Z_BIN_WIDTH_A = 0.1` Angstrom
  - `DEFAULT_THETA_BIN_DEG = 5.0` degree
  - `DEFAULT_WATER_OH_CUTOFF_A = 1.25` Angstrom

## Analysis 层输出契约（新增）

### `water_mass_density_z_distribution_analysis(...)` CSV

- 输出文件默认参数来自：`scripts/structure/Analysis/config.py`
- 界面定义：使用直接面向非金属环境的 interface 对（每侧一层）
- CSV 列定义：
  - `path_fraction_center`：界面 -> 两界面中点路径的归一化坐标（`[0, 1]`）
  - `distance_A`：对应平均物理距离（Angstrom）
  - `rho_ensemble_avg_g_cm3`：系综平均水质量密度（`g/cm^3`）

### `water_orientation_weighted_density_z_distribution_analysis(...)` CSV

- 输出文件默认参数来自：`scripts/structure/Analysis/config.py`
- 界面定义：使用直接面向非金属环境的 interface 对（每侧一层）
- CSV 列定义：
  - `path_fraction_center`：界面 -> 两界面中点路径的归一化坐标（`[0, 1]`）
  - `distance_A`：对应平均物理距离（Angstrom）
  - `orientation_ensemble_avg_1_A3`：系综平均取向加权密度（`1/Angstrom^3`）

### `ad_water_orientation_analysis(...)` 输出契约

- 输出文件默认参数来自：`scripts/structure/Analysis/config.py`
- 吸附层主峰定义：密度分布最高 bin 的距离位置
- 吸附层区间：
  - 下界：主峰前最后一个近零点
  - 上界：主峰后第一个局部极小值（平滑后）
- 结果文件：
  - profile CSV：`distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_1_A3,is_adsorbed_layer_bin`
  - range TXT：`adsorbed_layer_start_A, adsorbed_layer_end_A, main_peak_distance_A, ...`

### `compute_adsorbed_water_theta_distribution(...)` 输出契约

- 统计对象：吸附层区间内水分子的取向角 `theta`
- 角域：`0-180` degree
- 输出 CSV：`theta_degree,pdf_degree_inv`

### `plot_water_three_panel_analysis(...)` 输出契约

- 入口定位：集成三联图的推荐入口（同时落盘各中间结果，便于复现）
- 输出目录：`output_dir`（默认 `Path.cwd()`）
- 输出文件（默认命名来自：`scripts/structure/Analysis/config.py`；PNG 名可覆盖）：
  - PNG：`DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME`（参数 `output_png_name` 可覆盖）
  - density CSV：`DEFAULT_WATER_MASS_DENSITY_CSV_NAME`
  - orientation CSV：`DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME`
  - adsorbed profile CSV：`DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME`
  - adsorbed range TXT：`DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME`
  - theta distribution CSV：`DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME`
- CSV 列（header 与顺序固定）：
  - density CSV：`path_fraction_center,distance_A,rho_ensemble_avg_g_cm3`
  - orientation CSV：`path_fraction_center,distance_A,orientation_ensemble_avg_1_A3`
  - adsorbed profile CSV：`distance_A,rho_ensemble_avg_g_cm3,orientation_ensemble_avg_1_A3,is_adsorbed_layer_bin`
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
