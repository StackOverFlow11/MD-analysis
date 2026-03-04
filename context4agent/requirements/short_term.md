# 近期诉求（持续更新）

> 维护要求：每次提出"最近要做什么/优先级变化/卡点"，都更新这里。

## 当前阶段目标（与仓库现状对齐）

- **目标**：围绕"周期性金属-水界面"体系，提供可复现的水/电势分析（CSV/PNG）与可复用 API。
- **当前已覆盖**（对应 `src/md_analysis/`）：
  - 单帧工具（`utils/`）：金属界面层识别、H2O 拓扑识别、密度/取向/角度 PDF、cube 文件解析、slab-averaged potential
  - 水分析（`water/`）：从选定界面到两界面中点的系综平均（A 口径）、吸附层自动识别、吸附层角度分布、三联图输出
  - 电势分析（`potential/`）：center slab potential、Fermi energy、electrode potential U vs SHE、φ(z) overlay、thickness sensitivity
  - 集成入口：CLI（`md-analysis` 命令）、编程入口（`main.py`）
- Bader 电荷解析（`utils/BaderParser.py`）：从 VASP Bader 输出（ACF.dat + POTCAR）读取原始电子数与净电荷，附加到 ASE Atoms
- **当前未覆盖**：基于 Bader 电荷的表面电荷密度分析、按层/按元素电荷转移统计；Mulliken 电荷分析仍未实现。

## 已确认的体系前提（用户声明，值得记录）

- 三基矢**正交**的周期性体系
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中总会存在**两个表面**（两个界面）

## 当前已实现能力（按"可跑通的入口"）

- **CLI 入口**：
  - `md-analysis water --xyz md-pos-1.xyz --md-inp md.inp`
  - `md-analysis potential --cube-pattern "md-POTENTIAL-*.cube"`
  - `md-analysis all --xyz ... --md-inp ... --cube-pattern ...`
- **编程入口**：
  - `md_analysis.main.run_water_analysis(xyz_path, md_inp_path, ...)`
  - `md_analysis.main.run_potential_analysis(cube_pattern=..., md_out_path=..., ...)`
  - `md_analysis.main.run_all(...)`
- **水分析**：
  - `plot_water_three_panel_analysis(xyz_path, md_inp_path, ...)`
    - 输出：密度/取向 CSV、吸附层 profile CSV、吸附层 range TXT、吸附层角度分布 CSV、三联图 PNG
  - `water_mass_density_z_distribution_analysis(...)`
  - `water_orientation_weighted_density_z_distribution_analysis(...)`
  - `ad_water_orientation_analysis(...)`
  - `compute_adsorbed_water_theta_distribution(...)`
- **电势分析**：
  - `center_slab_potential_analysis(...)`
  - `fermi_energy_analysis(...)`
  - `electrode_potential_analysis(...)`
  - `thickness_sensitivity_analysis(...)`
  - `phi_z_planeavg_analysis(...)`

## 近期任务清单（仍待补齐）

- **I/O（读取与标准化）**
  - 统一记录单位、时间步、采样间隔等元数据（当前仅解析 `md.inp` 的 `ABC [angstrom]`）
- **Analysis（扩展分析量）**
  - Bader 电荷下游分析（近期重点）：
    - 表面电荷密度：结合 `detect_interface_layers` 识别表面层，计算每层净电荷密度；内部计算用 e/Å²，最终输出转换为 μC/cm²（与实验量纲对齐）
    - 按层/按元素电荷转移统计：分层聚合 `bader_net_charge`，输出每层各元素的平均净电荷
    - 典型工作流：CP2K MD → 提取结构帧 → VASP 单点 → Bader 分析 → `load_bader_atoms` → 表面电荷/电荷转移
  - Mulliken 电荷：按元素/分组/分层统计（优先级低于 Bader，待后续明确需求）
- **工程化（可复现与易用性）**
  - 依赖与环境说明（固定最小依赖集合与安装方式）

## 关键口径：已在当前实现中落地（不是待讨论）

- **界面参考面**：取"直接面向非金属环境"的金属界面层（每侧固定 1 层，共 2 层），使用其 c 分数坐标的圆均值作为界面位置。
- **法向/方向**：默认沿晶胞 `c` 轴（`normal="c"`，参考方向为 `+c_unit`）。
- **水取向定义**：`theta` 定义为 H-O-H 角平分线与 `+c_unit` 的夹角；剖面使用 `cos(theta) * m_water` 做取向加权，单位为 `g/cm^3`。
- **取向加权密度单位**：`g/cm^3`（与质量密度同单位；公式为 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）。
- **CSV 列名**：取向列为 `orientation_ensemble_avg_g_cm3`（旧列名 `orientation_ensemble_avg_1_A3` 已废弃）。
- **电极电势**：`U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP`（cSHE 方法）。
