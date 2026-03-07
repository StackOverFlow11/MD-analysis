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
  - Bader 电荷下游分析（`charge/BaderAnalysis.py`）：
    - 单帧表面电荷密度 `compute_frame_surface_charge(method=...)`，支持 `"counterion"`（反离子/溶质）和 `"layer"`（界面层净电荷）两种计算方法
    - 按帧指定原子索引提取净电荷（`trajectory_indexed_atom_charges`）
- **当前未覆盖**：按层/按元素电荷转移统计；Mulliken 电荷分析仍未实现。

## 已确认的体系前提（用户声明，值得记录）

- 三基矢**正交**的周期性体系
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中总会存在**两个表面**（两个界面）

## 当前已实现能力（按"可跑通的入口"）

- **CLI 入口**：
  - `md-analysis water --xyz md-pos-1.xyz --md-inp md.inp`
  - `md-analysis potential --cube-pattern "md-POTENTIAL-*.cube"`
  - `md-analysis charge --root-dir . --charge-method counterion|layer` （Bader 表面电荷密度时序分析，CSV + PNG 输出）
  - `md-analysis all --xyz ... --md-inp ... --cube-pattern ...`
- **编程入口**：
  - `md_analysis.main.run_water_analysis(xyz_path, md_inp_path, ...)`
  - `md_analysis.main.run_potential_analysis(cube_pattern=..., md_out_path=..., ...)`
  - `md_analysis.main.run_charge_analysis(output_dir=..., root_dir=..., ...)`
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
  - Bader 电荷下游分析（已实现表面电荷密度，`charge/BaderAnalysis.py`）：
    - ✅ 表面电荷密度（双方法）：
      - `method="counterion"`：排除水分子和金属原子，仅反离子/溶质物种净电荷贡献 σ
      - `method="layer"`：界面层金属原子净电荷求和 / 面积
      - CLI 通过 `--charge-method counterion|layer` 选择；输出目录按方法分离 `<outdir>/charge/<method>/`
    - ✅ 单帧原子净电荷提取：`frame_indexed_atom_charges` 传入 `(N,)` 索引，返回 `(N, 2)` 的索引+净电荷数组
    - ✅ 轨迹原子净电荷提取：`trajectory_indexed_atom_charges` 按帧传入 `(t, N)` 索引矩阵，返回 `(t, N, 2)` 的索引+净电荷数组（内部调用 `frame_indexed_atom_charges`）
    - ✅ 轨迹表面电荷密度时序：`trajectory_surface_charge` 逐帧计算表面电荷密度，返回 `(t, 2)` 的 μC/cm² 数组
    - ✅ 端到端表面电荷分析：`surface_charge_analysis` 输出 CSV（含累积平均）+ PNG，CLI 子命令 `md-analysis charge`
    - 按层/按元素电荷转移统计：分层聚合 `bader_net_charge`，输出每层各元素的平均净电荷（待实现）
    - 典型工作流：CP2K MD → 提取结构帧 → `generate_bader_workdir` 生成 VASP 工作目录 → VASP 单点 → Bader 分析 → `load_bader_atoms` → 表面电荷/电荷转移
  - Bader 工作目录生成（`scripts/BaderGen.py`）：
    - ✅ `generate_bader_workdir()`：从单帧 Atoms 生成完整 VASP 工作目录（POSCAR + INCAR + KPOINTS + POTCAR + script.sh）
    - ✅ POSCAR 通过 IndexMapper 生成，保留 XYZ↔POSCAR 索引映射
    - ✅ POTCAR 通过 vaspkit 103 自动生成（可选）
    - ✅ 提交脚本路径支持持久化配置（`~/.config/md_analysis/config.json`）
    - ✅ 多帧批量生成：`batch_generate_bader_workdirs(xyz_path, cell_abc, output_dir, *, frame_start/end/step, ...)`
    - ✅ CLI 支持：401（单帧）+ 402（批量），cell 来源支持 `.restart` 和 `md.inp`
    - ✅ RestartParser：`parse_abc_from_restart()` 从 CP2K `.restart` 文件解析正交 cell 参数
  - 持久化用户配置（`config.py`）：
    - ✅ `load_config`、`save_config`、`get_config`、`set_config`
    - ✅ CLI 设置菜单（901/902）支持查看和修改配置
  - Mulliken 电荷：按元素/分组/分层统计（优先级低于 Bader，待后续明确需求）
- **工程化（可复现与易用性）**
  - 依赖与环境说明（固定最小依赖集合与安装方式）

## 关键口径：已在当前实现中落地（不是待讨论）

- **界面参考面**：取"直接面向非金属环境"的金属界面层（每侧固定 1 层，共 2 层），使用其分数坐标的圆均值（`center_frac`）作为界面位置。界面层标签为 `"normal_aligned"` / `"normal_opposed"`。
- **法向/方向**：默认沿晶胞 `c` 轴（`normal="c"`，参考方向为 `+c_unit`）。自定义向量法向不支持。
- **水取向定义**：`theta` 定义为 H-O-H 角平分线与 `+c_unit` 的夹角；剖面使用 `cos(theta) * m_water` 做取向加权，单位为 `g/cm^3`。
- **取向加权密度单位**：`g/cm^3`（与质量密度同单位；公式为 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）。
- **CSV 列名**：取向列为 `orientation_ensemble_avg_g_cm3`（旧列名 `orientation_ensemble_avg_1_A3` 已废弃）。
- **电极电势**：`U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP`（cSHE 方法）。
