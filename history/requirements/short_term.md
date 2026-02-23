# 近期诉求（持续更新）

> 维护要求：每次提出"最近要做什么/优先级变化/卡点"，都更新这里。

## 当前阶段目标（与仓库现状对齐）

- **目标**：围绕"周期性金属-水界面"体系，提供可复现的水相关剖面分析（CSV/PNG）与可复用 API。
- **当前已覆盖**（对应 `scripts/structure/**`）：
  - 单帧：金属界面层识别、H2O 拓扑识别、密度/取向/角度 PDF 等统计
  - 多帧：从选定界面到两界面中点的系综平均（A 口径）、吸附层自动识别、吸附层角度分布、三联图输出
- **当前未覆盖**：电势（cube）与电荷（Mulliken）等分析仍未实现（仅保留数据样例与测试脚本命名）。

## 已确认的体系前提（用户声明，值得记录）

- 三基矢**正交**的周期性体系
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中总会存在**两个表面**（两个界面）

## 当前已实现能力（按"可跑通的入口"）

- **推荐一键入口**：
  - `scripts.structure.Analysis.plot_water_three_panel_analysis(xyz_path, md_inp_path, ...)`
    - 输出：密度/取向 CSV、吸附层 profile CSV、吸附层 range TXT、吸附层角度分布 CSV、三联图 PNG
- **基础分析入口**：
  - `water_mass_density_z_distribution_analysis(...)`
  - `water_orientation_weighted_density_z_distribution_analysis(...)`
  - `ad_water_orientation_analysis(...)`
  - `compute_adsorbed_water_theta_distribution(...)`

## 近期任务清单（仍待补齐）

- **I/O（读取与标准化）**
  - 继续补齐 CP2K 常见输出读取：`md-*.ener` / `CHARGE.mulliken` / `md-POTENTIAL-*.cube`
  - 统一记录单位、时间步、采样间隔等元数据（当前仅解析 `md.inp` 的 `ABC [angstrom]`）
- **Analysis（扩展分析量）**
  - 电势相关：从 Hartree cube 得到沿界面法向的平均电势 profile（含对齐与零点约定）
  - 电荷相关：按元素/分组/分层统计 Mulliken 电荷（明确"输出什么"和"局限性说明"）
- **工程化（可复现与易用性）**
  - CLI/批处理入口（对整目录/多轨迹跑同一套分析）
  - 依赖与环境说明（固定最小依赖集合与安装方式）

## 关键口径：已在当前实现中落地（不是待讨论）

- **界面参考面**：取"直接面向非金属环境"的金属界面层（每侧固定 1 层，共 2 层），使用其 c 分数坐标的圆均值作为界面位置。
- **法向/方向**：默认沿晶胞 `c` 轴（`normal="c"`，参考方向为 `+c_unit`）。
- **水取向定义**：`theta` 定义为 H-O-H 角平分线与 `+c_unit` 的夹角；剖面使用 `cos(theta) * m_water` 做取向加权，单位为 `g/cm^3`。
- **取向加权密度单位**：`g/cm^3`（与质量密度同单位；公式为 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）。
- **CSV 列名**：取向列为 `orientation_ensemble_avg_g_cm3`（旧列名 `orientation_ensemble_avg_1_A3` 已废弃）。

## 仍需明确/未实现的问题（保留为后续决策）

- 电势零点/对齐策略是什么？（真空区、体相水、或平均为 0）
