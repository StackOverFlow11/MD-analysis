# 近期诉求（持续更新）

> 维护要求：每次提出"最近要做什么/优先级变化/卡点"，都更新这里。

## 当前阶段目标（与仓库现状对齐）

- **目标**：围绕"周期性金属-水界面"体系，提供可复现的水/电势分析（CSV/PNG）与可复用 API。
- **当前已覆盖**（对应 `src/md_analysis/`）：
  - 单帧工具（`utils/`）：金属界面层识别、H2O 拓扑识别、密度/取向/角度 PDF、cube 文件解析、slab-averaged potential
  - 水分析（`water/`）：从选定界面到两界面中点的系综平均（A 口径）、吸附层自动识别、吸附层角度分布、三联图输出
  - 电势分析（`electrochemical/potential/`）：center slab potential、Fermi energy、electrode potential U vs SHE、φ(z) overlay、thickness sensitivity
  - 增强抽样（`enhanced_sampling/`）：慢增长自由能绘图（quick / publication）+ CSV 导出 + 约束 TI 收敛诊断与自由能积分 + CLI 集成
  - 集成入口：CLI（`md-analysis` 命令）、编程入口（`md_analysis.workflows`，30 个 `run_*` 函数返回 `WorkflowResult`；`main.py` 是同名薄 re-export facade）、Agent 入口（`agent/`：dispatch + JSON Schema + TaskResult + Tools-layer `TaskContract`；入口重构 Phase 3 / 5B 删了 6 个 legacy task，**剩余 8 个任务全部带完整 contract**）
- Bader 电荷解析（`utils/formats/bader.py`）：从 VASP Bader 输出（ACF.dat + POTCAR）读取原始电子数与净电荷，附加到 ASE Atoms
  - Bader 电荷下游分析（`electrochemical/charge/Bader/`）：
    - 核心数据结构 `BaderTrajectoryData` + `load_bader_trajectory()` — 加载轨迹并通过 IndexMap remap 回 XYZ 原子序
    - 单帧表面电荷密度 `compute_frame_surface_charge(method=...)`，支持 `"counterion"`（反离子/溶质）和 `"layer"`（界面层净电荷）两种计算方法
    - 按帧指定原子索引提取净电荷（`trajectory_indexed_atom_charges`）
    - 指定 XYZ 原子电荷追踪（`tracked_atom_charge_analysis`）— 含时演化 + 系综平均
    - Counterion 逐帧自动检测追踪（`counterion_charge_analysis`）— 含时演化 + 系综平均
- **当前未覆盖**：按层/按元素电荷转移统计；Mulliken 电荷分析仍未实现。

## 已确认的体系前提（用户声明，值得记录）

- 三基矢**正交**的周期性体系
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中总会存在**两个表面**（两个界面）

## 当前已实现能力（按"可跑通的入口"）

- **CLI 入口**（`md-analysis` 启动 VASPKIT 风格交互式编号菜单，无 argparse 参数）：
  - 1xx：Water（密度/取向/吸附层/三联图）
  - 21x：Potential（center potential / Fermi / electrode potential / φ(z) / thickness sensitivity / full）
  - 22x：Charge（221 counterion σ / 222 layer σ / 223 full σ / 224 single-side σ+φ / 225 tracked atoms / 226 counterion tracking）
  - 23x：Calibration（CSV/手动标定 + 预测）
  - 30x：Slow-Growth（301 quick plot / 302 publication plot）
  - 31x：Constrained TI（311 single-point diagnostics / 312 full analysis / 313 constant-potential correction）
  - 41x：Bader 工作目录生成（411 单帧 / 412 批量）
  - 42x：TI 工作目录生成（421 单帧 / 422 批量）
  - 9xx：Settings（配置管理）
- **编程入口**（canonical 模块：`md_analysis.workflows`；同名 re-export 在 `md_analysis.main`；旧 `run_*_analysis` / `run_all` 已在入口重构 Phase 7a 移除）：
  - 水：`run_water_three_panel`（composite）+ `run_water_density` / `run_water_orientation` / `run_ad_water_orientation` / `run_ad_water_theta` 4 个单步入口
  - 电势：`run_potential_full`（composite）+ `run_center_potential` / `run_fermi_energy` / `run_electrode_potential` / `run_phi_z_profile` / `run_thickness_sensitivity` 5 个单步入口
  - 表面电荷：`run_surface_charge(output_dir, root_dir, method=..., ...) -> WorkflowResult`
  - 追踪原子电荷：`run_tracked_charge(output_dir, root_dir, atom_indices_xyz, ...) -> WorkflowResult`
  - 反离子电荷：`run_counterion_charge(output_dir, root_dir, ...) -> WorkflowResult`
  - composite 水+电势：`run_interface_analysis(xyz_path, md_inp_path, output_dir, ...) -> WorkflowResult`
  - 标定：`run_calibration_fit(...)` / `run_calibration_predict(...)`
  - 增强采样：`run_slowgrowth_quick_plot` / `run_slowgrowth_publication_plot` / `run_ti_single_diagnostics` / `run_ti_full_analysis` / `run_ti_constant_potential_correction`
  - 脚本生成：`run_bader_single/batch` / `run_ti_single/batch` / `run_potential_single/batch` / `run_sp_single/batch`
  - 共计 30 个 `run_*` + `WorkflowResult`，详见 `src/md_analysis/workflows/__init__.py`
- **Agent 入口**（`md_analysis.agent`，非交互式，面向 AI agent / MCP Server）：
  - `dispatch(task, params)` → 统一任务执行，返回 `TaskResult`
  - `list_tasks()` → 枚举已注册任务
  - `get_task_schema(task)` → 有 contract 时从 `TaskContract.to_agent_schema()` 生成（权威），否则从 `target_fn` 签名推导（legacy 路径）。返回形状始终为 OpenAI function-calling 兼容的 `{name, description, parameters}`
  - Tools-layer 契约（`_contracts.py`）：`TaskContract`（`inputs` + 三分法 `outputs_artifacts/metrics/raw_model` + `preconditions` + `side_effects` + `exceptions`）、`FieldSpec`（`json_schema` 权威，`unit`/`shape`/`path_kind`/`category` 领域标注）、`ExceptionMapping`（FQN + `error_type` 重写 dispatch 分类，有序先具体后父类）
  - 入口重构 Phase 3 / 5B 删除了 6 个 legacy task（water_three_panel, potential_full, charge_surface, charge_tracked, charge_counterion, run_all）—— 对应业务直接通过 `md_analysis.workflows.run_*` 调用；剩余 **8 个任务全部带完整 contract**：calibration_fit_csv, calibration_predict, slowgrowth_quick, ti_full_analysis（composite wrapper, 由 `run_ti_full_from_root` 支撑）, bader_gen_batch, ti_gen_batch, sp_gen_batch, config_show（read-only）
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
  - Bader 电荷下游分析（`electrochemical/charge/Bader/`）：
    - ✅ 表面电荷密度（双方法）：
      - `method="counterion"`：排除水分子和金属原子，仅反离子/溶质物种净电荷贡献 σ
      - `method="layer"`：界面层金属原子净电荷求和 / 面积（`n_surface_layers` 参数控制每侧取几层，默认 1）
      - CLI 通过交互式菜单选择（221=counterion / 222=layer / 223=prompted）；输出目录按方法分离 `<outdir>/electrochemical/charge/<method>/`
    - ✅ 单帧原子净电荷提取：`frame_indexed_atom_charges` 传入 `(N,)` 索引，返回 `(N, 2)` 的索引+净电荷数组
    - ✅ 轨迹原子净电荷提取：`trajectory_indexed_atom_charges` 按帧传入 `(t, N)` 索引矩阵，返回 `(t, N, 2)` 的索引+净电荷数组（内部调用 `frame_indexed_atom_charges`）
    - ✅ 轨迹表面电荷密度时序：`trajectory_surface_charge` 逐帧计算表面电荷密度，返回 `(t, 2)` 的 μC/cm² 数组
    - ✅ 端到端表面电荷分析：`surface_charge_analysis` 输出 CSV（含累积平均）+ PNG；若存在标定文件（`calibration.json`）则自动追加外推电势列及 PNG 右轴
    - 按层/按元素电荷转移统计：分层聚合 `bader_net_charge`，输出每层各元素的平均净电荷（待实现）
    - 典型工作流：CP2K MD → 提取结构帧 → `generate_bader_workdir` 生成 VASP 工作目录 → VASP 单点 → Bader 分析 → `load_bader_atoms` → 表面电荷/电荷转移
  - Bader 工作目录生成（`scripts/BaderGen.py`）：
    - ✅ `generate_bader_workdir()`：从单帧 Atoms 生成完整 VASP 工作目录（POSCAR + INCAR + KPOINTS + POTCAR + script.sh）
    - ✅ POSCAR 通过 IndexMapper 生成，保留 XYZ↔POSCAR 索引映射
    - ✅ POTCAR 通过 vaspkit 103 自动生成（可选）
    - ✅ 提交脚本路径支持持久化配置（`~/.config/md_analysis/config.json`）
    - ✅ 多帧批量生成：`batch_generate_bader_workdirs(xyz_path, cell_abc, output_dir, *, frame_start/end/step, ...)`
    - ✅ CLI 支持：411（单帧）+ 412（批量），cell 来源支持 `.restart` 和 `md.inp`
    - ✅ `utils/formats/cp2k_cell.py`：`parse_abc_from_restart()` / `parse_abc_from_md_inp()` 从 CP2K `.restart` 或 `md.inp` 解析正交 cell 参数
    - ✅ `utils/formats/cp2k_colvar.py`：解析 CP2K COLVAR restart 元数据（COLLECTIVE、CONSTRAINT、FIXED_ATOMS）和 LagrangeMultLog（单约束/多约束自动检测），重建 ξ(t) 目标序列；溢出值 `***` 自动处理为 `np.nan`。**当前 canonical 入口在 `md_analysis.engines.cp2k`**（`CP2KParser.parse_metadata` / `parse_lambda_series`，配 `ConstraintMetadata` / `LambdaSeries` engine-neutral dataclass），`utils/formats/cp2k_colvar` 是单文件解析层
  - 持久化用户配置（`config.py`）：
    - ✅ `load_config`、`save_config`、`get_config`、`set_config`、`delete_config`
    - ✅ `CONFIGURABLE_DEFAULTS` 注册表：`layer_tol_A`、`z_bin_width_A`、`theta_bin_deg`、`water_oh_cutoff_A`
    - ✅ CLI 设置菜单（901-907）支持查看/修改配置和分析参数默认值
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
