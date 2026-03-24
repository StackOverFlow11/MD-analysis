# 代码架构（当前仓库结构与数据流）

用途：记录"仓库现在是什么样、模块怎么分、数据怎么流动"，用于协作对齐与避免口径漂移。

## 适用体系前提（当前实现依赖的假设）

- 三基矢**正交**的周期性体系（`md.inp` 中通过 `ABC [angstrom] a b c` 提供正交盒长）
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中存在**两个界面**（两侧表面）

## 仓库根目录结构（当前）

- `src/`：标准 src-layout，可安装包 `md_analysis`（`pip install .`）
- `test/`：pytest 单元/集成测试
- `data_example/`：最小可复现实例数据（`data_example/potential/`）
- `context4agent/`：架构/契约/决策/需求上下文（本目录）
- `README.md`：快速入口与最小跑通方式
- `CLAUDE.md`：Claude Code 协作指引

## `src/md_analysis/` 包内分层（当前）

### 1) `md_analysis.utils`（单帧/底层）

输入：单帧 `ase.Atoms`、cube 文件、或原始数组

- `config.py`
  - 全局常量：`TRANSITION_METAL_SYMBOLS`、`DEFAULT_METAL_SYMBOLS`、默认分箱宽度/阈值
  - 单位换算：`HA_TO_EV`、`BOHR_TO_ANG`
  - cSHE 常量：`DP_A_H3O_W_EV`、`MU_HPLUS_G0_EV`、`DELTA_E_ZP_EV`
- `CubeParser.py`
  - Gaussian cube 文件 I/O、plane-averaged φ(z)、slab-averaged potential
- `StructureParser/` 子包（结构解析：周期聚类、金属层识别、水分子拓扑）
  - `ClusterUtils.py`：1D 周期性聚类 + 最大间隙检测 + 间隙中点计算
  - `LayerParser.py`
    - 金属原子筛选 → 沿法向投影做 1D 聚类成"层"
    - 在周期水–金属–水体系下标记**两侧**直接面向环境的界面层（每侧固定 1 层，共 2 层）
    - 给界面层附带 `normal_unit`（从金属指向环境一侧）
  - `WaterParser.py`
    - 以 MIC 距离 + O–H 截断推断水分子拓扑，输出 `(n_water, 3)` 的 `[O, H1, H2]`
    - 基于氧索引计算：
      - 水质量密度 $\rho(z)$（`g/cm^3`）
      - 取向加权密度（`g/cm^3`，按 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）
      - 指定 c 分数窗口内的 $\theta$ PDF（`degree^-1`）
- `RestartParser/` 子包（CP2K restart 文件解析）
  - `CellParser.py`：CP2K cell 参数解析（`.restart` + `md.inp`）
  - `ColvarParser.py`：CP2K COLVAR restart + LagrangeMultLog 解析

> 低层 shape/单位/窗口规则以 `context4agent/architecture/modules/data_contract.md` 为准。

### 2) `md_analysis.water`（多帧/水分析工作流）

输入：`md-pos-*.xyz` 轨迹 + `md.inp`（用于解析正交盒长）

- `water/config.py`：水分析层默认参数（输出目录、文件名常量）
- `water/Water.py`
  - 集成三联图入口 `plot_water_three_panel_analysis(...)`
  - 将密度、取向、吸附层与 $\theta$ PDF 组合输出 PNG
- `water/WaterAnalysis/_common.py`（私有）
  - 从 `md.inp` 解析 `ABC [angstrom]`
  - 用 ASE `iread` 逐帧读取 xyz，并为每帧设置 cell + PBC
  - 每帧自动检测两侧界面位置（沿 c 分数坐标）
  - 计算"选定界面 → 两界面中点"的**半路径**上的密度/取向剖面
  - 将不同帧的剖面重采样到统一归一化网格后做等权系综平均（A 口径）
- `water/WaterAnalysis/WaterDensity.py`
  - 调用 `_common` 计算系综平均密度，并导出 CSV
- `water/WaterAnalysis/WaterOrientation.py`
  - 调用 `_common` 计算系综平均取向加权密度，并导出 CSV
- `water/WaterAnalysis/AdWaterOrientation.py`
  - 从密度剖面自动识别吸附层区间
  - 导出吸附层 profile（密度/取向/掩码）与区间文本
  - 二次遍历轨迹统计吸附层内的 $\theta$ 分布并导出 CSV

### 3) `md_analysis.electrochemical`（电化学分析分组包）

分组命名空间，包含 `potential`、`charge` 和 `calibration` 三个子包。自身不含业务逻辑，仅做 re-export。

顶层 `md_analysis` 通过 `from .electrochemical import potential, charge` 提供向后兼容的便捷访问。`calibration` 需通过 `from md_analysis.electrochemical import calibration` 访问。

### 3a) `md_analysis.electrochemical.charge`（多帧/电荷分析工作流）

输入：根目录下 `bader_t*_i*` 子目录，每帧各含 POSCAR + ACF.dat + POTCAR

- `electrochemical/charge/config.py`：单位换算常量（`E_PER_A2_TO_UC_PER_CM2`）、默认文件名
- `electrochemical/charge/Bader/` 子包：
  - `_frame_utils.py`：帧目录发现、排序、step/time 提取
  - `BaderData.py`：`BaderTrajectoryData` frozen dataclass + `load_bader_trajectory()` — 加载轨迹并通过 IndexMap remap 回 XYZ 原子序
  - `SurfaceCharge.py`：
    - `compute_frame_surface_charge(method=...)`：单帧表面电荷密度（`method="counterion"` 或 `"layer"`；结果存入 `atoms.info`）
    - `trajectory_surface_charge(method=...)`：多帧表面电荷密度时序 → `(t, 2)` ndarray
    - `surface_charge_analysis(method=...)`：端到端表面电荷密度分析（CSV + PNG 输出，含累积平均；若存在标定文件则自动追加外推电势列和右轴）
  - `AtomCharges.py`：
    - `frame_indexed_atom_charges()`：单帧指定原子索引提取净电荷 → `(N, 2)` ndarray（POSCAR 序）
    - `trajectory_indexed_atom_charges()`：按帧指定原子索引提取净电荷 → `(t, N, 2)` ndarray（POSCAR 序）
    - `tracked_atom_charge_analysis()`：指定 XYZ 原子追踪 Bader 净电荷含时演化 + 系综平均（CSV + PNG）
    - `counterion_charge_analysis()`：逐帧自动检测 counterion，输出含时演化 + 系综平均（CSV + PNG）

### 3b) `md_analysis.electrochemical.potential`（多帧/电势分析工作流）

输入：cube 文件 glob pattern + CP2K `md.out` + 可选 xyz 轨迹

- `electrochemical/potential/config.py`：电势分析输出文件名常量
- `electrochemical/potential/CenterPotential.py`
  - `center_slab_potential_analysis()`：逐帧 slab-averaged Hartree potential
  - `fermi_energy_analysis()`：从 md.out 提取 Fermi 能级时间序列
  - `electrode_potential_analysis()`：完整 U vs SHE 流程（调用以上二者并合并）
  - `thickness_sensitivity_analysis()`：扫描 slab 厚度，双轴图（U vs SHE + spatial std φ(z)）
- `electrochemical/potential/PhiZProfile.py`
  - `phi_z_planeavg_analysis()`：所有帧的 φ(z) overlay 可视化

### 3c) `md_analysis.electrochemical.calibration`（σ→φ 标定映射）

输入：(σ, φ) 数据点对（CSV 或手动输入），电势始终 V vs SHE

- `calibration/config.py`：Nernst 常量、默认文件名、拟合方法常量
- `calibration/_data.py`：`CalibrationData` 数据类、CSV/JSON I/O（自动检测表头）
- `calibration/_mapper.py`：`ChargePotentialMapper` ABC + `LinearMapper`/`PolynomialMapper`/`SplineMapper`
- `calibration/_plot.py`：标定拟合可视化（散点 + 拟合曲线 + R²/RMSE 标注）
- `calibration/CalibrationWorkflow.py`：
  - `calibrate()`：完整标定流程（加载数据→拟合→保存 JSON→出图）
  - `predict_potential()`：从保存的标定预测电势
  - `convert_reference()`：SHE ↔ RHE ↔ PZC 转换

### 5) `md_analysis.exceptions`（异常基类）

- `exceptions.py`：`MDAnalysisError(Exception)` — 所有领域异常的公共基类
  - 调用方可用 `except MDAnalysisError` 统一捕获全部领域异常

### 6) `md_analysis.config`（持久化用户配置）

- `config.py`：用户偏好持久化（`~/.config/md_analysis/config.json`）
  - `load_config()`、`save_config()`、`get_config()`、`set_config()`、`delete_config()`
  - 配置键常量：`KEY_VASP_SCRIPT_PATH`、`KEY_LAYER_TOL_A`、`KEY_Z_BIN_WIDTH_A`、`KEY_THETA_BIN_DEG`、`KEY_WATER_OH_CUTOFF_A`
  - 电势输出配置键：`KEY_POTENTIAL_REFERENCE`、`KEY_POTENTIAL_PH`、`KEY_POTENTIAL_TEMPERATURE_K`、`KEY_POTENTIAL_PHI_PZC`
  - `CONFIGURABLE_DEFAULTS`：可配置分析参数注册表（键 → 硬编码默认值 + 标签）

### 7) `md_analysis.main` / `md_analysis.cli`（集成入口）

- `main.py`：编程入口 `run_water_analysis()`、`run_potential_analysis()`、`run_charge_analysis()`、`run_all()`
- `cli/`：VASPKIT 风格交互式 CLI 包，注册为 `md-analysis` console script
  - `__init__.py`：`main()` 入口 + banner + 顶层菜单分发
  - `_prompt.py`：可复用的输入提示辅助函数（含 `_get_effective_default` 用于读取用户配置或硬编码默认值）
  - `_framework.py`：核心框架（`MenuNode`、`MenuGroup`、`MenuCommand`、`lazy_import()`）
  - `_params.py`：参数采集层次（`K` 键常量、`ParamCollector` ABC、泛型参数类）
  - `_water.py`：水分析命令类（101-105）：`WaterDensityCmd`、`WaterOrientationCmd`、`AdWaterOrientationCmd`、`AdWaterThetaCmd`、`WaterThreePanelCmd`
  - `_potential.py`：电势分析命令类（211-216）：`CenterPotentialCmd`、`FermiEnergyCmd`、`ElectrodePotentialCmd`、`PhiZProfileCmd`、`ThicknessSensitivityCmd`、`FullPotentialCmd`
  - `_charge.py`：电荷分析命令类（221-223）：`SurfaceChargeCmd`（通过 `method` 参数区分）、`_print_ensemble_summary()`
  - `_enhanced_sampling.py`：慢增长命令类（30=SG 301-302）：`SGQuickPlotCmd`、`SGPublicationPlotCmd`
  - `_constrained_ti.py`：约束 TI 命令类（31=TI 311-312）：`TISingleDiagCmd`、`TIFullAnalysisCmd`
  - `_scripts.py`：脚本/工具命令类（41=Bader 411-412，42=TI 421-422）：`BaderSingleCmd`、`BaderBatchCmd`、`TISingleCmd`、`TIBatchCmd`
  - `_settings.py`：设置命令类（901-909）：`SetVaspScriptCmd`、`SetCp2kScriptCmd`、`ShowConfigCmd`、`SetAnalysisDefaultCmd`、`ResetDefaultsCmd`、`SetPotentialReferenceCmd`

### 8) `md_analysis.enhanced_sampling`（增强抽样分析工作流）

- `enhanced_sampling/__init__.py`：空包（re-export 预留）
- `enhanced_sampling/slowgrowth/`：慢增长热力学积分子包
  - `SlowGrowth.py`：数据类（`Slowgrowth`、`SlowgrowthFull`、`SlowgrowthSegment`）+ 自由能积分（midpoint rule，CP2K 约定取负号）
  - `SlowGrowthPlot.py`：双轴绘图（quick / publication）+ CSV 导出 + 统一入口 `slowgrowth_analysis`
  - `config.py`：输出文件名常量
- `enhanced_sampling/constrained_ti/`：约束 TI 收敛诊断子包
  - `workflow.py`：编排器（`analyze_standalone`、`analyze_ti`、`standalone_diagnostics`、CSV 导出）
  - `io.py`：自动发现约束点目录（`discover_ti_points`、`load_ti_series`）
  - `plot.py`：2×2 诊断图 + 自由能曲线图
  - `integration.py`：梯形积分权重、SEM targets、自由能积分
  - `analysis/`：四步诊断引擎（ACF、F&P block average、running average、Geweke）
- 依赖方向：`enhanced_sampling` → `utils`（`ColvarParser`、`config`、`_io_helpers`）；无反向依赖
- CLI 集成：30=Slow-Growth（301-302，`cli/_enhanced_sampling.py`），31=Constrained TI（311-312，`cli/_constrained_ti.py`）

### 9) `md_analysis.scripts`（自动化脚本工具）

- `scripts/__init__.py`：re-exports `BaderGenError`, `generate_bader_workdir`, `batch_generate_bader_workdirs`
- `scripts/BaderGen.py`：从单帧 MD 结构生成 VASP Bader 工作目录（POSCAR + INCAR + KPOINTS + POTCAR + script.sh）
- `scripts/template/`：VASP 模板文件（INCAR、KPOINTS），通过 `importlib.resources` 访问
- `scripts/utils/IndexMapper.py`
  - CP2K XYZ ↔ VASP POSCAR 双射索引映射
  - `compute_index_map()`：计算排列映射
  - `write_poscar_with_map()`：写入带映射注释的 POSCAR
  - `read_index_map_from_poscar()`：从 POSCAR 注释行解码映射
  - `encode_comment_line()` / `decode_comment_line()`：注释行编解码
  - `remap_array()`：按映射重排 per-atom 数组

## 核心数据流（当前端到端）

### 水分析

1. 解析 `md.inp` 获取正交盒长 $a,b,c$
2. 逐帧读取 xyz → 设置 cell + PBC
3. 逐帧检测两侧界面 c 分数坐标（来自金属界面层）
4. 逐帧识别水分子 → 选定一侧界面起点 → 在半路径上分箱统计
5. 将每帧剖面映射到统一归一化坐标并等权平均
6. 导出 CSV；若走三联图入口，再额外计算吸附层与 $\theta$ 分布并输出 PNG

### 电势分析

1. 从 cube 文件 glob 获取帧列表（按 step 排序）
2. 逐帧读取 cube → plane-average → slab-average（以界面或 cell 中心为参考）
3. 从 `md.out` 提取 Fermi 能级 → 合并计算 U vs SHE
4. 可选：thickness sensitivity 扫描、φ(z) overlay
5. 导出 CSV + PNG

## 文档镜像与落位

- `context4agent/architecture/modules/src/` 必须镜像 `src/md_analysis/` 的目录结构，并在每个镜像目录维护：
  - `interface_exposure.md`：对外公开符号/导入方式/兼容承诺
  - `implementation_guidelines.md`：职责边界/依赖方向/契约同步要求
