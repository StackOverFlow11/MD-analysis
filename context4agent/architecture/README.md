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
- `AGENTS.md`：Claude Code 协作指引

## `src/md_analysis/` 包内分层（当前）

### 1) `md_analysis.utils`（单帧/底层）

输入：单帧 `ase.Atoms`、cube 文件、或原始数组

- `constants.py`
  - 全局常量：`TRANSITION_METAL_SYMBOLS`、`DEFAULT_METAL_SYMBOLS`、默认分箱宽度/阈值
  - 单位换算：`HA_TO_EV`、`BOHR_TO_ANG`
  - cSHE 常量：`DP_A_H3O_W_EV`、`MU_HPLUS_G0_EV`、`DELTA_E_ZP_EV`
- `structure/` 子包（结构解析：周期聚类、金属层识别、水分子拓扑）
  - `cluster.py`：1D 周期性聚类 + 最大间隙检测 + 间隙中点计算
  - `layer.py`
    - 金属原子筛选 → 沿法向投影做 1D 聚类成"层"
    - 在周期水–金属–水体系下标记**两侧**直接面向环境的界面层（每侧固定 1 层，共 2 层）
    - 给界面层附带 `normal_unit`（从金属指向环境一侧）
  - `water.py`
    - 以 MIC 距离 + O–H 截断推断水分子拓扑，输出 `(n_water, 3)` 的 `[O, H1, H2]`
    - 基于氧索引计算：
      - 水质量密度 $\rho(z)$（`g/cm^3`）
      - 取向加权密度（`g/cm^3`，按 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）
      - 指定 c 分数窗口内的 $\theta$ PDF（`degree^-1`）
- `formats/` 子包（单文件格式解析，按引擎家族分子包；不做目录发现与业务逻辑）
  - `common/cube.py`：Gaussian cube 文件 I/O、plane-averaged φ(z)、slab-averaged potential
  - `cp2k/cell.py`：CP2K cell 参数解析（`.restart` + `md.inp`）
  - `cp2k/colvar.py`：CP2K COLVAR restart + LagrangeMultLog 解析（输出 CP2K raw dataclass，canonical 类型在 `engines.models`）
  - `cp2k/inp.py`：CP2K `.inp` parse/generate 封装（cp2k-input-tools cit 树；全包唯一直接依赖 cit/pint 的模块，惰性加载。cell/colvar 读侧与 scripts 写侧"parse → 树编辑 → generate"的共用底座）
  - `cp2k/stdout.py`：CP2K `md.out` / `sp.out` Fermi 能级提取（step/time/fermi_raw）
  - `cp2k/xyz.py`：CP2K xyz 轨迹流式读取（`i = <step>` 注释解析，按 step 构建 `ase.Atoms`）
  - `vasp/report.py` / `vasp/outcar.py` / `vasp/locpot.py`：VASP REPORT / OUTCAR / LOCPOT 占位解析器（全部 raise `NotImplementedError`）
  - `bader/acf.py`：Bader ACF.dat 解析并将净电荷附着到 ASE `Atoms`
  - `bader/potcar.py`：VASP POTCAR `(element, ZVAL)` 解析
  - `bader/_errors.py`：`BaderParseError` 共享异常
- `io/` 子包（目录发现与通用 I/O 脚手架，不做格式解析）
  - `_frame_discovery.py`：帧工作目录发现与排序（`bader_t*_i*` / `potential_t*_i*` 命名约定）
  - `_io_helpers.py`：私有共享辅助（CSV 写出、累积平均）

> 低层 shape/单位/窗口规则以 `context4agent/architecture/modules/data_contract.md` 为准。

### 1a) `md_analysis.engines`（引擎适配层）

计算引擎（CP2K / VASP）门面层：把各引擎的原始文件组合（CP2K `*.restart` + `*.LagrangeMultLog`、`md.out` + cube 等）解析成 engine-neutral 的 frozen dataclass；上层业务（`water` / `electrochemical` / `enhanced_sampling` / `scripts`）只面向中性类型与 Protocol 编程，不触碰引擎细节。

- `__init__.py`：公共 facade — re-export Protocol / registry / 中性 dataclass / CP2K 门面函数，并在包 import 时默认注册 CP2K parser（`register_parser("cp2k", CP2KParser)`）
- `protocols.py`：`ConstraintMDParser` Protocol（`is_constraint_directory` / `parse_metadata` / `parse_lambda_series`）+ parser 注册表与目录嗅探（`register_parser` / `get_parser` / `infer_parser` / `resolve_parser`）+ `ParserInferenceError`；新增引擎 = 实现 Protocol + 注册，上层零改动
- `models.py`：engine-neutral frozen dataclass
  - 约束 MD：`ConstraintInfo` / `ColvarInfo` / `ConstraintMetadata` / `LambdaSeries` / `ConstraintRun`（metadata + λ(t) 复合视图，派生 target/time series 走 `@property`）
  - 电势：`PotentialFrame` / `CenterPotentialScalarFrame` / `FermiRecord`
  - cell：`CellSpec`（3×3 矩阵 + `abc_ang` / `is_orthorhombic`）
- `cp2k.py`：`CP2KParser`（Protocol 实现）+ 模块级 `read_*` 门面函数（`read_constraint_metadata` / `read_lambda_series` / `read_constraint_run` / `read_fermi_series` / `read_cell` / `read_continuous_potential_frames` / `read_distributed_potential_frames` 等）
- `vasp.py`：`VASPParser` 占位类——实现 Protocol 形状但所有方法 raise `NotImplementedError`，不在默认注册表中（避免 `infer_parser` 静默分派到 stub）
- 依赖方向：`engines` → `utils/{formats, constants, io}` + `exceptions`；禁止反向 import 业务包

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
- `calibration/_mapper.py`：`ChargePotentialMapper` ABC + `LinearMapper`/`PolynomialMapper`/`SplineMapper`/`DifferentialCapacitanceMapper`
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

### 7) `md_analysis.workflows` / `md_analysis.main` / `md_analysis.cli`（集成入口）

- `workflows/`：编程入口的 **canonical** 位置。全部 `run_*` 函数 + `WorkflowResult` frozen dataclass（`artifacts` / `metadata` / `extra`）。**完整 `run_*` 清单以 `src/md_analysis/workflows/__init__.py.__all__` 为权威**。按域拆子模块（`water.py` / `potential.py` / `charge.py` / `calibration.py` / `enhanced_sampling.py` / `scripts.py` / `composite.py` + `models.py` + `__init__.py`）。每个 `run_*` 主动 `Path()` 规范化路径参数；artifacts 字典值始终是 `Path`，metadata 轻量 scalars，extra 携带强类型 report（TI / calibration 等）。
- `main.py`：import-only **re-export facade**。导出与 `workflows/` 完全一致的公开符号（`WorkflowResult` + 全部 `run_*`，清单以 `workflows/__init__.py.__all__` 为权威）。**入口重构期间删除了所有 legacy 名**：`run_water_analysis` / `run_potential_analysis` / `run_charge_analysis` / `run_tracked_charge_analysis` / `run_counterion_charge_analysis` / `run_all`；旧返回类型 `dict[str, Path]` → 新返回类型 `WorkflowResult`，**未保留 alias**（避免同名静默换返回类型）。新调用方式 `for k, p in result.artifacts.items()`。
- `cli/`：VASPKIT 风格交互式 CLI 包，注册为 `md-analysis` console script
  - `__init__.py`：`main()` 入口 + banner + 顶层菜单分发
  - `_prompt.py`：可复用的输入提示辅助函数
  - `_framework.py`：核心框架（`MenuNode`、`MenuGroup`、`MenuCommand`、`lazy_import()`）
  - `_params.py`：参数采集层次（`K` 键常量、`ParamCollector` ABC、泛型参数类）
  - `_water.py`：水分析命令类（101-105）：`WaterDensityCmd`、`WaterOrientationCmd`、`AdWaterOrientationCmd`、`AdWaterThetaCmd`、`WaterThreePanelCmd`
  - `_potential.py`：电势分析命令类（211-216）：`CenterPotentialCmd`、`FermiEnergyCmd`、`ElectrodePotentialCmd`、`PhiZProfileCmd`、`ThicknessSensitivityCmd`、`FullPotentialCmd`
  - `_charge.py`：电荷分析命令类（221-226）：`SurfaceChargeCmd`（通过 `method` 参数区分）、`SingleSideChargeCmd`、`TrackedChargeCmd`、`CounterionChargeCmd`、`_print_ensemble_summary()`
  - `_enhanced_sampling.py`：慢增长命令类（30=SG 301-302）：`SGQuickPlotCmd`、`SGPublicationPlotCmd`
  - `_constrained_ti.py`：约束 TI 命令类（31=TI 311-313）：`TISingleDiagCmd`、`TIFullAnalysisCmd`、`TIConstPotCorrectionCmd`
  - `_scripts.py`：脚本/工具命令类（41=Bader 411-412，42=TI 421-422，43=SP Potential 431-432，44=DeePMD SP 441-442）：`BaderSingleCmd`、`BaderBatchCmd`、`TISingleCmd`、`TIBatchCmd`、`PotentialSingleCmd`、`PotentialBatchCmd`、`SpGenSingleCmd`、`SpGenBatchCmd`
  - `_settings.py`：设置命令类（900-931）：`SetVaspScriptCmd`、`SetCp2kScriptCmd`、`SetSpInpTemplateCmd`、`SetDpSpInpTemplateCmd`、`ShowConfigCmd`、`SetAnalysisDefaultCmd`、`ResetDefaultsCmd`、`SetPotentialReferenceCmd`

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
- 依赖方向：`enhanced_sampling` → `utils`（`formats.cp2k_colvar`、`config`、`_io_helpers`）；无反向依赖
- CLI 集成：30=Slow-Growth（301-302，`cli/_enhanced_sampling.py`），31=Constrained TI（311-313，`cli/_constrained_ti.py`）

### 9) `md_analysis.scripts`（自动化脚本工具）

- `scripts/__init__.py`：re-exports `BaderGenError`, `generate_bader_workdir`, `batch_generate_bader_workdirs`, `PotentialGenError`, `generate_potential_workdir`, `batch_generate_potential_workdirs`, `SpGenError`, `generate_sp_workdir`, `batch_generate_sp_workdirs`, `TIGenError`, `generate_ti_workdir`, `batch_generate_ti_workdirs`
- `scripts/BaderGen.py`：从单帧 MD 结构生成 VASP Bader 工作目录（POSCAR + INCAR + KPOINTS + POTCAR + script.sh）
- `scripts/template/`：内置模板（VASP INCAR、KPOINTS；CP2K 默认电势模板 `sp_potential.inp`，PotentialGen 三级回退的第三级），通过 `importlib.resources` 访问
- `scripts/_inp_tree_edit.py`：CP2K inp 树编辑 helper（工作流语义层，建于 `utils.formats.cp2k.inp` 格式原语之上；TIGen / PotentialGen / SpGen 写侧共享）
- `scripts/_potential_template.py`：内置电势模板常量 + 树校验器（`extract_kind_elements` / `has_smear` / `has_v_hartree_cube`）
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
