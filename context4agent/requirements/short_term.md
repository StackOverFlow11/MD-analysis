# 近期诉求（持续更新）

> 维护要求：每次提出"最近要做什么/优先级变化/卡点"，都更新这里。
> 本文只记「能力分组摘要 + 关键不变量 + 待办意图」；逐函数/逐菜单号
> **不手写**，以代码权威为准（见下「权威来源」）。

## 权威来源（不在本文复制，避免漂移）

- 编程入口完整清单：`src/md_analysis/workflows/__init__.py` 的 `__all__`
  （全部 `run_*` + `WorkflowResult`；`md_analysis.main` 是同名薄 re-export）
- CLI 菜单号与行为：`src/md_analysis/cli/` + `test/unit/cli/`
- 模块职责边界 / 数据流：`context4agent/architecture/`

## 当前阶段目标

- **目标**：围绕"周期性金属-水界面"体系，提供可复现的水/电势/电荷/
  增强采样分析（CSV/PNG）与可复用 API（CLI + workflows 双入口；后续走
  CLI + skill 路线）。
  > 注：非交互 `agent`（dispatch/TaskContract）入口已在 development
  > 移除，保留于 `development_agent` 分支（未来 MCP 路线）。
- **当前已覆盖**（按能力域，详清单见代码权威位置）：
  - 单帧底层（`utils/`）：界面层识别、H2O 拓扑、密度/取向/角度 PDF、
    cube 解析、slab-averaged potential
  - 水分析（`water/`）：界面→中点系综平均、吸附层识别/角度、三联图
  - 电势（`electrochemical/potential/`）：center potential、Fermi、
    electrode potential U vs SHE、φ(z) overlay、thickness sensitivity
  - 电荷（`electrochemical/charge/Bader/`）：表面电荷密度（counterion/
    layer 双方法）、原子电荷追踪、反离子逐帧检测；附 σ→φ 标定外推
  - 增强采样（`enhanced_sampling/`）：慢增长自由能（quick/publication）
    + 约束 TI 收敛诊断与自由能积分 + 恒电势修正
  - 脚本生成（`scripts/`）：Bader/TI/Potential/SP 工作目录批量生成
  - 双入口：CLI、`md_analysis.workflows`（`run_*`+`WorkflowResult`）
- **当前未覆盖**：按层/按元素电荷转移统计；Mulliken 电荷分析。

## 已确认的体系前提（用户声明，值得记录）

- 三基矢**正交**的周期性体系
- 始终包含**金属/水界面**
- 由于周期性边界条件，体系中总会存在**两个表面**（两个界面）

## 入口现状要点（只记策略与例外，不记逐项清单）

- **入口分层**：`cli` → `workflows` → 业务 → `engines` → `utils`。
  CLI 菜单命令全部走 `workflows.*` facade。
- **入口重构期间**移除了 legacy 名（旧 `run_*_analysis`/`run_all`），
  未保留 alias；业务经 `workflows.run_*` 调用。
- **关键坑**（影响后续 skill / 下游决策）：`WorkflowResult` **无**
  `.workdirs`；批量类 `run_*`（如 `run_ti_batch`）的工作目录路径在
  `WorkflowResult.artifacts`，强类型 report 在 `.extra`。
- **agent 路线例外**（`slowgrowth_quick` 直调、`config_show` read-only
  等）随 agent 层迁至 `development_agent` 分支，development 不再维护。

## 已批准行为变化（用户拍板，长效记录）

> 格式：default-preserving migration（参数默认保旧行为）/
> approved tightening（更合理但改旧行为，须用户批准）。

- **6.5 `strict`（default-preserving）**：`run_ti_full_analysis` 等
  workflow 入口 `strict: bool=True`（损坏点目录抛 `FileNotFoundError`）；
  CLI 显式传 `False` 维持菜单路径旧的宽松行为（WARN+skip）。
- **6.6 `overwrite`（default-preserving）**：`run_ti_batch` 的
  `overwrite: bool=False`（保留 collision guard），`True` 仅跳 guard
  逐文件覆盖**不清目录**；CLI 422 显式传 `True`。
- **6.5 TI dt 不一致（approved tightening，user 拍板接受）**：旧行为
  = warning + continue；新行为 = `ValueError`。理由：混用不同 dt 会让
  autocorrelation / N_eff / time-range 收敛诊断产生误导，继续算等于在
  不可信诊断上出结果，报错更合理。已在对应 commit message 标
  `BEHAVIOR CHANGE`。

## 近期任务清单（仍待补齐）

- **I/O**：统一记录单位/时间步/采样间隔等元数据（当前仅解析 `md.inp`
  的 `ABC [angstrom]`）
- **Analysis**：
  - 按层/按元素电荷转移统计（分层聚合 `bader_net_charge`，输出每层各
    元素平均净电荷）——待实现
  - Mulliken 电荷：按元素/分组/分层统计（优先级低于 Bader，需求待明确）
- **工程化**：固定最小依赖集合与安装方式说明
- **文档 follow-up**：中性化旧子阶段标签 `Phase 7b` / `Phase 7b1` /
  `Phase 7b2`（formats/engines 重构历史语境），避免与 canonical Phase 7
  混淆；独立计划处理。
- **旧 shim 弃用计划**：`batch_generate_ti_workdirs` 自入口重构 TI batch
  迁移后 CLI 无调用方（现行 backend = `generate_ti_batch_with_report`），
  但它仍是 **public Python API**（`scripts/__init__.py.__all__` + 直接
  单测）。**计划弃用**；真正删除须为**独立一笔 commit + 用户明确批准 +
  commit message 标 `API BREAK` + 同步 scripts public docs/tests**。
  本阶段（Phase 7）不删，仅登记此意图。

## 关键口径：已在当前实现中落地（不是待讨论）

- **界面参考面**：取"直接面向非金属环境"的金属界面层（每侧固定 1 层，共 2 层），使用其分数坐标的圆均值（`center_frac`）作为界面位置。界面层标签为 `"normal_aligned"` / `"normal_opposed"`。
- **法向/方向**：默认沿晶胞 `c` 轴（`normal="c"`，参考方向为 `+c_unit`）。自定义向量法向不支持。
- **水取向定义**：`theta` 定义为 H-O-H 角平分线与 `+c_unit` 的夹角；剖面使用 `cos(theta) * m_water` 做取向加权，单位为 `g/cm^3`。
- **取向加权密度单位**：`g/cm^3`（与质量密度同单位；公式为 $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$）。
- **CSV 列名**：取向列为 `orientation_ensemble_avg_g_cm3`（旧列名 `orientation_ensemble_avg_1_A3` 已废弃）。
- **电极电势**：`U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP`（cSHE 方法）。
