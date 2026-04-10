# md_analysis 包根 — 开发备忘

## 定位

包的顶层入口。管理 re-export、版本号、logging 初始化、用户配置持久化和编程入口点。

## 约定

- **Re-export 策略**：`__init__.py` 导出 `utils`、`water`、`electrochemical`，以及 `potential`/`charge`（从 electrochemical 提升）和 `MDAnalysisError`
- **不 re-export 的包**：`enhanced_sampling`、`scripts`、`agent` — 使用者需直接 `from md_analysis.agent import dispatch` 或 `from md_analysis.enhanced_sampling.slowgrowth import ...`
- **两个 config.py**：
  - `md_analysis/config.py` — 用户持久化配置（`~/.config/md_analysis/config.json`），管理 `KEY_VASP_SCRIPT_PATH` 等及电势输出参考配置（`KEY_POTENTIAL_REFERENCE`/`KEY_POTENTIAL_PH`/`KEY_POTENTIAL_TEMPERATURE_K`/`KEY_POTENTIAL_PHI_PZC`）
  - `md_analysis/utils/constants.py` — 物理常量和硬编码默认值（`HA_TO_EV`、`DEFAULT_LAYER_TOL_A` 等）
  - 混淆这两个是常见错误
- **NullHandler**：`__init__.py` 在 `md_analysis` logger 上设置 `NullHandler()`（PEP 282），CLI 或应用程序负责配置实际 handler
- **异常层次**：所有领域异常继承 `MDAnalysisError`（在 `exceptions.py` 定义），调用方可 `except MDAnalysisError` 统一捕获
- **编程入口**：`main.py` 提供 `run_water_analysis()`、`run_potential_analysis()`、`run_charge_analysis()`、`run_tracked_charge_analysis()`、`run_counterion_charge_analysis()`、`run_all()`
- **Agent 入口**：`agent/` 提供 `dispatch(task, params)` 统一调度、`get_task_schema()` 自动 JSON Schema、`list_tasks()` 任务枚举。薄适配层，不含分析逻辑。当前 12 个任务，含 `ti_full_analysis`（CLI 312，自定义 handler 编排 discover→load→analyze→plot+CSV）和 `sp_gen_batch`（CLI 442，批量生成 CP2K SP 目录用于 DeePMD 训练）
- **`run_*` 目录契约**：每个 `run_*_analysis()` 的 `output_dir` 参数即**最终写入目录**（不再自动前置 `water/`、`electrochemical/potential/` 等），仅保留必要的内部子目录（如 charge 的 `<method>/`）。调用方需自行提供完整路径；`run_all` 会按下述标准布局分派路径。
- **标准输出目录结构**（`run_all` 以及 CLI 菜单路径均按此推导）：
  - `<outdir>/water/`
  - `<outdir>/electrochemical/potential/<sub>/`（`<sub>`：`center`/`fermi`/`electrode`/`phi_z`/`thickness_sensitivity`）
  - `<outdir>/electrochemical/charge/<method>/`
  - `<outdir>/electrochemical/charge/tracked/` 或 `.../counterion_tracking/`
  - `<outdir>/electrochemical/calibration/{fit,predict}/`
  - `<outdir>/enhanced_sampling/slowgrowth/`
  - `<outdir>/enhanced_sampling/constrained_ti/`
- **导入规范**：包内用相对导入（`.`/`..`/`...`），测试用绝对导入

## 陷阱与历史 Bug

- `config.py` 的 `save_config()` 会自动创建父目录，但首次 `get_config()` 前配置文件可能不存在 — 返回 `None` 而非报错
- `main.py` 中的 import 全部延迟到函数体内（避免启动时加载 numpy/matplotlib）
- `CONFIGURABLE_DEFAULTS` 注册表的键必须与 `utils/constants.py` 中的默认常量一一对应

## 子目录

| 目录 | 用途 |
|---|---|
| `cli/` | 交互式 CLI → `cli/CLAUDE.md` |
| `agent/` | Agent-friendly 非交互式编程入口（dispatch + JSON Schema + TaskResult） |
| `utils/` | 底层工具 → `utils/CLAUDE.md` |
| `water/` | 水分析 → `water/CLAUDE.md` |
| `electrochemical/` | 电化学 → `electrochemical/CLAUDE.md` |
| `enhanced_sampling/` | 增强采样 → `enhanced_sampling/CLAUDE.md` |
| `scripts/` | 自动化脚本 → `scripts/CLAUDE.md` |
