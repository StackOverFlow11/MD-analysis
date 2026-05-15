# md_analysis 包根 — 开发备忘

## 定位

包的顶层入口。管理 re-export、版本号、logging 初始化、用户配置持久化和编程入口点。

## 约定

- **Re-export 策略**：`__init__.py` 导出 `utils`、`water`、`electrochemical`，以及 `potential`/`charge`（从 electrochemical 提升）和 `MDAnalysisError`
- **不 re-export 的包**：`enhanced_sampling`、`scripts` — 使用者需直接 `from md_analysis.enhanced_sampling.slowgrowth import ...`
- **两个 config.py**：
  - `md_analysis/config.py` — 用户持久化配置（`~/.config/md_analysis/config.json`），管理 `KEY_VASP_SCRIPT_PATH` 等及电势输出参考配置（`KEY_POTENTIAL_REFERENCE`/`KEY_POTENTIAL_PH`/`KEY_POTENTIAL_TEMPERATURE_K`/`KEY_POTENTIAL_PHI_PZC`）
  - `md_analysis/utils/constants.py` — 物理常量和硬编码默认值（`HA_TO_EV`、`DEFAULT_LAYER_TOL_A` 等）
  - 混淆这两个是常见错误
- **NullHandler**：`__init__.py` 在 `md_analysis` logger 上设置 `NullHandler()`（PEP 282），CLI 或应用程序负责配置实际 handler
- **异常层次**：所有领域异常继承 `MDAnalysisError`（在 `exceptions.py` 定义），调用方可 `except MDAnalysisError` 统一捕获
- **编程入口**：`workflows/` 提供 `run_*()` 函数集，全部返回 `WorkflowResult`（`artifacts` / `metadata` / `extra`）。`main.py` 是同一批名字的薄 re-export facade（import-only）。**完整清单以 `workflows/__init__.py.__all__` 为权威**（亦见 `workflows/CLAUDE.md`）；composite `run_interface_analysis` 取代旧 `run_all`。
- **`run_*` 目录契约**：每个 `run_*()` 的 `output_dir` 参数即**最终写入目录**（不再自动前置 `water/`、`electrochemical/potential/` 等），仅保留必要的内部子目录（如 surface charge 的 `<method>/`）。调用方需自行提供完整路径；`run_interface_analysis` 会按下述标准布局分派路径。
- **标准输出目录结构**（`run_interface_analysis` 以及 CLI 菜单路径均按此推导）：
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
- `workflows/*.py` 把重业务模块（matplotlib、ase、底层分析）`import` 延迟到函数体内（避免启动时加载 numpy/matplotlib）；`main.py` 是 import-only facade，模块顶层只 import `workflows` 包内的 `run_*` 名字，本身不引入业务依赖
- `CONFIGURABLE_DEFAULTS` 注册表的键必须与 `utils/constants.py` 中的默认常量一一对应

## 子目录

| 目录 | 用途 |
|---|---|
| `cli/` | 交互式 CLI → `cli/CLAUDE.md` |
| `workflows/` | 程序化入口 facade（`run_*` + `WorkflowResult`；清单见 `workflows/__init__.py.__all__`） |
| `engines/` | CP2K/VASP 引擎门面 + engine-neutral dataclass → `engines/CLAUDE.md` |
| `utils/` | 底层工具（formats / structure / io 三层）→ `utils/CLAUDE.md` |
| `water/` | 水分析 → `water/CLAUDE.md` |
| `electrochemical/` | 电化学 → `electrochemical/CLAUDE.md` |
| `enhanced_sampling/` | 增强采样 → `enhanced_sampling/CLAUDE.md` |
| `scripts/` | 自动化脚本 → `scripts/CLAUDE.md` |

## 依赖方向（utils/engines 重构后）

```
cli / scripts
   ↓
main.py / workflows (upper)
   ↓
water / electrochemical / enhanced_sampling   ← 业务工作流
   ↓
engines                                       ← CP2K / VASP 门面 + engine-neutral models
   ↓
utils/{formats, structure, io}, constants     ← 单文件解析 + 几何 helper + 路径发现
   ↓
exceptions
```

约束：
- `engines/` 只允许 import `utils/`、`exceptions` 和自身内部；**不允许**反向 import `electrochemical` / `water` / `enhanced_sampling` / `cli` / `scripts`
- `utils/` 运行时**不**依赖 `engines/`（只允许 `TYPE_CHECKING` 块下的字符串注解）
- 业务工作流 import `engines.cp2k.read_*` facade，而不是直接走 `utils.formats.cp2k.*`
