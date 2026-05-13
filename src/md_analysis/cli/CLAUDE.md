# cli — 开发备忘

## 定位

VASPKIT 风格交互式编号菜单。无 argparse，所有输入通过 `input()` 提示。

## 约定

### 菜单编号方案
- `1xx`：Water (101-105)
- `21x`：Potential (211-216), `22x`：Charge (221-226), `23x`：Calibration (231-233)
- `30x`：Slow-Growth (301-302), `31x`：Constrained TI (311-313)
- `41x`：Bader (411-412), `42x`：TI (421-422), `43x`：SP Potential (431-432)
- `90x`：Show/Reset (900, 909), `91x`：Script Paths (911-913), `92x`：Analysis Defaults (921-924), `93x`：Potential Output (931)
- 编号前缀必须与父 MenuGroup 编号匹配

### 框架模式
- `MenuGroup`（非叶）/ `MenuCommand`（叶），在 `__init__.py` 的 `build_menu_tree()` 中组装
- `build_flat_index()` 使任意叶节点编号可从根直达（如直接输入 `421`）
- **flat index 只索引 `MenuCommand`**，不索引 `MenuGroup` — 输入子菜单编号（如 `42`）必须先进入父菜单

### 输出目录自动推导（output_name 机制）
- 每个 `MenuNode` 拥有 `output_name: str` 属性，贡献一段路径
- `MenuGroup` 通过构造参数 `output_name=` 设置（如 `MenuGroup("2", ..., output_name="electrochemical")`）
- `MenuCommand` 通过类属性或 `__init__` 中设置 `self.output_name`
- `MenuCommand.output_subdir` 是 `@property`，自动遍历父链拼接所有 `output_name` 段
- `MenuGroup.add()` 自动建立 `node.parent` 引用
- 移动命令到不同 `MenuGroup` 时输出路径自动更新
- 无 `output_name` 的节点（Scripts、Settings）不参与路径推导

### lazy_import
- 所有 `execute()` 方法通过 `lazy_import()` 延迟加载分析模块（36+ 处使用）
- 目的：CLI 启动只加载框架代码，不触发 numpy/matplotlib/ase
- **numpy 也必须延迟导入**：CLI 模块顶层不得出现 `import numpy`，需要时在函数体内导入
- 常量（如 `AU_TIME_TO_FS`）同理：从 `utils.constants` 导入时放在使用它的函数体内

### Workflows facade 迁移状态（入口重构）

入口重构后，`md_analysis.workflows.*` 是 canonical 程序化入口。CLI 的 `execute()`
应优先通过 `lazy_import("md_analysis.workflows.<domain>", "run_*")` 调度，
读取 `WorkflowResult.artifacts` 拿文件路径，并从 `metadata` / `extra` 取摘要指标；
不应再回到底层科学模块拼装流程。

**已迁到 workflows facade**：

| CLI | 命令 | workflow facade |
|---|---|---|
| 105 | `WaterThreePanelCmd` | `workflows.water.run_water_three_panel` |
| 216 | `FullPotentialCmd` | `workflows.potential.run_potential_full` |
| 225 | `TrackedChargeCmd` | `workflows.charge.run_tracked_charge` |
| 226 | `CounterionChargeCmd` | `workflows.charge.run_counterion_charge` |
| 301 | `SGQuickPlotCmd` | `workflows.enhanced_sampling.run_slowgrowth_quick_plot` |
| 302 | `SGPublicationPlotCmd` | `workflows.enhanced_sampling.run_slowgrowth_publication_plot` |
| 411 | `BaderSingleCmd` | `workflows.scripts.run_bader_single` |
| 412 | `BaderBatchCmd` | `workflows.scripts.run_bader_batch` |
| 421 | `TISingleCmd` | `workflows.scripts.run_ti_single` |
| 431 | `PotentialSingleCmd` | `workflows.scripts.run_potential_single` |
| 432 | `PotentialBatchCmd` | `workflows.scripts.run_potential_batch` |
| 441 | `SpGenSingleCmd` | `workflows.scripts.run_sp_single` |
| 442 | `SpGenBatchCmd` | `workflows.scripts.run_sp_batch` |

**保留底层直调（workflows facade 当前覆盖不全；都已记录原因）**：

| CLI | 命令 | 保留原因 |
|---|---|---|
| 101 / 102 / 103 / 104 | water 单步命令 | `workflows.water` 当前只导出 `run_water_three_panel`（composite），没有 density / orientation / adsorbed / theta 的单步 `run_*` |
| 211 / 212 / 213 / 214 / 215 | potential 单步命令 | `workflows.potential` 当前只导出 `run_potential_full`（composite），没有 center / fermi / electrode / phi_z / thickness_sensitivity 的单步 `run_*` |
| 221 / 222 / 223 | `SurfaceChargeCmd` | `workflows.charge.run_surface_charge` 不接受 `potential_reference` / `potential_pH` / `potential_temperature_K` / `potential_phi_pzc` 这 4 个底层 `surface_charge_analysis` 已支持的参数；CLI 当前需要它们才能输出 φ 列 |
| 224 | `SingleSideChargeCmd` | `workflows.charge.run_surface_charge` 不接受 `target_side`；底层 `surface_charge_analysis` 支持 |
| 231 / 232 / 233 | calibration | `workflows.calibration.run_calibration_fit` / `run_calibration_predict` 要求 `calibration_json_path` 必填，CLI 当前允许传 None 走全局默认（`~/.config/md_analysis/calibration.json`）。让 CLI 显式提前 resolve 默认值不在 Phase 4 范围内 |
| 311 / 312 / 313 | constrained TI | CLI 有 Python-slice 切片 UI、逐点 equilibration override、约束点交互列表等富交互流程；workflow `run_ti_full_analysis` / `run_ti_constant_potential_correction` 是单次端到端调用，目前不能完整覆盖这套交互细节 |
| 422 | `TIBatchCmd` | `workflows.scripts.run_ti_batch`（底层 `generate_ti_batch_with_report`）不接受 `colvar_id`（MVP 限定 primary CV）；CLI 仍 expose colvar_id |

任何后续往 workflows 收敛的工作，前提是先扩展 workflow 签名覆盖这些 gap；
不要在 CLI 端 hack 绕过。当 workflow 覆盖完整后再统一迁这些命令。

### 参数采集
- `K` 类：字符串键常量，防止拼写错误（含 `INP_TEMPLATE`、`GEN_POTCAR` 等）
- `ParamCollector` ABC：`collect(ctx)` 提示用户 + `apply_default(ctx)` 静默填充
- `params` 元组：总是提示；`advanced_params` 元组：用户选择"修改高级参数"时才提示
- **`_ConfigBackedParam(ParamCollector)`**：内部基类，提供 `_get_config_value()` 方法和基础 `apply_default()` 实现，统一配置读取逻辑
- **`ConfigDefaultParam(_ConfigBackedParam)`**：从 `~/.config/md_analysis/config.json` 读取用户覆盖值，fallback 到 `CONFIGURABLE_DEFAULTS` 注册表中的硬编码默认值
- **`ConfigStrParam(_ConfigBackedParam)`**：字符串提示参数，带 config-backed 默认值。直接通过 `get_config(key)` 读取（不经过 `CONFIGURABLE_DEFAULTS` 注册表）。预定义实例：`vasp_script`、`cp2k_script`、`sp_inp_template`
- **`BoolParam(ParamCollector)`**：yes/no 提示。字段：`key`、`label`、`default`。预定义实例：`gen_potcar`
- **`DisplayAction(ParamCollector)`**：在参数采集过程中执行副作用（如打印轨迹信息），不在 ctx 中存储值。`apply_default` 为 no-op
- 所有 Potential 命令（211-216）的 `params` 元组首位为 `input_mode`（`ChoiceParam`："continuous"/"distributed"），后接模式相关参数（`sp_root_dir`、`sp_dir_pattern`、`sp_cube_filename`、`sp_out_filename`）。这 4 个 sp_* 参数通过 `ConditionalParam` 包装，仅当 `ctx[K.INPUT_MODE] == "distributed"` 时提示用户，否则静默应用默认值。`execute()` 通过 `_is_distributed(ctx)` 分派调用
- **`ConditionalParam(inner, predicate)`**：通用包装器，仅当 `predicate(ctx)` 为真时调用 `inner.collect(ctx)`，否则调用 `inner.apply_default(ctx)`。要求 predicate 依赖的 ctx 键在 params 元组中先于此参数出现
- Scripts 命令（411-412、431-432）已转换为声明式 `params` 元组（使用 `ConfigStrParam`、`BoolParam`、`DisplayAction`），不再覆写 `_collect_all_params`。TI 命令（421-422）仍使用 `_resolve_cp2k_script_path()` + `_collect_all_params` 覆写模式

### 错误处理
- `MenuCommand.run()` 的 inline try-except 捕获 `MDAnalysisError`/`FileNotFoundError`/`ValueError`/`RuntimeError` → 打印简洁消息
- 未知异常 → `logger.error(..., exc_info=True)` 记录完整 traceback + 打印简洁消息到 stdout

### 测试钩子
- `_prompt.py` 的 `set_input_source(fn)` 可注入自定义输入函数

## 新增命令检查清单

1. 在 `_<module>.py` 中创建 `MenuCommand` 子类，定义 `params`/`advanced_params`/`output_name`/`execute()`
2. 在 `__init__.py` 中 import 并在 `build_menu_tree()` 中注册到对应 `MenuGroup`
3. 如需新参数键 → 在 `_params.py` 的 `K` 类中添加常量
4. 如需新参数类型 → 创建 `ParamCollector` 子类或使用现有泛型类

## Constrained TI (312/313) 交互流程

`_collect_ti_base_params` 采集共享参数，包含：
- `K.TI_ROOT_DIR`、`K.EQUILIBRATION`、`K.EPSILON_TOL_EV`、`K.TI_REVERSE`
- `K.AUTO_EQUILIBRATION`：可选自动预平衡迭代（二分砍前半直到收敛）

CLI 不再 prompt directory pattern：discover_ti_points 默认 `parser="auto" + dir_filter=None`（嗅探 + 内容过滤），目录命名完全自由。

`_run_ti_core` 执行共享 TI 分析，流程：
1. 发现约束点 → 带索引列表显示 `[0] ξ = ...`
2. **Python 切片选择**（可选）：用户输入如 `3:8`、`::2`、`:8` 等，空回车 = 全部
3. 可选逐点 equilibration 覆盖
4. 加载数据 → dt 一致性检查
5. `analyze_ti(... auto_equilibration=ctx[K.AUTO_EQUILIBRATION])` → 控制台摘要表 → 写文件

## 陷阱与历史 Bug

- 菜单码重编号（bace527）：旧代码中 401/402 已改为 411/412
- `CellAbcParam.collect()` 允许一次重试（.restart 失败 → 切 md.inp），第二次失败才报错
- `_discover_restart_file()` 排除 `_\d+.restart` 检查点文件（正则过滤）
- SG 命令会检测 LagrangeMultLog 中的 overflow（NaN 步），并在终端打印警告
- SG 命令 301/302 的 `output_name` 由父 `MenuGroup("30", output_name="slowgrowth")` 提供，`_SlowgrowthPlotCmd` 自身不定义 `output_name`（否则路径重复拼接为 `slowgrowth/slowgrowth`）
- `K.TI_DIR_PATTERN` 已删除（2026-05-10 IO 重构）；TI discover 现走 parser-driven 自动模式。`K.DIR_PATTERN`（Bader 用）保留
