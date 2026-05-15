# `md_analysis.cli` implementation guidelines

> Source: `src/md_analysis/cli/`

## Role

Interactive CLI package providing a VASPKIT-style numbered menu interface. Replaces the former `CLI.py` argparse module.

## Design principles

- Pure interactive: no command-line arguments, all input via `input()` prompts
- `_framework.py` 提供核心基础设施：`MenuNode`、`MenuGroup`、`MenuCommand`、`lazy_import()`
- `_params.py` 提供声明式参数采集：`ParamCollector` ABC + 泛型参数类（`StrParam`、`FloatParam`、`IntParam`、`ChoiceParam`、`ConditionalParam`、`BoolParam`、`DisplayAction`、`ConfigStrParam` 等）
- 每个子菜单模块通过 `MenuCommand` 子类实现，定义 `params`、`advanced_params`、`output_name` 和 `execute(self, ctx)` 方法
- `output_subdir` 由 `@property` 自动遍历父链（`parent.output_name`）拼接，无需硬编码
- `_prompt.py` 的共享 helper（`prompt_str`、`prompt_float` 等）由 `_params.py` 的参数类内部调用
- 错误处理由 `MenuCommand.run()` 的 inline try-except 统一处理
- `KeyboardInterrupt` / `EOFError` 在顶层 `main()` 中捕获以实现干净退出

## Dependencies

- `cli` -> `workflows` (for integrated workflow facades like `run_water_three_panel` and `run_potential_full`; the entrance refactor removed the CLI dependency on the legacy `md_analysis.main:run_*_analysis` shims)
- `cli` -> `water`, `electrochemical.potential`, `electrochemical.charge`, `electrochemical.calibration` (for individual analysis functions)
- `cli` -> `electrochemical.potential.config` (for default constants)
- `cli` -> `scripts` (for `generate_bader_workdir`, `batch_generate_bader_workdirs`)
- `cli` -> `enhanced_sampling.slowgrowth` (via `lazy_import` for `slowgrowth_analysis`)
- `cli` -> `enhanced_sampling.constrained_ti` (via `lazy_import` for `standalone_diagnostics`, `analyze_ti`, `discover_ti_points`, etc.)
- `cli` -> `utils.formats.cp2k_cell` (for `parse_abc_from_restart`, `parse_abc_from_md_inp`)
- `cli` -> `utils.formats.cp2k_colvar` (via `lazy_import` for `ColvarMDInfo`, used by `_enhanced_sampling.py` info display)
- `cli` -> `config` (for persistent user configuration, `CONFIGURABLE_DEFAULTS` registry, and `delete_config`)
- No reverse dependencies: no other module imports from `cli`

## Logging configuration

- `main()` configures a `StreamHandler` on the `"md_analysis"` logger at `INFO` level, outputting to stderr with format `"%(levelname)s: %(message)s"`
- Configuration is idempotent: only adds the handler if no non-NullHandler handlers exist
- The library itself uses `NullHandler` (set in `md_analysis/__init__.py`), so logging is silent unless the CLI (or an application) explicitly configures a handler
- `MenuCommand.run()` 的 except 分支对未知异常以 `ERROR` level 加 `exc_info=True` 记录完整 traceback，同时在 stdout 打印简洁消息给用户

## 命令类架构

所有子菜单模块使用 `MenuCommand` 子类实现（取代了旧的 `_cmd_<code>()` 函数模式）：

| 模块 | 命令类 | 菜单码 |
|------|--------|--------|
| `_water.py` | `WaterDensityCmd`, `WaterOrientationCmd`, `AdWaterOrientationCmd`, `AdWaterThetaCmd`, `WaterThreePanelCmd` | 101-105 |
| `_potential.py` | `CenterPotentialCmd`, `FermiEnergyCmd`, `ElectrodePotentialCmd`, `PhiZProfileCmd`, `ThicknessSensitivityCmd`, `FullPotentialCmd` (all support `input_mode`: continuous/distributed) | 211-216 |
| `_charge.py` | `SurfaceChargeCmd`（通过 `method` 参数区分 counterion/layer/prompted）, `SingleSideChargeCmd`, `TrackedChargeCmd`, `CounterionChargeCmd` | 221-226 |
| `_calibration.py` | `CalibrateFromCSVCmd`, `CalibrateManualCmd`, `PredictPotentialCmd` | 231-233 |
| `_enhanced_sampling.py` | `SGQuickPlotCmd`, `SGPublicationPlotCmd`（共享基类 `_SlowgrowthPlotCmd`） | 301-302 (sub-group 30) |
| `_constrained_ti.py` | `TISingleDiagCmd`, `TIFullAnalysisCmd`, `TIConstPotCorrectionCmd` | 311-313 (sub-group 31) |
| `_scripts.py` | `BaderSingleCmd`, `BaderBatchCmd`, `TISingleCmd`, `TIBatchCmd`, `PotentialSingleCmd`, `PotentialBatchCmd` | 411-412 (sub-group 41), 421-422 (sub-group 42), 431-432 (sub-group 43) |
| `_settings.py` | `ShowConfigCmd`, `ResetDefaultsCmd` (sub-group 90); `SetVaspScriptCmd`, `SetCp2kScriptCmd`, `SetSpInpTemplateCmd` (sub-group 91); `SetAnalysisDefaultCmd`×4 (sub-group 92); `SetPotentialReferenceCmd` (sub-group 93) | 900-931 |

`SurfaceChargeCmd.execute()` 使用 `surface_charge_analysis()` 返回的 `SurfaceChargeResult` 数据类直接打印系综统计摘要（不再依赖独立的 `_print_ensemble_summary()` 函数，已删除）。

221/222（固定 method）通过 `output_name` 由框架自动解析输出子目录（`charge/counterion/`、`charge/layer/`）。223（动态 method）`output_name` 为空，`execute()` 中手动追加 `ctx[K.METHOD]` 到输出路径以避免不同方法输出覆盖同一文件。

## Cell parameter acquisition

All sub-menus requiring cell parameters (water 101-105, scripts 411-412) use the shared `_prompt_cell_abc()` helper from `_prompt.py`:

1. Prompt cell source: `.restart` (default) or `md.inp`
2. Parse the chosen file (`parse_abc_from_restart` or `parse_abc_from_md_inp`)
3. On failure, offer one retry with a different file
4. On second failure, raise `CellParseError` (caught by `MenuCommand.run()` 的 try-except)
5. Return `(a, b, c)` tuple, passed as `cell_abc` keyword argument to analysis functions

## Parameter flow

1. User selects analysis code from sub-menu
2. Required parameters prompted (if any)
3. "Modify advanced parameters? (y/N)" gate for optional parameters
4. Handler calls the appropriate analysis function
5. Results printed, control returns to parent menu

## Configurable analysis defaults

Settings sub-group 92 allows users to persistently override algorithm defaults from `utils/constants.py`:

- 921: `layer_tol_A` (layer clustering tolerance)
- 922: `z_bin_width_A` (z-axis bin width)
- 923: `theta_bin_deg` (theta bin width)
- 924: `water_oh_cutoff_A` (water O-H cutoff)
- 909: reset all analysis defaults (sub-group 90)
- 931: potential output reference (SHE/RHE/PZC) — `SetPotentialReferenceCmd` (sub-group 93)

`_ConfigBackedParam` 内部基类（在 `_params.py` 中定义）提供 `_get_config_value()` 和基础 `apply_default()` 实现。两个子类：
- `ConfigDefaultParam`：通过 `CONFIGURABLE_DEFAULTS` 注册表 fallback 到硬编码默认值。Analysis sub-menus (`_potential.py`, `_water.py`, `_charge.py`) 使用此参数类来填充提示默认值。
- `ConfigStrParam`：直接通过 `get_config(key)` 读取字符串配置值（不经过 `CONFIGURABLE_DEFAULTS` 注册表）。预定义实例 `vasp_script`、`cp2k_script`、`sp_inp_template` 供 Scripts 命令使用。

`BoolParam`：yes/no 布尔提示参数（字段：`key`、`label`、`default`）。预定义实例 `gen_potcar` 供 Bader Scripts 命令使用。

`DisplayAction`：在参数采集过程中执行副作用（如打印轨迹信息），不在 ctx 中存储值。`apply_default` 为 no-op。Scripts 命令（411-412、431-432）使用 `DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ]))` 在采集完 XYZ 路径后立即显示轨迹元信息。

Library function signatures remain unchanged — persistence only affects CLI prompt defaults.

## Conditional parameter collection (`ConditionalParam`)

`ConditionalParam(inner, predicate)` wraps another `ParamCollector` and prompts only when `predicate(ctx)` returns true. Otherwise `inner.apply_default(ctx)` is called silently. Used by Potential commands (211-216) to skip distributed-mode parameters (`sp_root_dir`, `sp_dir_pattern`, `sp_cube_filename`, `sp_out_filename`) when the user selects `input_mode = "continuous"`. The predicate is evaluated against the live ctx, so any keys it reads must be collected earlier in the `params` tuple.

### Potential output reference (909)

`SetPotentialReferenceCmd` allows users to configure the default output potential reference scale for `surface_charge_analysis` extrapolation:

- **SHE** (default): no conversion
- **RHE**: prompts for pH and temperature (K), applies Nernst shift φ_RHE = φ_SHE + (RT/F)·ln(10)·pH
- **PZC**: prompts for φ_PZC (V vs SHE), applies φ_PZC = φ_SHE − φ_pzc

Config keys: `KEY_POTENTIAL_REFERENCE`, `KEY_POTENTIAL_PH`, `KEY_POTENTIAL_TEMPERATURE_K`, `KEY_POTENTIAL_PHI_PZC`. `ShowConfigCmd` displays these in a separate "Potential Output" section; `ResetDefaultsCmd` also clears them. `SurfaceChargeCmd.execute()` reads these from config and passes to `surface_charge_analysis()`.

## Enhanced sampling CLI (`_enhanced_sampling.py` + `_constrained_ti.py`)

菜单 3 下分两个子组：`30)` Slow-Growth、`31)` Constrained TI Analysis。

两个模块均遵循延迟导入约定：顶层无 `import numpy` 或 `from ..utils.constants import ...`，numpy 和常量均在使用它们的函数体内导入。

### Slow-Growth（`_enhanced_sampling.py`，sub-group 30）

301/302 共享基类 `_SlowgrowthPlotCmd`（注意：`output_name` 由父 `MenuGroup("30", output_name="slowgrowth")` 提供，命令自身无 `output_name`）：

- **文件发现**：
  - `_discover_restart_file(workdir)`：glob `*.restart`，排除 `_\d+\.restart` 检查点文件
  - `_discover_log_file(workdir)`：glob `*.LagrangeMultLog`
  - 两者均要求恰好 1 个匹配；否则回退到用户手动输入
- **轨迹信息预览**：`_print_sg_info()` 通过 `lazy_import` 获取 `ColvarMDInfo`，显示步数、时间步、CV 范围，并检测 NaN（溢出）步。numpy 和 `AU_TIME_TO_FS` 常量在此函数体内导入
- **参数采集**：restart path、log path、initial/final step、colvar ID、output dir
- **执行**：通过 `lazy_import` 调用 `slowgrowth_analysis`，`_plot_style` 由子类决定（`"quick"` 或 `"publication"`）

### Constrained TI（`_constrained_ti.py`，sub-group 31）

模块级共享辅助函数（312/313 去重）：
- **`_collect_ti_base_params(ctx)`**：采集 TI_ROOT_DIR、TI_DIR_PATTERN（含 `_VALID_PATTERNS` 校验循环）、EQUILIBRATION、EPSILON_TOL_EV、TI_REVERSE。不采集 OUTDIR
- **`_run_ti_core(ctx)`**：共享 TI 分析流水线：discover → 交互式逐点 equilibration → load → dt 一致性校验 → analyze → 终端摘要 → 写出文件。返回 `(ti_report, point_defs, xi_values, outdir)`。numpy 在此函数体内导入
- **`_VALID_PATTERNS`**、**`_VALID_SIDES`**、**`_VALID_METHODS`** 为模块级常量

命令类：
- **311 `TISingleDiagCmd`**：单点收敛诊断。复用 SG 的 `_discover_restart_file` / `_discover_log_file`。覆写 `_collect_all_params()`（SEM target 为可空 float，用 `prompt_str` + 手动转换）。调用 `standalone_diagnostics()`。
- **312 `TIFullAnalysisCmd`**：调用 `_collect_ti_base_params()` + `_run_ti_core()`，打印 ΔA + 文件列表
- **313 `TIConstPotCorrectionCmd`**：调用 `_collect_ti_base_params()` + `_run_ti_core()`（Phase 1），再执行 Phase 2 恒电势修正（Nørskov）。从各 `ti_target_*/bader/` 提取系综平均 σ，通过 calibration mapper 外推 Φ，计算修正项并输出修正后自由能曲线。需要 calibration.json（硬错误）；缺 bader/ 时 WARN 并跳过修正。额外参数：`K.TARGET_SIDE`（aligned/opposed）、`K.CALIBRATION_JSON`。
