# `md_analysis.cli` implementation guidelines

> Source: `src/md_analysis/cli/`

## Role

Interactive CLI package providing a VASPKIT-style numbered menu interface. Replaces the former `CLI.py` argparse module.

## Design principles

- Pure interactive: no command-line arguments, all input via `input()` prompts
- `_framework.py` 提供核心基础设施：`MenuNode`、`MenuGroup`、`MenuCommand`、`lazy_import()`
- `_params.py` 提供声明式参数采集：`ParamCollector` ABC + 泛型参数类（`StrParam`、`FloatParam`、`IntParam`、`ChoiceParam` 等）
- 每个子菜单模块通过 `MenuCommand` 子类实现，定义 `params`、`advanced_params`、`output_name` 和 `execute(self, ctx)` 方法
- `output_subdir` 由 `@property` 自动遍历父链（`parent.output_name`）拼接，无需硬编码
- `_prompt.py` 的共享 helper（`prompt_str`、`prompt_float` 等）由 `_params.py` 的参数类内部调用
- 错误处理由 `MenuCommand` 框架统一应用 `_handle_cmd_error` 装饰器
- `KeyboardInterrupt` / `EOFError` 在顶层 `main()` 中捕获以实现干净退出

## Dependencies

- `cli` -> `main` (for integrated workflow functions like `run_water_analysis`)
- `cli` -> `water`, `electrochemical.potential`, `electrochemical.charge`, `electrochemical.calibration` (for individual analysis functions)
- `cli` -> `electrochemical.potential.config` (for default constants)
- `cli` -> `scripts` (for `generate_bader_workdir`, `batch_generate_bader_workdirs`)
- `cli` -> `enhanced_sampling.slowgrowth` (via `lazy_import` for `slowgrowth_analysis`)
- `cli` -> `utils.CellParser` (for `parse_abc_from_restart`, `parse_abc_from_md_inp`)
- `cli` -> `utils.RestartParser.ColvarParser` (via `lazy_import` for `ColvarMDInfo`, used by `_enhanced_sampling.py` info display)
- `cli` -> `config` (for persistent user configuration, `CONFIGURABLE_DEFAULTS` registry, and `delete_config`)
- No reverse dependencies: no other module imports from `cli`

## Logging configuration

- `main()` configures a `StreamHandler` on the `"md_analysis"` logger at `INFO` level, outputting to stderr with format `"%(levelname)s: %(message)s"`
- Configuration is idempotent: only adds the handler if no non-NullHandler handlers exist
- The library itself uses `NullHandler` (set in `md_analysis/__init__.py`), so logging is silent unless the CLI (or an application) explicitly configures a handler
- `_handle_cmd_error` logs unexpected exceptions at `ERROR` level with `exc_info=True` for full traceback in logs, while printing a concise message to stdout for the user

## 命令类架构

所有子菜单模块使用 `MenuCommand` 子类实现（取代了旧的 `_cmd_<code>()` 函数模式）：

| 模块 | 命令类 | 菜单码 |
|------|--------|--------|
| `_water.py` | `WaterDensityCmd`, `WaterOrientationCmd`, `AdWaterOrientationCmd`, `AdWaterThetaCmd`, `WaterThreePanelCmd` | 101-105 |
| `_potential.py` | `CenterPotentialCmd`, `FermiEnergyCmd`, `ElectrodePotentialCmd`, `PhiZProfileCmd`, `ThicknessSensitivityCmd`, `FullPotentialCmd` | 211-216 |
| `_charge.py` | `SurfaceChargeCmd`（通过 `method` 参数区分 counterion/layer/prompted） | 221-223 |
| `_enhanced_sampling.py` | `SGQuickPlotCmd`, `SGPublicationPlotCmd`（共享基类 `_SlowgrowthPlotCmd`） | 301-302 |
| `_scripts.py` | `BaderSingleCmd`, `BaderBatchCmd`, `TISingleCmd`, `TIBatchCmd` | 411-412 (sub-group 41), 421-422 (sub-group 42) |
| `_settings.py` | `SetVaspScriptCmd`, `ShowConfigCmd`, `SetAnalysisDefaultCmd`（通过 `config_key` 参数复用）, `ResetDefaultsCmd` | 901-907 |

`_charge.py` 中 `_print_ensemble_summary()` 作为独立辅助函数保留，在 `SurfaceChargeCmd.execute()` 结束时调用。

## Cell parameter acquisition

All sub-menus requiring cell parameters (water 101-105, scripts 411-412) use the shared `_prompt_cell_abc()` helper from `_prompt.py`:

1. Prompt cell source: `.restart` (default) or `md.inp`
2. Parse the chosen file (`parse_abc_from_restart` or `parse_abc_from_md_inp`)
3. On failure, offer one retry with a different file
4. On second failure, raise `CellParseError` (caught by `@_handle_cmd_error`)
5. Return `(a, b, c)` tuple, passed as `cell_abc` keyword argument to analysis functions

## Parameter flow

1. User selects analysis code from sub-menu
2. Required parameters prompted (if any)
3. "Modify advanced parameters? (y/N)" gate for optional parameters
4. Handler calls the appropriate analysis function
5. Results printed, program exits

## Configurable analysis defaults

Settings menu 903-907 allow users to persistently override algorithm defaults from `utils/config.py`:

- 903: `layer_tol_A` (layer clustering tolerance)
- 904: `z_bin_width_A` (z-axis bin width)
- 905: `theta_bin_deg` (theta bin width)
- 906: `water_oh_cutoff_A` (water O-H cutoff)
- 907: reset all analysis defaults

The `_get_effective_default(key)` helper in `_prompt.py` reads the user config first, falling back to the hardcoded default from `CONFIGURABLE_DEFAULTS` registry. Analysis sub-menus (`_potential.py`, `_water.py`, `_charge.py`) use this helper to populate prompt defaults. Library function signatures remain unchanged — persistence only affects CLI prompt defaults.

## Enhanced sampling CLI (`_enhanced_sampling.py`)

`_enhanced_sampling.py` 通过 `MenuCommand` 子类（`SGQuickPlotCmd`、`SGPublicationPlotCmd`）实现 301/302 菜单项，共享基类 `_SlowgrowthPlotCmd`：

- **文件发现**：
  - `_discover_restart_file(workdir)`：glob `*.restart`，排除 `_\d+\.restart` 检查点文件
  - `_discover_log_file(workdir)`：glob `*.LagrangeMultLog`
  - 两者均要求恰好 1 个匹配；否则回退到用户手动输入
- **轨迹信息预览**：`_print_sg_info()` 通过 `lazy_import` 获取 `ColvarMDInfo`，显示步数、时间步、CV 范围，并检测 NaN（溢出）步
- **参数采集**：restart path、log path、initial/final step、colvar ID、output dir
- **执行**：通过 `lazy_import` 调用 `slowgrowth_analysis`，`_plot_style` 由子类决定（`"quick"` 或 `"publication"`）
