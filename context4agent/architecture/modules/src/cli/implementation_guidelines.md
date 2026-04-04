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
- 错误处理由 `MenuCommand.run()` 的 inline try-except 统一处理
- `KeyboardInterrupt` / `EOFError` 在顶层 `main()` 中捕获以实现干净退出

## Dependencies

- `cli` -> `main` (for integrated workflow functions like `run_water_analysis`)
- `cli` -> `water`, `electrochemical.potential`, `electrochemical.charge`, `electrochemical.calibration` (for individual analysis functions)
- `cli` -> `electrochemical.potential.config` (for default constants)
- `cli` -> `scripts` (for `generate_bader_workdir`, `batch_generate_bader_workdirs`)
- `cli` -> `enhanced_sampling.slowgrowth` (via `lazy_import` for `slowgrowth_analysis`)
- `cli` -> `enhanced_sampling.constrained_ti` (via `lazy_import` for `standalone_diagnostics`, `analyze_ti`, `discover_ti_points`, etc.)
- `cli` -> `utils.CellParser` (for `parse_abc_from_restart`, `parse_abc_from_md_inp`)
- `cli` -> `utils.RestartParser.ColvarParser` (via `lazy_import` for `ColvarMDInfo`, used by `_enhanced_sampling.py` info display)
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
| `_settings.py` | `SetVaspScriptCmd`, `SetCp2kScriptCmd`, `ShowConfigCmd`, `SetAnalysisDefaultCmd`（通过 `config_key` 参数复用）, `ResetDefaultsCmd`, `SetPotentialReferenceCmd`, `SetSpInpTemplateCmd` | 901-910 |

`_charge.py` 中 `_print_ensemble_summary()` 作为独立辅助函数保留，在 `SurfaceChargeCmd.execute()` 结束时调用。

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

Settings menu 903-907 allow users to persistently override algorithm defaults from `utils/config.py`:

- 903: `layer_tol_A` (layer clustering tolerance)
- 904: `z_bin_width_A` (z-axis bin width)
- 905: `theta_bin_deg` (theta bin width)
- 906: `water_oh_cutoff_A` (water O-H cutoff)
- 907: reset all analysis defaults
- 909: potential output reference (SHE/RHE/PZC) — `SetPotentialReferenceCmd`

`ConfigDefaultParam` 类（在 `_params.py` 中定义）通过 `apply_default()` 方法读取用户持久化配置，fallback 到 `CONFIGURABLE_DEFAULTS` 注册表中的硬编码默认值。Analysis sub-menus (`_potential.py`, `_water.py`, `_charge.py`) 使用此参数类来填充提示默认值。Library function signatures remain unchanged — persistence only affects CLI prompt defaults.

### Potential output reference (909)

`SetPotentialReferenceCmd` allows users to configure the default output potential reference scale for `surface_charge_analysis` extrapolation:

- **SHE** (default): no conversion
- **RHE**: prompts for pH and temperature (K), applies Nernst shift φ_RHE = φ_SHE + (RT/F)·ln(10)·pH
- **PZC**: prompts for φ_PZC (V vs SHE), applies φ_PZC = φ_SHE − φ_pzc

Config keys: `KEY_POTENTIAL_REFERENCE`, `KEY_POTENTIAL_PH`, `KEY_POTENTIAL_TEMPERATURE_K`, `KEY_POTENTIAL_PHI_PZC`. `ShowConfigCmd` displays these in a separate "Potential Output" section; `ResetDefaultsCmd` also clears them. `SurfaceChargeCmd.execute()` reads these from config and passes to `surface_charge_analysis()`.

## Enhanced sampling CLI (`_enhanced_sampling.py` + `_constrained_ti.py`)

菜单 3 下分两个子组：`30)` Slow-Growth、`31)` Constrained TI Analysis。

### Slow-Growth（`_enhanced_sampling.py`，sub-group 30）

301/302 共享基类 `_SlowgrowthPlotCmd`（注意：`output_name` 由父 `MenuGroup("30", output_name="slowgrowth")` 提供，命令自身无 `output_name`）：

- **文件发现**：
  - `_discover_restart_file(workdir)`：glob `*.restart`，排除 `_\d+\.restart` 检查点文件
  - `_discover_log_file(workdir)`：glob `*.LagrangeMultLog`
  - 两者均要求恰好 1 个匹配；否则回退到用户手动输入
- **轨迹信息预览**：`_print_sg_info()` 通过 `lazy_import` 获取 `ColvarMDInfo`，显示步数、时间步、CV 范围，并检测 NaN（溢出）步
- **参数采集**：restart path、log path、initial/final step、colvar ID、output dir
- **执行**：通过 `lazy_import` 调用 `slowgrowth_analysis`，`_plot_style` 由子类决定（`"quick"` 或 `"publication"`）

### Constrained TI（`_constrained_ti.py`，sub-group 31）

- **311 `TISingleDiagCmd`**：单点收敛诊断。复用 SG 的 `_discover_restart_file` / `_discover_log_file`。覆写 `_collect_all_params()`（SEM target 为可空 float，用 `prompt_str` + 手动转换）。调用 `standalone_diagnostics()`。
- **312 `TIFullAnalysisCmd`**：多点 TI 分析。调用 `discover_ti_points()` → `load_ti_series()` → `analyze_ti()`，含 dt 一致性校验。TI_DIR_PATTERN 限制为 `"ti_target"/"xi"/"auto"`。支持 `K.TI_REVERSE`（反向积分，ξ 降序，初态 = max ξ）。对所有约束点生成诊断图。
- **313 `TIConstPotCorrectionCmd`**：恒电势自由能修正（Nørskov）。先运行完整 312 分析，再从各 `ti_target_*/bader/` 提取系综平均 σ，通过 calibration mapper 外推 Φ，计算修正项并输出修正后自由能曲线。需要 calibration.json（硬错误）；缺 bader/ 时 WARN 并跳过修正。新增参数：`K.TARGET_SIDE`（aligned/opposed）、`K.CALIBRATION_JSON`。
