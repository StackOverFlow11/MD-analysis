# MD-Analysis 架构审查报告

> 审查日期：2026-04-04
> 最后更新：2026-04-06（P2–P7 全部完成）
> 审查范围：整体系统设计

## 当前设计优点

- src-layout 包结构清晰，依赖方向单向（CLI → 业务模块 → utils），无循环依赖
- CLI 框架 MenuGroup/MenuCommand 树 + flat index 快捷跳转 + ParamCollector 声明式参数采集，抽象级别恰当
- lazy_import 延迟加载 numpy/matplotlib/ase，CLI 启动零开销
- `_frame_source.py` 的 PotentialFrame 帧抽象统一了 continuous/distributed 两种输入模式
- MDAnalysisError 异常层次体系（11 个领域异常），CLI 统一捕获
- frozen dataclass 不可变数据模式贯穿全局
- Logging 符合 PEP 282（NullHandler + CLI 配置 StreamHandler）
- output_name 通过父链自动推导输出目录路径
- CSV 输出统一通过 `_write_csv`（dict rows）和 `_write_csv_from_arrays`（numpy arrays）两个入口
- matplotlib 绑定已从 6 处分析模块中剥离到专用 `_plot.py` / `plot.py` 文件

## 已解决问题

### ~~P1: CSV 写入方式不一致~~ (已修复 ffee59d)

新增 `_write_csv_from_arrays(path, columns)` 函数，替换 water 模块 7 处 `np.savetxt` + PhiZProfile 1 处自定义 `csv.writer`。所有 CSV 输出现通过 `_io_helpers` 两个统一入口。

### ~~P2: utils/__init__.py 巨型 re-export hub~~ (已修复)

清空 `utils/__init__.py` 的 re-export（`__all__ = []`）。所有调用方（包内/测试/外部）统一采用子模块直接导入路径（`from md_analysis.utils.constants import HA_TO_EV`、`from md_analysis.utils.StructureParser.LayerParser import ...` 等）。同步更新 `context4agent/architecture/modules/src/utils/` 文档与 `utils/CLAUDE.md` 约定。

### ~~P3: 两个 config.py 同名易混淆~~ (已修复)

`utils/config.py` → `utils/constants.py`（`git mv`）。所有 25 处 `utils.config` 导入（src + tests + docs）同步重命名为 `utils.constants`。用户持久化配置保留 `md_analysis/config.py` 名称，命名彻底解耦。

### ~~P4: Potential 命令的 sp_* 参数总是全量采集~~ (已修复)

新增 `ConditionalParam(inner, predicate)` 泛型包装器，将所有 6 个 Potential 命令的 4 个 `sp_*` 参数包裹起来，仅当 `ctx[K.INPUT_MODE] == "distributed"` 时提示用户，continuous 模式下静默应用默认值。新增 6 个单元测试覆盖 `ConditionalParam` 行为与 Potential 参数。

### ~~P5: 帧目录发现逻辑重复~~ (已修复)

新增 `utils/_frame_discovery.py` 共享私有模块，提供 `FRAME_DIR_STEP_TIME_RE` 编译正则、`extract_step_time_from_dirname()` 与 `discover_frame_dirs()` helper。`charge/Bader/_frame_utils.py` 与 `potential/_frame_source.py` 改为委托此模块（`constrained_ti/io.py` 基于 xi 值排序，语义不同，保留独立实现）。新增 13 个单元测试覆盖共享 helper。

### ~~P6: plotting 逻辑嵌入分析函数~~ (已修复)

matplotlib 绑定从 6 个分析模块拆分到 4 个新的 `_plot.py` 文件：
- `electrochemical/charge/Bader/_plot.py`（`plot_surface_charge`/`plot_single_side_charge`/`plot_tracked_charges`/`plot_counterion_charges`）
- `electrochemical/potential/_plot.py`（`plot_series_with_cumavg`/`plot_thickness_sensitivity`/`plot_phi_z_profile`）
- `water/_plot.py`（`plot_three_panel` + `_savgol_smooth_window5`）
- `enhanced_sampling/constrained_ti/plot.py`（新增 `plot_corrected_free_energy_profile`，`correction.py` 仅保留 re-export）

剩余 6 个 matplotlib 导入点均在专用 `_plot.py` / `plot.py` / `SlowGrowthPlot.py` 文件中，与计算逻辑彻底解耦。

### ~~P7: main.py 的 _nest 参数是设计气味~~ (已修复)

移除 `run_water_analysis`、`run_potential_analysis`、`run_charge_analysis`、`run_tracked_charge_analysis`、`run_counterion_charge_analysis` 五个函数的 `_nest`/`_nest_water` 参数。新契约：`output_dir` 即最终写入目录，不再自动前置 `water/`、`electrochemical/potential/` 等。`run_all` 显式按标准布局分派路径，CLI 移除 `_nest=False` 传参。同步更新 `md_analysis/CLAUDE.md` 与集成测试。

### ~~P8: Water.py 顶层 import matplotlib~~ (已修复 ffee59d)

matplotlib.ticker imports 从模块顶层移入 `plot_water_three_panel_analysis()` 函数体内，与项目惯例一致。（此次 P6 进一步将整段 plotting 逻辑迁入 `water/_plot.py`。）

## 优先级矩阵

| 编号 | 问题 | 影响 | 工作量 | 状态 |
|------|------|------|--------|------|
| P1 | CSV 写入不一致 | 中 | 低 | 已修复 |
| P2 | utils/__init__.py 巨型 re-export | 中 | 中 | 已修复 |
| P3 | 两个 config.py 同名 | 中 | 低 | 已修复 |
| P4 | sp_* 参数全量采集 | 中低 | 低 | 已修复 |
| P5 | 帧目录发现逻辑重复 | 中低 | 中 | 已修复 |
| P6 | plotting 嵌入分析函数 | 低 | 高 | 已修复 |
| P7 | _nest 参数设计气味 | 低 | 低 | 已修复 |
| P8 | Water.py 顶层 import matplotlib | 低 | 低 | 已修复 |
