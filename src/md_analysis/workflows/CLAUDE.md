# workflows — 开发备忘

## 定位

程序化入口的 **canonical** 模块。每个 `run_*` 函数封装一种分析或工作目录生成，
返回统一的 `WorkflowResult`（`artifacts` / `metadata` / `extra`）。`md_analysis.main`
是同一批名字的薄 re-export facade（79 行 import-only）；CLI / notebook / agent
handler 都通过本包调用业务流程，业务模块本身不持有 workflow 编排逻辑。

入口重构计划见 `context4agent/requirements/entrance_reconstruction.md`；
当前共 21 个 `run_*` + `WorkflowResult`，分 8 个子模块。

## 模块布局

| 文件 | 内容 |
|---|---|
| `models.py` | `WorkflowResult` frozen dataclass、`MissingArtifactError`、`require_artifacts_exist(result)` 校验 helper |
| `water.py` | `run_water_three_panel` |
| `potential.py` | `run_potential_full` |
| `charge.py` | `run_surface_charge` / `run_tracked_charge` / `run_counterion_charge` |
| `calibration.py` | `run_calibration_fit` / `run_calibration_predict` |
| `enhanced_sampling.py` | `run_slowgrowth_quick_plot` / `run_slowgrowth_publication_plot` / `run_ti_single_diagnostics` / `run_ti_full_analysis` / `run_ti_constant_potential_correction` |
| `scripts.py` | `run_bader_single/batch` / `run_ti_single/batch` / `run_potential_single/batch` / `run_sp_single/batch`（8 个工作目录生成入口） |
| `composite.py` | `run_interface_analysis`（water + potential，取代已删除的 `run_all`） |

## `WorkflowResult` 契约

- `artifacts: dict[str, Path]` —— 实际写盘的文件 / 工作目录。键稳定，
  `require_artifacts_exist(result)` 校验每个路径在磁盘上存在；批量
  workflow 用 `workdir_<i>` 索引键（非 zero-padded，按 enumeration 顺序）。
- `metadata: dict[str, Any]` —— 轻量 scalars（n_frames / method / 模式标
  志 / 子分析列表等），必须 JSON 友好。
- `extra: Any | None` —— 强类型 report 对象（`TIFullAnalysisReport` /
  `CalibrationFitReport` / `BaderGenBatchReport` 等）。如果 `extra` 实现
  了 `to_dict()`，`WorkflowResult.to_dict()` 会委托调用；否则原样透传。

## 约定

- **Path-normalize 在入口**：每个 `run_*` 函数体首部把 `xyz_path` / `output_dir`
  / `root_dir` / `md_out_path` / `calibration_json_path` 等 `Path | str | None`
  参数 `Path(...)` 包一遍；artifacts 字典里的值始终是 `Path`。
- **Lazy import 重业务**：模块顶层只 import 常量、`WorkflowResult` 和其他
  workflow facade；matplotlib / ase / 真正的业务模块（`water.*`、
  `electrochemical.*`、`enhanced_sampling.*`、`scripts.*`）一律在函数体
  内 lazy import，避免 CLI 启动时触发 numpy/matplotlib 加载。
- **不引入新业务**：workflows 只组织参数、输出路径、artifact 收集和
  metadata；底层的科学计算保留在 `water/` / `electrochemical/` /
  `enhanced_sampling/` / `scripts/`。如果发现需要在 workflow 层加算法
  逻辑，先回到底层模块加 helper，再让 workflow 调用。
- **不提交作业**：`scripts.py` 中的 8 个生成入口只写文件 / 工作目录，
  不调度 PBS / SLURM / VASP / CP2K。提交动作由集群侧脚本完成。
  `run_bader_*` 的 `generate_potcar=True` 会本地调 `vaspkit 103` 生成
  POTCAR 文件，仍属于纯文件生成。
- **Artifact 合并 + collision check**：composite workflow（`run_interface_analysis`）
  合并 leaf workflow 的 artifacts dict；如果 key 冲突直接 `RuntimeError`，
  不允许悄悄覆盖。
- **批量 workflow metadata 契约**：`scripts.py` 的所有 `run_*_batch`
  metadata 一律包含 `n_successful` / `n_skipped` / `n_failed` 三计数
  + `workdir_paths` 列表；当前底层 succeed-all-or-raise，所以
  `n_skipped == n_failed == 0`，但 key 占位以保未来扩展。

## 依赖方向（必须遵守）

```
cli / agent
   ↓
workflows
   ↓
water / electrochemical / enhanced_sampling / scripts
   ↓
engines / utils / exceptions
```

- workflows 不允许 import `cli` 或 `agent`（架构守卫扫描，零命中）
- 业务模块（water / electrochemical / 等）不允许 import workflows

## 测试位置

- 单测：`test/unit/workflows/`（10 个 test 文件，覆盖 leaf workflow
  artifact + metadata + extra 契约 + 错误路径 + mock orchestration）
- 集成：`test/integration/test_main.py` 通过 `md_analysis.main` re-export
  跑端到端验收（依赖 `data_example/potential/dense/` 真实数据）

## 陷阱与历史 Bug

- **Bug a773461**（Phase 2 fix）：`run_counterion_charge` 早期总把
  `counterion_charge_png` 加入 artifacts，但底层在 0 counterion 帧时不写
  PNG。修复：用 `*_with_report` wrapper，根据 `result.png_path is not None`
  条件加入。
- **Bug e54ef5b**（Phase 4.2 fix）：`run_ti_constant_potential_correction`
  的 `point_slice` 切片条件原写 `if point_slice is not None`，但底层
  `run_ti_full_from_root` 同时把 `""` 和 `None` 视为 "不切片"。修复：
  对齐为 `if point_slice is not None and point_slice != ""`，并加 mock
  级 orchestration 测试覆盖空串路径。
- `extra` 字段必须包含 `to_dict()` 方法才能 JSON 序列化；如果底层 report
  没有 `to_dict()`，把它放进 `extra` 时调用方需要自己处理或写到磁盘前
  序列化。
