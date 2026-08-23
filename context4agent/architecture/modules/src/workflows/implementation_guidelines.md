# `md_analysis.workflows` 内部实现准则（当前实现口径）

> 适用范围：`src/md_analysis/workflows/`。详细开发备忘（陷阱与历史 Bug）见 `src/md_analysis/workflows/CLAUDE.md`。

## 1. 职责边界

- workflows 是 **facade 层**：每个 `run_*` 只做参数组织、路径规范化、输出目录创建、artifact 收集与 metadata 组装，最后返回 `WorkflowResult`。
- **业务逻辑不下沉到 workflows**：底层的科学计算（数据解析、统计算法、绘图细节）保留在 `water/` / `electrochemical/` / `enhanced_sampling/` / `scripts/`。若发现需要在 workflow 层加算法逻辑，先回到底层模块加 helper，再让 workflow 调用。
- **不提交作业**：`scripts.py` 的 8 个生成入口只写文件 / 工作目录，不调度 PBS / SLURM / VASP / CP2K；`run_bader_*` 的 `generate_potcar=True` 本地调 `vaspkit 103` 生成 POTCAR，仍属纯文件生成。

## 2. 模块布局

| 文件 | 内容 |
|---|---|
| `models.py` | `WorkflowResult` frozen dataclass、`MissingArtifactError`、`require_artifacts_exist(result)` |
| `water.py` | 水分析 5 个 `run_*`（三联图 + 4 个单步） |
| `potential.py` | 电势分析 6 个 `run_*`（`run_potential_full` 复合 + 5 个单步；全部支持 continuous/distributed 双输入模式） |
| `charge.py` | 电荷分析 3 个 `run_*`（surface / tracked / counterion） |
| `calibration.py` | 标定 2 个 `run_*`（fit / predict） |
| `enhanced_sampling.py` | 增强抽样 5 个 `run_*`（SG 2 + constrained TI 3） |
| `scripts.py` | 工作目录生成 8 个 `run_*`（Bader / TI / Potential / SP 各 single+batch） |
| `composite.py` | `run_interface_analysis`（water + potential，取代已删除的 `run_all`） |
| `__init__.py` | 按域 re-export 全部公开符号；`__all__` 是接口权威 |

当前共 30 个 `run_*` + 3 个模型符号。

## 3. 依赖方向约束（必须遵守）

```
cli
   ↓
workflows
   ↓
water / electrochemical / enhanced_sampling / scripts
   ↓
engines / utils / exceptions
```

- workflows **不允许 import `cli`**（重构期架构守卫扫描确认零命中）。
- 业务模块（water / electrochemical / 等）**不允许 import workflows**。
- workflows 顶层只 import 常量（各域 `config`）、`WorkflowResult` 与兄弟 workflow facade（`composite.py` 组合 `water.py` / `potential.py`；`enhanced_sampling.py` 的 correction 组合同包 `run_ti_full_analysis`）；唯一例外是 `calibration.py` 因签名类型标注需要顶层 `import numpy`。

## 4. 实现约定

1. **Path-normalize 在入口**：每个 `run_*` 函数体首部把 `xyz_path` / `output_dir` / `root_dir` / `md_out_path` / `calibration_json_path` 等 `Path | str | None` 参数 `Path(...)` 包一遍；artifacts 字典的值始终是 `Path`。
2. **Lazy import 重依赖**：matplotlib / ase / 真正的业务模块（`water.*`、`electrochemical.*`、`enhanced_sampling.*`、`scripts.*`）一律在函数体内 lazy import，避免 CLI 启动时触发 numpy/matplotlib 加载。
3. **artifacts 键稳定**：键由 workflow 定义并在 docstring / 镜像接口文档中登记；批量 workflow 用 `workdir_<i>`（十进制枚举索引，非 zero-padding）。`test/unit/cli/test_workflow_wiring.py` 锁定了 CLI 读取的键名。
4. **metadata 轻量且 JSON 友好**：帧数、方法、模式标志等 scalars；批量入口固定带 `n_successful` / `n_skipped` / `n_failed` + `workdir_paths`。
5. **extra 强类型 report**：底层有 `*Report` / `*Result` 对象时放 `extra`；实现 `to_dict()` 的对象可被 `WorkflowResult.to_dict()` 序列化，否则由调用方处理。
6. **条件 artifact**：底层可能不产出的文件（counterion PNG、TI 单点诊断产物）按实际产出与否条件加入 artifacts，避免 `require_artifacts_exist` 把合法空跑判为失败。
7. **composite 合并 + 冲突检查**：`run_interface_analysis` 合并 leaf artifacts 时若键冲突直接 `RuntimeError`，不允许悄悄覆盖。
8. **参数透传语义**：
   - `run_ti_full_analysis` / `run_ti_constant_potential_correction` 的 `strict=True` 为 agent-facing 默认；`run_ti_constant_potential_correction` 必须把 `strict` 同时透传给内部 `run_ti_full_analysis` 和二次 `discover_ti_points`，否则两次点集错位、per-point Bader 对齐漂移。
   - `run_ti_batch` 的 `colvar_id` / `overwrite` 进 metadata；`verbose` 是 UI-only，不进 metadata。
   - `run_surface_charge` 在任何文件系统副作用**之前**校验 `target_side`，避免非法值泄漏空目录。

## 5. 测试位置

- 单测：`test/unit/workflows/`（8 个测试文件，覆盖 models 契约 + 各域 artifact/metadata/extra 契约 + 错误路径 + mock orchestration）
- CLI 契约：`test/unit/cli/test_workflow_wiring.py`（monkey-patch `lazy_import`，锁定 CLI→workflows 的模块路径、kwargs 透传与 artifact 键读取）
- 集成：`test/integration/test_main.py` 经 `md_analysis.main` re-export 跑端到端验收（依赖 `data_example/potential/dense/` 真实数据）

## 6. 同步要求

改动本包时必须同步更新：

- 同目录 `interface_exposure.md`（`__all__` 清单、artifacts 键、目录布局约定）
- `src/md_analysis/main.py` re-export 清单（公开符号增删时）
- `src/md_analysis/workflows/CLAUDE.md`（开发备忘，含陷阱记录）
- CLI 菜单接线变化时：`context4agent/architecture/modules/src/cli/` 两份文档
