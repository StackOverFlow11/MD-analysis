# `md_analysis.agent` 实现指南（当前实现）

> 对应代码：`src/md_analysis/agent/`

> ⚠️ **入口重构主体已完成。** 当前文档描述完成后的 contract-backed agent 实现。
> 后续大规模破坏性重构按 `context4agent/requirements/overall_reconstruction_plan.md`
> 推进，迁移到新架构后需同步更新本文档。

## 1. 职责边界

- **是**：参数转换（str→Path）、任务注册与分发、JSON Schema 生成、统一异常包装
- **不是**：分析逻辑（委托给 `main.py` 或子模块函数）

## 2. 模块结构

| 文件 | 职责 |
|------|------|
| `__init__.py` | re-export public API + 触发 `_handlers` 注册 |
| `_core.py` | TaskResult, TaskHandler Protocol, TaskDef（含 `contract` / `reference_fn`）, 注册表 |
| `_contracts.py` | Tools-layer 结构化契约：FieldSpec, ExceptionMapping, TaskContract（双 schema 出口） |
| `_dispatch.py` | dispatch(), get_task_schema(), contract-aware coercion & 异常分类, legacy `_coerce_params()` |
| `_handlers.py` | _make_handler() 工厂 + 任务注册（contract 实例紧邻 register 调用） |

## 3. 依赖方向

```
agent/ → main.py → water/ | electrochemical/ | enhanced_sampling/
agent/ → utils/io/cell_resolver.py → utils/formats/cp2k_cell.py
agent/ → exceptions.py (MDAnalysisError)
```

- `agent/` 不被任何其他包反向依赖
- `agent/` 不被 `md_analysis/__init__.py` re-export（避免 import 开销）

## 4. 关键设计决策

### 4.1 Schema 生成：contract-first，签名推导为 fallback

新任务应优先使用 **`TaskContract`**（`_contracts.py`）声明输入/输出契约：
- `get_task_schema()` 在 `task_def.contract is not None` 时调用 `contract.to_agent_schema()`，返回 OpenAI `{name, description, parameters}` 形状（**公开 API 形状不变**）
- 未来 MCP server 层另行调用 `contract.to_mcp_tool_schema()` 产出 `inputSchema` 形状
- `FieldSpec.json_schema`（draft-07 片段）是机器可读权威；`FieldSpec.type` 仅作人读标注，**不参与**生成/转换

旧任务（无 contract）走 legacy 路径：`inspect.signature()` + `typing.get_type_hints()` 从 `target_fn` 推导，`param_descriptions`/`param_choices` 补充元数据。两条路径共存，新任务**必须**走 contract。

**理由**：避免在函数签名和参数定义之间维护双份声明（审查 Round 1 核心问题）；同时解决 composite handler 签名与真实参数不一致导致的 schema 错误，并为未来 MCP Resources/Prompts 层提供结构化信号（`outputs_metrics` 暴露业务判定所需字段）。

### 4.2 Handler 工厂模式

`_make_handler(task_name, target_fn_path)` 生成直通 handler，内部通过 `importlib.import_module()` 延迟导入目标函数。

**规则**：handler 内部禁止顶层 import 分析模块（numpy/matplotlib/ase），必须在函数体内延迟导入。

### 4.3 分层异常捕获

异常处理集中在 `dispatch()` 层，handler 内部不做 try-except。分两路：

**有 contract**：按 `contract.exceptions` 元组**有序**匹配（先具体子类后父类），用 `ExceptionMapping.error_type` 覆盖默认分类。`exception_fqn` 是 fully-qualified path（如 `"md_analysis.scripts.TIGen.TIGenError"`、`"builtins.FileNotFoundError"`），通过 `_resolve_exception_class()` 懒加载；解析失败的条目跳过（记录 DEBUG）。

**无 contract 或未匹配**：fallback 到默认分类：
- `ValueError/TypeError` → `ERROR_VALIDATION`
- `FileNotFoundError` → `ERROR_FILE_NOT_FOUND`
- `PermissionError` → `ERROR_FILE_NOT_FOUND`
- `MDAnalysisError` → `ERROR_ANALYSIS`
- `Exception` → `ERROR_INTERNAL`（含 `logger.error(exc_info=True)`）

**典型用法**：`InsufficientSamplingError`（继承 `MDAnalysisError`）对代码而言是分析错误，但对 agent 而言语义更贴近"输入不满足前置条件"；此时 contract 可声明 `error_type="validation"` 重写分类。

### 4.4 handler 调用层级策略

优先调用最高可用编排函数：有 `main.py` 的 `run_*()` 就用它，无则直调子模块函数。

### 4.5 `_normalize_outputs` 返回类型分支

`_normalize_outputs`（`_handlers.py`）将 handler 返回值转为 `dict[str, str]`：

| 返回类型 | 转换结果 |
|---|---|
| `dict` | 直接字符串化 values |
| `Path` | `{"output": str(path)}` |
| `list` | `{"workdir_0": ..., "workdir_1": ..., ...}`（批量脚本 `list[Path]`，如 `sp_gen_batch`） |
| 带 `workdirs` 属性（tuple/list）| `{"workdir_0": ..., ...}`（如 `TIGenBatchReport`） |
| 带 `csv_path` 属性 | `{"csv": ..., "png": ...}`（若同名 PNG 存在） |
| 其他 | `{}` + debug 日志 |

添加新的返回类型分支时需确保不影响既有任务的转换行为。

## 5. 新增任务检查清单

1. 在 `_handlers.py` 中：
   - 简单直通：`_make_handler(task_name, "module:fn")`
   - 多步编排：手写 handler（参考 `_make_ti_gen_batch_handler()` / `_ti_full_analysis_handler` / `_make_bader_gen_batch_handler()`）
   - **Composite 原则**：如果真实业务函数签名和 contract 不一致（如老 `ti_full_analysis` 的 `analyze_ti` 接受 `xi_values` / `lambda_series_list` / `dt`），**先抽一个签名与 contract 匹配的真实 wrapper 函数**，让 `target_fn` 指向它；避免用 handler 内部重新组装参数
   - **已落地范例**：
     - `ti_gen_batch`（CLI 422）：pass-through 脚本准备型。wrapper
       `scripts.TIGen:generate_ti_batch_with_report` 做预检 + collision check
       + 结构化 report。
     - `bader_gen_batch`（CLI 412）：同类 pass-through 脚本准备型。wrapper
       `scripts.BaderGen:generate_bader_batch_with_report` 做 xyz/script 预检
       + 帧元数据收集（frame_indices / steps / times_fs）；`generate_potcar=True`
       会调 `vaspkit 103` 作为本地副作用。**不提交作业、不解析 Bader 输出**
       —— 作业提交归 pbs-auto / 用户，结果解析归 `charge_*` 类 API。
     - `ti_full_analysis`（CLI 312）：composite 改造样例。wrapper
       `enhanced_sampling.constrained_ti.workflow:run_ti_full_from_root`，签名
       和 contract 一致；handler `_ti_full_analysis_handler` 只做
       artifacts/metrics 路由（report 的文件路径 → `TaskResult.outputs`、
       JSON 可序列化数值 → `TaskResult.summary`），不含业务逻辑。
2. 定义 `TaskContract`（紧邻 `register()` 调用写入同一段）：
   - `inputs`：每个 `FieldSpec` 必须填 `description` + `json_schema`（draft-07）
   - `outputs_artifacts` / `outputs_metrics` / `outputs_raw_model` 三分法
   - `preconditions`：仅输入级**硬约束**（违反必抛异常）；软约束归 Resources 层
   - `side_effects`：写哪些文件 / 改哪些状态
   - `exceptions`：有序 —— **先具体子类后父类**，`exception_fqn` 用 FQN，`error_type` 按 agent 视角选 validation/file_not_found/analysis/internal
3. 调用 `register(TaskDef(..., contract=..., target_fn=..., handler=...))`
4. 为不可 MCP 返回的内部对象（如 `TIReport`）保留 `outputs_raw_model` 声明 —— 供未来 Resources 层引用字段路径
5. 新增测试：
   - `test/unit/agent/test_contracts.py`：schema 形状、coercion 规则、异常映射顺序
   - 任务级测试：dispatch() 成功/失败路径
6. 更新本文档的任务清单 + `interface_exposure.md` 第 3 节

## 5a. md-analysis Agent Scope 边界（重要）

md-analysis 的 agent tasks 只覆盖 **分析 + 工作目录准备**：

- ✅ 本包内做：文件解析、数值分析、绘图、VASP/CP2K 输入目录生成、本地
  `vaspkit` 调用（POTCAR 生成）
- ❌ **不做**：提交 PBS / SLURM 作业、监控队列状态、解析作业输出、读取集群
  负载

作业提交归 **pbs-auto**（独立 MCP server，未来 Prompts 层串联）。`bader_gen_batch`
是典型例子：只写目录，不提交；`ti_bader_status`、"恒电势修正"这类跨阶段
task 需要多个 MCP server 协作，留给 Prompts 层。

## 6. Tools 层 vs Resources 层边界（重要）

`TaskContract` 是**机械性**契约 —— 只描述"调用发生了什么"，不含业务判断：

**不放入 `TaskContract` 的内容**：
- "调用成功但业务不合格" 的 failure mode 枚举（如 `NOT_CONVERGED_NEQ_DRIFT` / `CONVERGED_BUT_SHORT_TRAJ`）→ 未来 Resources 层
- 软约束 / 领域假设（如 "1–2 ps 收敛可能非遍历" / "σ_ref 应取 IS/FS 中点"）→ 未来 Resources 层
- 跨任务 workflow 决策树 → 未来 Prompts 层

**Tools 层要暴露的**：足够的 `outputs_metrics` 信号，让未来 Resources 层能基于它们判定业务失败并给出 `suggested_action`。例：`ti_full_analysis` 的 `outputs_metrics.per_point[*]` 要暴露 `time_total_fs` / `n_eff` / `failure_reasons`，否则 Resources 层无从判断"遍历性假阳"等情形。

## 7. 契约同步

本包变更时需同步：
- `context4agent/architecture/modules/src/agent/interface_exposure.md`
- `context4agent/architecture/README.md`（第 7a 节）
- `context4agent/requirements/short_term.md`（编程入口清单）
