# `md_analysis.agent` 接口暴露约定（当前实现）

> 对应代码：`src/md_analysis/agent/__init__.py`
>
> 本文档定义 `md_analysis.agent` 的符号级公开接口与暴露边界。

> ⚠️ **utils/engines 重构进行中**（见 `context4agent/requirements/current_reconstructions.md`）。
> Phase 2-8 期间 `agent` 层不作为公开接口承诺面：
> - `target_fn` 字符串会随模块迁移（`utils/RestartParser` → `utils/formats/cp2k_colvar` 等）同步更新
> - Phase 5 之后 dataclass 名称会有 rename：`ColvarRestart` → `ConstraintMetadata`、`LagrangeMultLog` → `LambdaSeries`
> - 本文档下方细节将在 Phase 9 集中重写
>
> 重构期间集成请以代码为准，不要基于本文档当前细节做下游对接。

## 1. 接口角色定义

- `md_analysis.agent` 是面向 AI agent 的非交互式编程入口。
- 与 `cli/`（面向人类的交互式菜单）和 `main.py`（面向脚本的编程 API）互补。
- 所有公开符号通过 `__init__.py` 的 `__all__` 导出。
- 不包含分析逻辑——仅做参数转换、任务分发、结果包装。

## 2. 当前公开接口清单

### 2.1 核心调度函数

- `dispatch(task: str, params: dict[str, Any] | None = None) -> TaskResult`
  - 执行指定任务，返回统一结果
  - 分层异常捕获：validation / file_not_found / analysis / internal
  - 参数类型转换：
    - 有 contract 时走 `_coerce_params_from_contract()`（规则仅依赖 `FieldSpec.path_kind` 和 `json_schema["type"]` 顶层；不解析 `FieldSpec.type` 字符串）
    - 无 contract 时走 legacy `_coerce_params()`：从 `target_fn` 注解推导（str→Path, list→tuple/set）
  - 异常分类：有 contract 时优先按 `contract.exceptions` 有序匹配（**先具体子类后父类**），用 `ExceptionMapping.error_type` 覆盖默认分类；无匹配或无 contract 时走默认（ValueError/TypeError→validation, FileNotFoundError/PermissionError→file_not_found, MDAnalysisError→analysis, others→internal）

- `list_tasks() -> list[dict[str, Any]]`
  - 枚举所有已注册任务，返回 `[{name, category, description, cli_codes}, ...]`

- `get_task_schema(task: str) -> dict[str, Any]`
  - 返回 OpenAI function-calling 兼容形状：`{name, description, parameters: {type, properties, required}}`
  - **Contract-first**：当 `TaskDef.contract is not None` 时，委托给 `contract.to_agent_schema(name, description)`；**返回形状不变**
  - Fallback（无 contract）：从 `target_fn` 签名自动推导 —— 使用 `typing.get_type_hints()` 解析 stringified annotations；通过 `TaskDef.param_descriptions` / `param_choices` 补充元数据

### 2.2 数据结构

- `TaskResult` (frozen dataclass)
  - 字段：`success`, `task`, `outputs: dict[str, str]`, `summary: dict[str, Any]`, `error_type: str | None`, `errors: list[str]`, `warnings: list[str]`
  - 方法：`to_dict() -> dict[str, Any]`（JSON-serializable）

- `TaskDef` (frozen dataclass)
  - 字段：`name`, `category`, `description`, `handler: TaskHandler`, `target_fn: str`, `cli_codes`, `param_descriptions`, `param_choices`, `contract: TaskContract | None = None`, `reference_fn: str | None = None`
  - `contract`：Tools 层结构化契约（优先于 `target_fn` 做 schema / coercion / 异常分类）
  - `reference_fn`：可选的 dotted path，指向签名与 `contract.inputs` 一致的真实函数（供静态签名校验使用；默认回退 `target_fn`）

- `TaskHandler` (Protocol, runtime_checkable)
  - 签名：`(params: dict[str, Any]) -> TaskResult`

### 2.5 Tools 层契约数据结构（`_contracts.py`）

- `FieldSpec` (frozen dataclass)
  - 字段：`description`, `json_schema: dict`（权威机器可读 schema，draft-07），`type: str`（人读标注，不参与生成/转换），`required`, `default`, `choices`, `category: Literal["artifact", "metric", "raw_model"]`, `unit`, `shape`, `path_kind: Literal["file", "dir", "glob"] | None`

- `ExceptionMapping` (frozen dataclass)
  - 字段：`exception_fqn: str`（FQN，如 `"md_analysis.scripts.TIGen.TIGenError"`）、`triggered_by`, `error_type: Literal["validation", "file_not_found", "analysis", "internal"]`, `user_visible`（预留给未来 MCP server 层，dispatch 当前忽略）

- `TaskContract` (frozen dataclass)
  - 字段：`inputs`, `outputs_artifacts`, `outputs_metrics`, `outputs_raw_model`, `preconditions`, `side_effects`, `exceptions`（有序：先具体子类后父类）
  - 方法：
    - `to_agent_schema(name, description) -> dict`：OpenAI function-calling 形状
    - `to_mcp_tool_schema(name, description) -> dict`：MCP `inputSchema` 形状（供未来 MCP server 封装层）

### 2.3 注册 API

- `register(task_def: TaskDef) -> None`
  - 向全局注册表添加任务（幂等）
  - 外部包可调用此函数注册自定义任务

### 2.4 错误类型常量（`_core.py`）

- `ERROR_VALIDATION = "validation"`
- `ERROR_FILE_NOT_FOUND = "file_not_found"`
- `ERROR_ANALYSIS = "analysis"`
- `ERROR_INTERNAL = "internal"`

## 3. 当前注册任务（14 个）

| 任务名 | 目标函数 | CLI 编号 | 类别 | 有 contract |
|--------|---------|---------|------|---|
| `water_three_panel` | `main:run_water_analysis` | 105 | water | — |
| `potential_full` | `main:run_potential_analysis` | 216 | potential | — |
| `charge_surface` | `main:run_charge_analysis` | 221-223 | charge | — |
| `charge_tracked` | `main:run_tracked_charge_analysis` | 225 | charge | — |
| `charge_counterion` | `main:run_counterion_charge_analysis` | 226 | charge | — |
| `run_all` | `main:run_all` | — | composite | — |
| `calibration_fit_csv` | `CalibrationWorkflow:calibrate` | 231 | calibration | — |
| `calibration_predict` | `CalibrationWorkflow:predict_potential` | 233 | calibration | — |
| `slowgrowth_quick` | `SlowGrowthPlot:slowgrowth_analysis` | 301 | enhanced_sampling | — |
| `ti_full_analysis` | `constrained_ti.workflow:run_ti_full_from_root` | 312 | enhanced_sampling | ✅ |
| `bader_gen_batch` | `scripts.BaderGen:generate_bader_batch_with_report` | 412 | scripts | ✅ |
| `ti_gen_batch` | `scripts.TIGen:generate_ti_batch_with_report` | 422 | scripts | ✅ |
| `sp_gen_batch` | `scripts.SpGen:batch_generate_sp_workdirs` | 442 | scripts | — |
| `config_show` | `config:load_config` | 900 | meta | — |

## 4. 推荐导入方式

```python
from md_analysis.agent import dispatch, list_tasks, get_task_schema
from md_analysis.agent import TaskResult, TaskDef, register
# 扩展新任务时：
from md_analysis.agent._contracts import (
    FieldSpec, ExceptionMapping, TaskContract,
)
```

## 5. 稳定性

- `dispatch` / `list_tasks` / `get_task_schema` / `TaskResult`：**Stable**（`get_task_schema()` 返回形状在引入 contract 后**保持不变**）
- `register` / `TaskDef` / `TaskHandler`：**Stable**（供外部扩展）
- `FieldSpec` / `ExceptionMapping` / `TaskContract`：**Evolving**（MVP 阶段，字段可能调整；在 Resources 层引入前保留弹性）
- `TaskDef.contract` / `TaskDef.reference_fn`：**Evolving**（迁移期；旧任务暂未迁移）
- 任务注册清单：**Evolving**（Phase 2 将扩展至 ~33 个任务）

## 6. Contract-backed 任务说明

当前 3 个任务带完整 `TaskContract`：

- **`ti_gen_batch`**（CLI 422）：pass-through 到 `scripts.TIGen.generate_ti_batch_with_report`；含文件预检 + collision check
- **`ti_full_analysis`**（CLI 312）：composite，backed by 真实 wrapper `constrained_ti.workflow.run_ti_full_from_root`（handler 薄化）。`summary.per_point` 是**未来 Resources 层**用于判定业务失败类别（非平衡漂移 / N_eff 太少 / 遍历性假阳等）的信号源 —— 字段：`point_index`、`xi`、`n_analyzed`、`time_start_fs`、`time_end_fs`、`time_total_fs`、`tau_corr`、`n_eff`、`sem_final_au`、`sem_max_au`、`geweke_z`、`drift_D`、`passed`、`failure_reasons`
- **`bader_gen_batch`**（CLI 412）：pass-through 到 `scripts.BaderGen.generate_bader_batch_with_report`；**只**准备 VASP Bader 工作目录（POSCAR/INCAR/KPOINTS [+ POTCAR via `vaspkit 103`] + script.sh），**不**提交任务、**不**解析 Bader 输出、**不**做恒电势修正。`summary` 暴露 `n_frames`/`frame_indices`/`steps`/`times_fs`/`generate_potcar`
