# `md_analysis.agent` 接口暴露约定（当前实现）

> 对应代码：`src/md_analysis/agent/__init__.py`
>
> 本文档定义 `md_analysis.agent` 的符号级公开接口与暴露边界。

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
  - 参数自动类型转换（str→Path, list→tuple/set）

- `list_tasks() -> list[dict[str, Any]]`
  - 枚举所有已注册任务，返回 `[{name, category, description, cli_codes}, ...]`

- `get_task_schema(task: str) -> dict[str, Any]`
  - 从目标函数签名自动生成 JSON Schema（OpenAI function calling 兼容）
  - 使用 `typing.get_type_hints()` 解析 stringified annotations
  - 通过 `TaskDef.param_descriptions` / `param_choices` 补充函数签名无法表达的元数据

### 2.2 数据结构

- `TaskResult` (frozen dataclass)
  - 字段：`success`, `task`, `outputs: dict[str, str]`, `summary: dict[str, Any]`, `error_type: str | None`, `errors: list[str]`, `warnings: list[str]`
  - 方法：`to_dict() -> dict[str, Any]`（JSON-serializable）

- `TaskDef` (frozen dataclass)
  - 字段：`name`, `category`, `description`, `handler: TaskHandler`, `target_fn: str`, `cli_codes`, `param_descriptions`, `param_choices`

- `TaskHandler` (Protocol, runtime_checkable)
  - 签名：`(params: dict[str, Any]) -> TaskResult`

### 2.3 注册 API

- `register(task_def: TaskDef) -> None`
  - 向全局注册表添加任务（幂等）
  - 外部包可调用此函数注册自定义任务

### 2.4 错误类型常量（`_core.py`）

- `ERROR_VALIDATION = "validation"`
- `ERROR_FILE_NOT_FOUND = "file_not_found"`
- `ERROR_ANALYSIS = "analysis"`
- `ERROR_INTERNAL = "internal"`

## 3. 当前注册任务（Phase 1：10 个）

| 任务名 | 目标函数 | CLI 编号 | 类别 |
|--------|---------|---------|------|
| `water_three_panel` | `main:run_water_analysis` | 105 | water |
| `potential_full` | `main:run_potential_analysis` | 216 | potential |
| `charge_surface` | `main:run_charge_analysis` | 221-223 | charge |
| `charge_tracked` | `main:run_tracked_charge_analysis` | 225 | charge |
| `charge_counterion` | `main:run_counterion_charge_analysis` | 226 | charge |
| `run_all` | `main:run_all` | — | composite |
| `calibration_fit_csv` | `CalibrationWorkflow:calibrate` | 231 | calibration |
| `calibration_predict` | `CalibrationWorkflow:predict_potential` | 233 | calibration |
| `slowgrowth_quick` | `SlowGrowthPlot:slowgrowth_analysis` | 301 | enhanced_sampling |
| `config_show` | `config:load_config` | 900 | meta |

## 4. 推荐导入方式

```python
from md_analysis.agent import dispatch, list_tasks, get_task_schema
from md_analysis.agent import TaskResult, TaskDef, register
```

## 5. 稳定性

- `dispatch` / `list_tasks` / `get_task_schema` / `TaskResult`：**Stable**
- `register` / `TaskDef` / `TaskHandler`：**Stable**（供外部扩展）
- 任务注册清单：**Evolving**（Phase 2 将扩展至 ~33 个任务）
