# `md_analysis.agent` 实现指南（当前实现）

> 对应代码：`src/md_analysis/agent/`

## 1. 职责边界

- **是**：参数转换（str→Path）、任务注册与分发、JSON Schema 生成、统一异常包装
- **不是**：分析逻辑（委托给 `main.py` 或子模块函数）

## 2. 模块结构

| 文件 | 职责 |
|------|------|
| `__init__.py` | re-export public API + 触发 `_handlers` 注册 |
| `_core.py` | TaskResult, TaskHandler Protocol, TaskDef, 注册表 |
| `_dispatch.py` | dispatch(), get_task_schema(), _coerce_params(), _annotation_to_schema() |
| `_handlers.py` | _make_handler() 工厂 + Phase 1 任务注册 |

## 3. 依赖方向

```
agent/ → main.py → water/ | electrochemical/ | enhanced_sampling/
agent/ → utils/cell_resolver.py → utils/RestartParser/CellParser.py
agent/ → exceptions.py (MDAnalysisError)
```

- `agent/` 不被任何其他包反向依赖
- `agent/` 不被 `md_analysis/__init__.py` re-export（避免 import 开销）

## 4. 关键设计决策

### 4.1 Schema 自动生成（非手写 ParamDef）

`get_task_schema()` 使用 `inspect.signature()` + `typing.get_type_hints()` 从目标函数签名自动推导 JSON Schema。`TaskDef` 仅保留 `param_descriptions` 和 `param_choices` 作为补充注解。

**理由**：避免在函数签名和参数定义之间维护双份声明（审查 Round 1 核心问题）。

### 4.2 Handler 工厂模式

`_make_handler(task_name, target_fn_path)` 生成直通 handler，内部通过 `importlib.import_module()` 延迟导入目标函数。

**规则**：handler 内部禁止顶层 import 分析模块（numpy/matplotlib/ase），必须在函数体内延迟导入。

### 4.3 分层异常捕获

异常处理集中在 `dispatch()` 层，handler 内部不做 try-except：
- `ValueError/TypeError` → `ERROR_VALIDATION`
- `FileNotFoundError` → `ERROR_FILE_NOT_FOUND`
- `PermissionError` → `ERROR_FILE_NOT_FOUND`
- `MDAnalysisError` → `ERROR_ANALYSIS`
- `Exception` → `ERROR_INTERNAL`（含 `logger.error(exc_info=True)`）

### 4.4 handler 调用层级策略

优先调用最高可用编排函数：有 `main.py` 的 `run_*()` 就用它，无则直调子模块函数。

### 4.5 `_normalize_outputs` 返回类型分支

`_normalize_outputs`（`_handlers.py`）将 handler 返回值转为 `dict[str, str]`：

| 返回类型 | 转换结果 |
|---|---|
| `dict` | 直接字符串化 values |
| `Path` | `{"output": str(path)}` |
| `list` | `{"workdir_0": ..., "workdir_1": ..., ...}`（批量脚本 `list[Path]`，如 `sp_gen_batch`） |
| 带 `csv_path` 属性 | `{"csv": ..., "png": ...}`（若同名 PNG 存在） |
| 其他 | `{}` + debug 日志 |

添加新的返回类型分支时需确保不影响既有任务的转换行为。

## 5. 新增任务检查清单

1. 在 `_handlers.py` 中用 `_make_handler()` 工厂（简单直通）或手写 handler（需多步编排，如 `_handle_ti_full_analysis`）
2. 调用 `register(TaskDef(...))` 注册
3. 可选：`param_descriptions` / `param_choices` 补充 schema 注解
4. 更新本文档的任务清单
5. 更新 `interface_exposure.md` 第 3 节

## 6. 契约同步

本包变更时需同步：
- `context4agent/architecture/modules/src/agent/interface_exposure.md`
- `context4agent/architecture/README.md`（第 7a 节）
- `context4agent/requirements/short_term.md`（编程入口清单）
