# 用户诉求（维护说明）

本目录用于记录"你希望这个 CP2K 后处理包实现什么"。

## 文件说明

- `short_term.md`：**近期诉求**（持续更新）。特点：明确、可执行、可拆解成任务。
- `long_term.md`：**远期诉求**（按需更新）。特点：路线图/愿景/未来可能引入的方向。
- `overall_reconstruction_plan.md`：**当前大规模破坏性重构的总体顺序**。Phase 0 → Phase 7 时序轴。
- `target_structure.md`：**目标模块职责边界**。配合 `overall_reconstruction_plan.md` 提供静态架构边界。

## 更新规则

- 每次提出新的近期目标、优先级变化、或"下一步做什么"，都更新 `short_term.md`。
- 只有当明确要求"更新远期诉求/路线图/愿景"时，才更新 `long_term.md`。
- 大规模重构期间，时序变更更新 `overall_reconstruction_plan.md`；边界变更更新 `target_structure.md`。
- 架构与实现细节记录在 `context4agent/architecture/` 中，不放在本目录。
- **本目录不维护完整 API 手册**：`short_term.md` 只写「能力分组摘要 +
  关键不变量」；public 函数清单 / 参数 / 菜单号一律**引用代码权威位置**
  （`src/md_analysis/workflows/__init__.py`、agent registry/`get_task()`、
  `src/md_analysis/cli/`），不手写易漂移长清单。已确认的科学口径（公式/
  单位/CSV 列名）属例外，必须在 `short_term.md` 保留原样。
