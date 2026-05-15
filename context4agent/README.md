# History / Context

本目录用于**记录与项目持续相关的上下文**，供后续开发/Agent 协作时快速对齐，避免因对话过长导致遗忘。

## 使用规则（重要）

- **近期诉求（Short-term）持续更新**：每次提出"最近要做什么/优先级变化/本周目标"，都更新 `requirements/short_term.md`。
- **远期诉求（Long-term）按需更新**：只有明确要求"更新远期诉求/路线图"时，才更新 `requirements/long_term.md`。
- **临时重要上下文先放 temp**：不确定归属时先写入 `temp/`，后续再搬运并在原条目标记"已迁移"。
- **新增文件先加索引**：在本 README 的"目录索引"里补一条，方便检索。
- **未协商不得固化约定**：任何**未与用户协商确定**的接口/口径/命名/数据契约/单位约定，**不得自行采用或写入** `context4agent/` 对应部分；只能先记录为"待讨论/待确认"。

## 文档分层（治理守则）

`context4agent/` 是「agent 协作与架构决策手册」，**不是**完整 API 文档、
函数参数手册、测试报告归档或提交流水账。读者是后续协作 agent / 未来的
自己，不是最终包用户。它只记录**后续 agent 做决策时必须知道、且不能仅靠
读当前代码稳定推断**的信息。

四层职责：

| 层 | 目录 | 角色 | 更新条件 |
|---|---|---|---|
| frozen | `architecture/` | 稳定架构边界 / 数据流方向 / 跨模块契约 / 已确认科学口径 | 仅当架构分层、数据契约、科学定义/单位变化，**且用户明确批准**；普通实现迁移/测试补强/bugfix **不更新** |
| active | `requirements/` | 当前阶段目标、时序、需求状态、关键不变量 | 每次提出新近期目标/优先级变化/下一步;长 API 清单只引代码权威位置，不手写 |
| historical | 已完成阶段在 active 内的**摘要短条目** | 已批准行为变化 / 长期例外 / 影响后续决策的历史坑 | 用短条目记录"为什么"，不保留每轮实现细节 |
| temp | `temp/`（git untracked） | 局部 plan 草稿 / handoff / 未定归属上下文 | 默认不入库；沉淀后压成短条目搬入正式文件，原处标「已迁移」 |

**权威来源原则**：public API 导出以 `src/md_analysis/workflows/__init__.py`
为准；CLI 行为以 `src/md_analysis/cli/` + 测试为准。
`context4agent/` 只记录**入口策略、长期例外、已批准行为变化、科学口径**，
不复制易漂移的逐函数/逐菜单号长清单。

**固化前三问**：(1) 这条信息 3 个月后还会影响 agent 决策吗？(2) 它无法
从当前代码/测试稳定推断吗？(3) 涉及约定时是否已被用户确认（非 agent
临时判断）？至少满足前两条、涉约定必满足第三条，才适合写入。

**行为变化三分类**:accidental drift（应阻止/修复，不记）、
default-preserving migration（参数默认值保持旧行为，如 `strict`/`overwrite`
——记短条目）、approved tightening（更合理但改旧行为，须用户批准
——记短条目）。只记后两类中会影响未来决策的事实，不贴完整计划/长测试输出。

## 目录索引

- `user_level.md`：用户当前水平、习惯、偏好（面向协作效率）
- `requirements/`
  - `short_term.md`：近期诉求（持续更新）
  - `long_term.md`：远期诉求（按需更新）
  - `README.md`：诉求拆分原则与维护说明
- `architecture/`
  - `README.md`：代码架构总览（目录/模块职责/数据流）
  - `modules/`：子模块说明（含 `src/` 目录层级对应的接口暴露/实现准则文档）
    - `modules/README.md`：模块文档治理硬约束（目录镜像与双文档）
    - `modules/data_contract.md`：跨模块核心数据契约（形状/单位/CSV 列头）
    - `modules/glossary_units.md`：术语与单位约定汇总
- `temp/`：重要上下文临时寄存区

## 约定

- Markdown 公式格式：行内使用 `$...$`，块级使用 `$$...$$`；不使用 `\(...\)` 与 `\[...\]` 作为公式定界符
- 术语与单位制：已在 `architecture/modules/glossary_units.md` 与 `architecture/modules/data_contract.md` 中确认
