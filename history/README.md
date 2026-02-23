# History / Context

本目录用于**记录与项目持续相关的上下文**，供后续开发/Agent 协作时快速对齐，避免因对话过长导致遗忘。

## 使用规则（重要）

- **近期诉求（Short-term）持续更新**：每次提出"最近要做什么/优先级变化/本周目标"，都更新 `requirements/short_term.md`。
- **远期诉求（Long-term）按需更新**：只有明确要求"更新远期诉求/路线图"时，才更新 `requirements/long_term.md`。
- **临时重要上下文先放 temp**：不确定归属时先写入 `temp/`，后续再搬运并在原条目标记"已迁移"。
- **新增文件先加索引**：在本 README 的"目录索引"里补一条，方便检索。
- **未协商不得固化约定**：任何**未与用户协商确定**的接口/口径/命名/数据契约/单位约定，**不得自行采用或写入** `history/` 对应部分；只能先记录为"待讨论/待确认"。

## 目录索引

- `user_level.md`：用户当前水平、习惯、偏好（面向协作效率）
- `requirements/`
  - `short_term.md`：近期诉求（持续更新）
  - `long_term.md`：远期诉求（按需更新）
  - `README.md`：诉求拆分原则与维护说明
- `architecture/`
  - `README.md`：代码架构总览（目录/模块职责/数据流）
  - `modules/`：子模块说明（含 `scripts/` 目录层级对应的接口暴露/实现准则文档）
    - `modules/README.md`：模块文档治理硬约束（目录镜像与双文档）
    - `modules/data_contract.md`：跨模块核心数据契约（形状/单位/CSV 列头）
    - `modules/glossary_units.md`：术语与单位约定汇总
  - `decisions/`：架构/技术决策记录（ADR）
- `temp/`：重要上下文临时寄存区

## 约定

- Markdown 公式格式：行内使用 `$...$`，块级使用 `$$...$$`；不使用 `\(...\)` 与 `\[...\]` 作为公式定界符
- 术语与单位制：已在 `architecture/modules/glossary_units.md` 与 `architecture/modules/data_contract.md` 中确认
