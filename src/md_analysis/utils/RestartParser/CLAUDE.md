# RestartParser — 开发备忘

## 定位

CP2K restart 文件和 LagrangeMultLog 日志的解析器。

## 约定

- **CellParser**：支持 `.restart`（正则提取 `&CELL` 块的 A/B/C 向量）和 `md.inp`（正则提取 `ABC [angstrom]` 行）
- **ColvarParser**：所有正则匹配为 case-insensitive
- **数据类**：全部 `frozen=True`（`ConstraintInfo`、`ColvarInfo`、`ColvarRestart`、`LagrangeMultLog`、`ColvarMDInfo`）
- **ColvarInfo 索引**：`colvars[colvar_id]` 按 `colvar_id` 查找（不是列表位置），`.primary` 返回第一个约束

## 陷阱与历史 Bug

### TARGET_GROWTH 单位（关键！— bug 7c69901）
- restart 中的 `TARGET_GROWTH` 是 **per a.u. time**，不是 per step
- 转换为 per-step：`growth_per_step = target_growth_au × (timestep_fs / AU_TIME_TO_FS)`
- 不做此转换会导致自由能积分偏差一个时间步因子

### Target Series 公式
- `ξ(k) = target_au + (k - step_start) × target_growth_au × dt_au`
- `target_au` 是 restart 快照时刻（`step_start`）的值，不是初始值
- `k` 是绝对步数（从 MD 开始计），不是相对偏移

### LagrangeMultLog Overflow（bug 58d0b58）
- CP2K 在 SHAKE/RATTLE 乘子过大时输出 `***`（星号串）
- `_safe_float()` 将 `***` 解析为 `float('nan')`，不报错
- NaN 向下游传播，允许绘图/分析优雅降级

### Fixed Atoms LIST 语法
- 支持 `\` 续行、逗号/空格分隔、`N..M` 范围展开
- 返回 1-indexed 原子编号的 sorted tuple（保持 CP2K 约定）

### 多 COLLECTIVE 块
- 一个 `&CONSTRAINT` 内可有多个 `&COLLECTIVE` 块（耦合 CV 场景）
- 每个块有独立的 `colvar_id`、`TARGET`、`TARGET_GROWTH`
