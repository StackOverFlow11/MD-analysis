# slowgrowth — 开发备忘

## 定位

慢增长（slow-growth）自由能计算：数据结构、midpoint 积分、双轴绘图、CSV 导出。

## 约定

### CP2K 符号约定（关键）
- 自由能积分取**负号**：ΔA_k = −Σ (λ_i + λ_{i+1})/2 × Δξ
- `target_growth_au` 在 `Slowgrowth` 中是 **per step**（已从 restart 的 per-a.u.-time 转换）
- Δξ = `target_growth_au`（per-step 增量）

### 类继承链
- `Slowgrowth`（frozen dataclass，基类）→ `SlowgrowthFull`（`from_paths()`/`from_md_info()`，保留 `md_info`）→ `SlowgrowthSegment`（`segment()` 截取，re-zero FE，保留 `parent` 引用）
- `reversed()`：反转数组 + 取反 growth + 取反并翻转 FE + 重置 step/time

### 绘图
- Quick plot：双轴，左=FE (eV)，右=Lagrange 乘子 (a.u.)，含 moving average
- Publication plot：增强格式化
- Bottom x-axis=CV (a.u.)，top x-axis=MD step（仅 quick plot）
- CV 递减时自动反转 x 轴

## 陷阱与历史 Bug

- **Bug bd3ed98**：`reversed()` 曾多一次取反 FE，导致符号错误。正确：negate → flip → re-zero（减去首元素）
- **Bug 7c69901**：TARGET_GROWTH per-a.u.-time → per-step 转换遗漏，导致积分偏差一个 dt_au 因子
- **Bug e3e116e**：积分曾用正号，应为负号（CP2K 约定）
- **Bug bc16b9a**：CV 递减时 x 轴未反转，视觉上自由能曲线方向反了
- `SlowgrowthFull.md_info` 保留原始 `ColvarMDInfo` — 可用于回溯 restart 参数
- NaN 步（overflow）在积分中传播 — 下游绘图需处理 NaN gap
