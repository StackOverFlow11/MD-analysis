# md_analysis.enhanced_sampling.slowgrowth — Implementation Guidelines

## 职责

慢增长（slow-growth）热力学积分的数据处理、自由能积分、绘图和 CSV 导出。

## 文件组织

- `SlowGrowth.py`：数据类（`Slowgrowth`、`SlowgrowthFull`、`SlowgrowthSegment`）+ 积分函数 `_integrate_midpoint`
- `SlowGrowthPlot.py`：绘图（quick / publication）+ CSV 导出 + 统一入口 `slowgrowth_analysis`
- `config.py`：输出文件名常量

## 积分公式

CP2K 约定：自由能变为约束力（Lagrange 乘子）对集体变量的**负积分**：

$$
\Delta A_k = -\sum_{i=0}^{k-1} \frac{\lambda_i + \lambda_{i+1}}{2} \Delta\xi
$$

其中 $\Delta\xi$ = `target_growth_au`（常数步长，a.u./step），$\lambda_i$ = `lagrange_shake[i]`。

**单位转换**：CP2K restart 中的 `TARGET_GROWTH` 单位为 per a.u. time，需乘以 `dt_au = timestep_fs / AU_TIME_TO_FS` 转换为 per-step。此转换在 `SlowgrowthFull.from_md_info()` 中完成。

返回长度为 $n$ 的数组，`result[0] = 0`。

## 步范围逻辑

`slowgrowth_analysis` 中：
- `initial_step <= final_step`（或 `final_step=None`）：正向切片 `[initial, final)`
- `initial_step > final_step`：先取 `[final, initial)` 切片，再调用 `.reversed()` 反转（交换初末态、翻转数组、自由能翻转并重新归零）

## 反转逻辑（`Slowgrowth.reversed()`）

- `target_au`、`lagrange_shake`：翻转顺序
- `free_energy_au`：翻转 → 减去新起点值（重新归零）。归零操作自然产生符号反转：`ΔA(B→A) = -ΔA(A→B)`
- `target_growth_au`：取负
- `steps`：重置为 `[0, 1, ..., n-1]`
- `times_fs`：重置为 `[0, dt, 2dt, ...]`

## 溢出处理

CP2K 输出中 `***` 值由 `ColvarParser._safe_float()` 转换为 `np.nan`。含 NaN 的步在积分中传播（midpoint rule 中任一侧为 NaN 则该段积分为 NaN）。CLI 中 `_print_sg_info()` 检测并警告用户使用 `initial_step` / `final_step` 避开这些索引。

## 依赖

- `utils.RestartParser.ColvarParser`：`ColvarMDInfo`（解析 restart + log）
- `utils.config`：`HA_TO_EV`（Hartree → eV 转换）
- `utils._io_helpers`：`_write_csv`（CSV 写入）
- `matplotlib`：延迟导入（`import matplotlib; matplotlib.use("Agg")`）

## 绘图细节

### Quick plot

- 左轴：自由能 (eV)，红色
- 右轴：Lagrange 乘子 (a.u.)，蓝色
- 顶轴：MD 步编号（通过 `secondary_xaxis` 线性映射）
- 移动平均线：Lagrange 乘子的滑动均值（`ma_window` 默认 50），深蓝色，末帧标注 MA 值
- 标注：barrier peak（$\Delta F^\dagger$）+ 总自由能变（$\Delta F$），箭头指向对应点
- CV 递减时自动翻转 x 轴（`invert_xaxis`），确保时间从左到右

### Publication plot

- 左轴：自由能 (eV)，红色；右轴：Lagrange 乘子 (a.u.)，蓝色
- 使用 `rc_context` 设置 serif 字体 + STIX 数学字体
- 能量值显示在图例中（invisible handles），无箭头标注
- CV 递减时自动翻转 x 轴（`invert_xaxis`），确保时间从左到右

### CSV 输出

列：`step,time_fs,target_au,lagrange_au,free_energy_au,free_energy_ev`
