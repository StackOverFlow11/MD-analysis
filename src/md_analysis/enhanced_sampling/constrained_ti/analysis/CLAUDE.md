# constrained_ti/analysis — 开发备忘

## 定位

四步诊断引擎的纯函数实现。每个模块只依赖 `numpy`（无 scipy 依赖），接收 ndarray 返回 frozen dataclass。

## 模块

| 文件 | 输入 | 输出 | 用途 |
|---|---|---|---|
| `_acf_core.py` | series | acf, tau, n_eff, sem | ACF 计算核心（Sokal 1997 自洽截断） |
| `autocorrelation.py` | series, sem_max | `AutocorrResult` | Step 2: 封装 _acf_core + pass/fail |
| `block_average.py` | series, sem_max | `BlockAverageResult` | Step 3: F&P pow2 block averaging |
| `running_average.py` | series, sem | `RunningAverageResult` | Step 1: cumulative mean drift |
| `geweke.py` | series | `GewekeResult` | Step 4: Geweke z-test |

## block_average.py 关键设计

### F&P 方法（Flyvbjerg & Petersen 1989）

- Block sizes: 2 的幂次，`1, 2, 4, ..., B_max`，其中 `B_max = N // min_blocks`
- `SEM(B) = std(block_means, ddof=1) / √n_b`
- `δSEM(B) = SEM(B) / √(2(n_b − 1))`（来自 χ² 分布）
- 平台检测：连续 `n_consecutive` 个 pow2 level 的 SEM 增量 < √(δSEM_i² + δSEM_{i+1}²)（合成不确定度）
- `n_b < 2` 时 SEM 和 δSEM 均为 NaN（跳过）

### 公开函数

- `_generate_block_sizes(n, min_blocks)` → pow2 ndarray
- `analyze_block_average(series, *, sem_max, min_blocks, n_consecutive)` → `BlockAverageResult`

### 已删除（2026-03-24）

- `_generate_block_sizes()` dense 策略（1..20 连续 + 1.25x 几何级数）
- `_arctan_fit.py` 整个模块（arctan SEM 外推）
- 旧 plateau 检测（最后 window 个点相对展度法）
- `dense_sampling`, `arctan_r2_min`, `arctan_min_points` 等参数

## 约定

- 纯函数，无 logging，无 I/O
- 所有 dataclass 定义在 `../models.py`
- config 常量在 `../config.py`
