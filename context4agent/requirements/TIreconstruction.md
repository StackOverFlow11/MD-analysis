# Arctan 外推集成到 Block Average 模块 — 重构计划 (Rev.2, 已实施)

> **临时文件声明**：本文件（`TIreconstruction.md`）和同目录下的 `DESIGN.md` 均为临时规划文档，后续可能会清理。所有代码实现细节（模块架构设计、数据接口约定、算法实现、数据流转）在实施完成后必须同步到 `context4agent/` 对应节点目录和 `src/` 对应模块的 `CLAUDE.md` 中，以保证这些位置作为唯一持久的权威参考。

## 0. 背景与动机

### 0.1 发现

在对比 `shakeerror` 脚本（jiangch v2025.05）与我们的 `constrained_ti` 模块时发现：arctan 拟合 `SEM(B) = A·arctan(B·x)` 外推到 `A·π/2` 得到的渐近 SEM，与 FFT-ACF 方法（Sokal 1997 自洽截断）得到的 SEM 仅差 **0.4%**。

验证数据：`data_example/ti/ti_target_0.302356`（1970 帧，τ\_corr = 8.2）

| 方法 | SEM | 隐含 τ\_corr | N\_eff |
|------|-----|-------------|--------|
| shakeerror 原始公式（`std/(n-1)`，有 bug） | 2.10e-04 | 4.5 | 219 |
| 修正 SEM + arctan 外推 | 2.84e-04 | **8.2** | 120 |
| ACF 自洽截断 (Sokal) | 2.84e-04 | **8.2** | 120 |

### 0.2 当前 plateau 检测的问题

当前 `block_average.py` 使用末尾 4 个 2 的幂次块大小（`B = 2^7, 2^8, 2^9, 2^10`）的相对离散度 < 20% 来判断 plateau。在实际数据中经常失败：

- `ti_target_0.302356`：1970 帧，`plateau_reached = False`，被迫回退到 `SEM_auto`
- 原因：2 的幂次采样点太稀疏（仅 ~9 个点），末尾几个点波动大

### 0.3 arctan 外推的优势

1. **不要求曲线达到 plateau** — 从上升段就能外推渐近值
2. **定量估计** — 给出连续的 `A·π/2` 值，而非二值的 "plateau reached / not reached"
3. **隐含 τ\_corr** — 从渐近值直接反推，可与 ACF 交叉验证
4. **密集采样改善拟合** — 连续整数 + 几何级数的混合策略提供 ~30 个点

### 0.4 模型性质声明

> **重要**：arctan 模型 `SEM(B) = A·arctan(B·x)` 是**经验拟合函数**，不是 block averaging 的精确理论公式。精确公式（指数衰减 ACF 情况）为 `SEM(B) = σ_∞·√(1 − (τ/B)·(1−exp(−B/τ)))`。arctan 与之形状定性一致（从零起步、单调递增、趋于渐近值），但不等价。对于多指数衰减 ACF（非简单 AR(1)），SEM(B) 可能出现阶梯状上升，此时 arctan 拟合质量下降——R² 阈值机制会自动捕获此情况并 fallback。

### 0.5 目标

将 arctan 外推作为 block average 的**主要** SEM 估计方法，plateau 检测降为辅助诊断，ACF 作为最终 fallback。

---

## 1. 架构设计

### 1.1 SEM 选择三级层级（核心变更）

```
当前:  plateau → ACF (fallback)
新:    arctan 外推 → plateau → ACF (fallback)
```

```python
# workflow.py 中的决策逻辑
if block_avg.arctan is not None and block_avg.arctan.reliable:
    sem_final = block_avg.arctan.sem_asymptote   # 主选：arctan 渐近值
elif block_avg.plateau_reached:
    sem_final = block_avg.sem_plateau            # 备选：plateau 平台值
else:
    sem_final = autocorr.sem_auto                # 最终回退：ACF 方法
```

### 1.2 模块职责划分

```
analysis/_arctan_fit.py (新建)    — 纯拟合函数 fit_arctan_sem()，零状态依赖
models.py (修改)                  — 新增 ArctanFitResult；BlockAverageResult 新增 arctan 嵌套字段
analysis/block_average.py (修改)  — 密集采样 + 调用 _arctan_fit + 返回扩展结果
workflow.py (修改)                — 三级 SEM 层级 + _BLOCK_KEYS 更新
plot.py (修改)                    — 叠加拟合曲线
config.py (修改)                  — 新增默认阈值
__init__.py (修改)                — import + __all__ 新增 ArctanFitResult 导出
```

### 1.3 依赖关系

```
models.py (改)  ← config.py (改)
     ↑                ↑
_arctan_fit.py (新)   |        ← 只含纯函数，从 models 导入 ArctanFitResult
     ↑                |
block_average.py (改) ←-----   ← 调用 _arctan_fit，传入 sigma/n 用于 tau 计算
     ↑
workflow.py (改)               ← 三级 SEM 选择
     ↑
plot.py (改)                   ← 叠加拟合曲线（只从 models + config 导入类型）
```

- `_arctan_fit.py` **不定义任何 dataclass、不含 logging**（遵循 `_acf_core.py` 模式）
- `ArctanFitResult` 定义在 `models.py`（遵循所有结果类型集中定义的惯例）
- `plot.py` 只从 `models.py` 和 `config.py` 导入（保持"不导入 workflow/analysis"的干净模式）
- scipy.optimize.curve_fit 通过函数体内延迟导入处理（遵循 `calibration/_mapper.py` 模式）

### 1.4 scipy 可选依赖处理

scipy 不在 `pyproject.toml` 的 `dependencies` 中。处理策略：

1. `_arctan_fit.py` 中延迟导入，ImportError 时返回 `None`：
   ```python
   def fit_arctan_sem(...) -> ArctanFitResult | None:
       try:
           from scipy.optimize import curve_fit
       except ImportError:
           return None  # scipy 不可用 → 整体跳过
   ```
2. `block_average.py` 透传 `None` → `BlockAverageResult(arctan=None)`
3. `pyproject.toml` 的 `[project.optional-dependencies]` 新增：
   ```toml
   [project.optional-dependencies]
   ti = ["scipy>=1.10"]
   ```

**设计决策**：`fit_arctan_sem()` 返回 `ArctanFitResult | None`。`None` = 无法尝试（scipy 不可用/点数不足）；`ArctanFitResult(reliable=False)` = 已尝试但不可靠。这样在 workflow 层，`arctan is None` 和 `arctan.reliable == False` 有清晰的语义区分。

---

## 2. 详细修改清单

### 2.1 修改 `config.py` — 新增常量

```python
# ---------------------------------------------------------------------------
# Step 3b: Arctan extrapolation
# ---------------------------------------------------------------------------

# Minimum R² for arctan fit to be considered reliable.
DEFAULT_ARCTAN_R2_MIN: float = 0.95

# Minimum number of (non-NaN) block sizes needed to attempt arctan fit.
DEFAULT_ARCTAN_MIN_POINTS: int = 5
```

> **Rev.2 变更**：删除 `DEFAULT_ARCTAN_CROSS_VALID_RTOL`（Round 2 发现无任何代码路径引用此常量，为死代码）。

### 2.2 修改 `models.py` — 嵌套 dataclass 设计

**新增** `ArctanFitResult`（独立 frozen dataclass）：

```python
@dataclass(frozen=True)
class ArctanFitResult:
    """Result of arctan extrapolation on the SEM(B) curve."""
    sem_asymptote: float         # A · π/2 — 外推的渐近 SEM
    A: float                     # 幅度参数
    B: float                     # 速率参数
    r2: float                    # 拟合 R²
    reliable: bool               # R² ≥ 阈值 且 A>0, B>0 且 pcov 非奇异
    tau_corr_implied: float      # 从渐近值反推的 τ_corr = N·SEM²/(2·σ²)
    fit_curve: np.ndarray | None # 拟合曲线（与输入 block_sizes 等长），绘图用
```

**修改** `BlockAverageResult`（现有 8 字段不变，新增 1 个嵌套字段）：

```python
@dataclass(frozen=True)
class BlockAverageResult:
    # --- 现有字段（全部保留，顺序不变） ---
    block_sizes: np.ndarray
    sem_curve: np.ndarray
    sem_plateau: float
    sem_at_max_B: float
    plateau_rtol: float
    plateau_reached: bool
    cross_valid_ok: bool            # 语义不变：仍为 plateau vs ACF 交叉验证
    passed: bool | None             # 语义不变：仅反映 plateau 子步骤结论

    # --- 新增：arctan 外推结果（嵌套） ---
    arctan: ArctanFitResult | None  # None = 未尝试（点数不足/scipy 不可用）
```

**设计决策说明**：
- `cross_valid_ok` 保持原有语义（plateau SEM vs ACF SEM），**不复用**于 arctan
- `passed` 保持原有语义（仅反映 plateau），arctan 的 pass/fail 由 workflow 层聚合
- `arctan: None` vs `arctan.reliable=False` 语义区分明确（见 §1.4）
- 访问路径：`report.block_avg.arctan.sem_asymptote`、`report.block_avg.arctan.reliable`

### 2.3 新建 `analysis/_arctan_fit.py`

纯函数模块，与 `_acf_core.py` 同级。**不含任何 dataclass 定义，不含 logging**。

```python
def fit_arctan_sem(
    block_sizes: np.ndarray,
    sem_curve: np.ndarray,
    *,
    n_total: int,                                         # 序列总长度（用于 tau 计算）
    sigma_series: float,                                  # 序列标准差 ddof=0（用于 tau 计算）
    r2_min: float = DEFAULT_ARCTAN_R2_MIN,
    min_points: int = DEFAULT_ARCTAN_MIN_POINTS,
) -> ArctanFitResult | None:
    """Fit SEM(B) = A · arctan(B · x) via WLS and extrapolate to A · π/2.

    Returns None if scipy is unavailable or if fewer than min_points valid
    (non-NaN) data points exist. Returns ArctanFitResult(reliable=False)
    if the fit was attempted but failed quality checks.

    WLS weighting: sigma_i = sqrt(block_sizes[i] / n_total) passed to
    scipy.optimize.curve_fit with absolute_sigma=False. This downweights
    large-B points where fewer blocks yield noisier SEM estimates.
    """
```

> **Rev.2 变更**：
> - 返回类型改为 `ArctanFitResult | None`（`None` = 跳过，`reliable=False` = 已尝试但质量不达标）
> - 新增 `n_total` 和 `sigma_series` 参数，函数内部一次性计算 `tau_corr_implied`，避免 `dataclasses.replace`
> - WLS 权重统一为 `sigma_i = sqrt(block_sizes[i] / n_total)`（即 `1/sqrt(n_blocks_i)`），文档注明这是简化近似
> - 删除"近乎水平检测"（Round 2 发现阈值 1.1 在实际 iid 数据中永不触发，是死代码；R² < 0.95 已能捕获此情况）

实现细节：
- 延迟导入 `scipy.optimize.curve_fit`，`ImportError` → 返回 `None`
- **预处理**：过滤 NaN 值（`mask = ~np.isnan(sem_curve)`），有效点数 < `min_points` → 返回 `None`
- **加权最小二乘**：`sigma=np.sqrt(block_sizes[mask] / n_total)`，`absolute_sigma=False`
  - 注释说明：此权重为简化近似；精确权重应为 `SEM(B)/sqrt(2*(n_blocks-1))`，但在 `absolute_sigma=False` 模式下只影响相对权重比，实际差异很小
- 参数边界 `bounds=([1e-15, 1e-15], [inf, inf])`（避免精确零导致 Jacobian 退化）
- 初始猜测：`A0 = 2·max(sem_valid)/π`，`B0 = 1/median(block_sizes_valid)`
- **收敛检查**：捕获 `RuntimeError`（不收敛）→ `ArctanFitResult(reliable=False, sem_asymptote=nan, ...)`
- **pcov 检查**：`np.any(np.isinf(np.diag(pcov)))` → `reliable=False`
- **可靠性判据**：`r2 >= r2_min` AND `A > 1e-15` AND `B > 1e-15` AND pcov 非奇异
- **tau_corr_implied 计算**（在此函数内完成）：
  ```python
  if reliable and sigma_series > 1e-30:
      tau_corr_implied = n_total * sem_asymptote**2 / (2 * sigma_series**2)
  else:
      tau_corr_implied = float("nan")
  ```

### 2.4 修改 `analysis/block_average.py` — 核心变更

**函数签名变更**：

```python
def analyze_block_average(
    series: np.ndarray,
    *,
    sem_max: float | None = None,
    sem_auto: float | None = None,
    min_blocks: int = DEFAULT_MIN_BLOCKS,
    plateau_window: int = DEFAULT_PLATEAU_WINDOW,
    plateau_rtol: float = DEFAULT_PLATEAU_RTOL,
    cross_valid_rtol: float = DEFAULT_CROSS_VALID_RTOL,
    dense_sampling: bool = True,        # 新增：True=混合策略，False=传统 powers-of-2
    arctan_r2_min: float = DEFAULT_ARCTAN_R2_MIN,         # 新增
    arctan_min_points: int = DEFAULT_ARCTAN_MIN_POINTS,   # 新增
) -> BlockAverageResult
```

**块大小采样策略**（新增 `_generate_block_sizes()` 私有函数）：

```python
def _generate_block_sizes(n: int, min_blocks: int = 4) -> np.ndarray:
    """混合策略：小块用连续整数，大块用几何级数。

    B = 1..20 连续整数 + 从 21 开始的 1.25x 几何级数。
    覆盖 B ≈ 2τ 的典型转折区域（MD 数据 τ 通常为 5~15）。
    """
    max_b = n // min_blocks if min_blocks > 0 else n
    if max_b < 1:
        return np.array([1], dtype=int)
    sizes = set(range(1, min(21, max_b + 1)))  # B = 1..20 连续
    b = 21.0                                    # 几何级数从 21 开始
    while int(b) <= max_b:
        sizes.add(int(b))
        b *= 1.25
    return np.array(sorted(sizes), dtype=int)
```

> **Rev.2 变更**：
> - 连续整数范围从 `1..16` 扩展到 `1..20`，几何级数从 21 开始（消除 16→20 间隙，额外仅 4 个点）
> - 添加 `max_b < 1` 保底逻辑（修复 `n < 4` 时返回空数组的 bug）

**主函数逻辑**：
1. 根据 `dense_sampling` 选择块大小生成策略（`True` → `_generate_block_sizes()`，`False` → 传统 powers-of-2）
2. 计算 SEM 曲线（逻辑不变）
3. plateau 检测（逻辑不变）
4. 交叉验证 `cross_valid_ok`（逻辑不变，仍比较 `sem_plateau` vs `sem_auto`）
5. **新增**：调用 arctan 拟合
   ```python
   sigma = float(np.std(series, ddof=0))  # 与 _acf_core 保持一致
   arctan_result = fit_arctan_sem(
       block_sizes, sem_curve,
       n_total=len(series), sigma_series=sigma,
       r2_min=arctan_r2_min, min_points=arctan_min_points,
   )
   ```
6. 返回 `BlockAverageResult(..., arctan=arctan_result)`

> **Rev.2 变更**：
> - 删除 `dataclasses.replace` 模式，`tau_corr_implied` 在 `fit_arctan_sem()` 内部一次性完成
> - `sigma` 和 `n_total` 作为参数传入 `fit_arctan_sem()`

### 2.5 修改 `workflow.py`

1. **SEM 选择从两级改为三级**（见 §1.1）
2. **failure_reasons 更新**：
   - arctan 为 None 或不可靠时：`"Arctan fit unavailable/unreliable (R²={r2:.3f}); trying plateau."`
   - plateau 也未达到时：`"Block-average plateau not reached; falling back to SEM_auto."`
3. **`_POINT_CSV_COLUMNS` 新增列**（追加在末尾，不插入中间）：
   ```python
   "sem_arctan",           # block_avg.arctan.sem_asymptote
   "arctan_r2",            # block_avg.arctan.r2
   "arctan_reliable",      # block_avg.arctan.reliable
   "tau_corr_implied",     # block_avg.arctan.tau_corr_implied
   ```
4. **`_point_to_row()` 新增对应字段**：
   ```python
   "sem_arctan": f"{r.block_avg.arctan.sem_asymptote:.8f}" if r.block_avg.arctan else "N/A",
   "arctan_r2": f"{r.block_avg.arctan.r2:.4f}" if r.block_avg.arctan else "N/A",
   "arctan_reliable": str(r.block_avg.arctan.reliable) if r.block_avg.arctan else "N/A",
   "tau_corr_implied": f"{r.block_avg.arctan.tau_corr_implied:.2f}" if r.block_avg.arctan else "N/A",
   ```
5. **`_BLOCK_KEYS` 更新**（新增 3 个键）：
   ```python
   _BLOCK_KEYS = {
       "min_blocks", "plateau_window", "plateau_rtol", "cross_valid_rtol",
       "dense_sampling", "arctan_r2_min", "arctan_min_points",  # 新增
   }
   ```

> **Rev.2 变更**：新增第 5 点 `_BLOCK_KEYS` 更新（Round 2 发现遗漏，不更新会导致 engine overrides 中的新参数被静默丢弃）。

### 2.6 修改 `plot.py` — 底部左图增强

在 SEM(B) 散点图上叠加：
- x 轴刻度：当 `block_avg.arctan is not None` 时使用 `log(base=10)`，否则保留 `log(base=2)`（兼容 `dense_sampling=False` 旧模式）
- arctan 拟合曲线（密集插值 x 的平滑红线，仅当 `arctan is not None` 时绘制）
- `sem_asymptote = A·π/2` 水平虚线（可靠=绿色 `C2`，不可靠=灰色 `0.6`）

> **Rev.2 变更**：x 轴 base 动态选择（Round 2 发现固定 log10 在 `dense_sampling=False` 模式下降低 powers-of-2 数据可读性）。

右下 summary table 新增：
- `SEM_arctan = xxx (R²=0.xx)`（arctan 存在时）
- `τ_implied = xxx`（arctan 可靠时）

标题更新：
```python
arctan_ok = block_avg.arctan is not None and block_avg.arctan.reliable
ax.set_title(
    f"Block Average (arctan={'yes' if arctan_ok else 'no'}, "
    f"plateau={'yes' if report.block_avg.plateau_reached else 'no'})"
)
```

### 2.7 修改 `__init__.py` — import 和导出更新

```python
from .models import (
    ArctanFitResult,          # 新增
    AutocorrResult,
    BlockAverageResult,
    ...
)

__all__ = [
    # Models
    "ArctanFitResult",          # 新增
    "AutocorrResult",
    "BlockAverageResult",
    "ConstraintPointInput",
    "ConstraintPointReport",
    "GewekeResult",
    "RunningAverageResult",
    "TIPointDefinition",
    "TIReport",
    # Exceptions
    "ConvergenceError",
    "InsufficientSamplingError",
]
```

> **Rev.2 变更**：明确列出 import 语句更新（Round 2 发现 Rev.1 只提了 `__all__` 未提 import）。

---

## 3. 边界情况处理

| 场景 | 处理 | `fit_arctan_sem` 返回值 |
|------|------|------------------------|
| scipy 不可用 | ImportError → 整体跳过 | `None` |
| 有效点数（非 NaN）< 5 | 跳过拟合 | `None` |
| `curve_fit` 不收敛 | 捕获 RuntimeError | `ArctanFitResult(reliable=False, sem_asymptote=nan, ...)` |
| pcov 包含 inf（奇异协方差矩阵） | 质量不达标 | `ArctanFitResult(reliable=False, ...)` |
| R² < 0.95 | 质量不达标 | `ArctanFitResult(reliable=False, ...)` |
| A ≈ 0 或 B ≈ 0（收敛到边界） | 质量不达标 | `ArctanFitResult(reliable=False, ...)` |
| SEM 曲线单调递减 | R² 很低，自然 fallback | `ArctanFitResult(reliable=False, ...)` |
| SEM 曲线含 NaN | 拟合前过滤，有效点 < 5 → 跳过 | `None` 或 `ArctanFitResult` |
| `sigma_series ≈ 0`（常数序列） | `tau_corr_implied = NaN` | `ArctanFitResult(tau_corr_implied=nan)` |
| 小 N（< 4）| `_generate_block_sizes` 返回 `[1]`，1 个点 < 5 → 跳过 | `None` |

> **Rev.2 变更**：
> - 统一了 `None` vs `ArctanFitResult(reliable=False)` 的语义（§1.4 设计决策）
> - 删除"近乎水平检测 max/min < 1.1"（Round 2 发现是死代码，R² 机制已覆盖）
> - 新增第三列，明确每种场景的返回值

---

## 4. 测试计划

### 4.1 新建 `test/unit/enhanced_sampling/constrained_ti/test_arctan_fit.py`

| # | 测试 | 断言 | 优先级 |
|---|------|------|--------|
| 1 | 合成 arctan 曲线 + 小噪声（seed=42） | 恢复 A, B 参数（rtol < 5%），R² > 0.99，`reliable=True` | P0 |
| 2 | AR(1)（N=10000, φ=0.9, seed=42）的真实 SEM 曲线 | `sem_asymptote` 与理论 SEM 偏差 < **15%** | P0 |
| 3 | 全零 SEM 曲线 | 返回 `ArctanFitResult(reliable=False)` | P0 |
| 4 | 单调递减数据 | 返回 `ArctanFitResult(reliable=False)` | P0 |
| 5 | 有效点数 < 5 | 返回 `None` | P0 |
| 6 | scipy 不可用（`unittest.mock.patch` 模拟 ImportError） | 返回 `None` | P0 |
| 7 | `curve_fit` 抛出 RuntimeError（mock `curve_fit`） | 返回 `ArctanFitResult(reliable=False)`，无异常泄漏 | P0 |
| 8 | SEM 曲线含 NaN 值（中间 3 个点为 NaN） | 过滤后拟合正常或返回 `None`（视剩余点数） | P1 |
| 9 | 极小 SEM 值（~1e-10 量级） | 数值稳定，不 overflow/underflow | P1 |
| 10 | R² 恰在阈值（0.949 → `reliable=False`; 0.951 → `reliable=True`） | 阈值判断无 off-by-one | P2 |
| 11 | pcov 对角线含 inf（mock `curve_fit` 返回 inf pcov） | `reliable=False` | P1 |
| 12 | `tau_corr_implied` 与已知理论值对比（合成数据） | 偏差 < 15% | P1 |
| 13 | `sigma_series=0` | `tau_corr_implied=nan`，不 crash | P1 |

> **Rev.2 变更**：
> - AR(1) 容差从 10% 放宽到 **15%**（Round 2 发现 10% 在统计波动下可能导致 flaky test）
> - 新增 #11（pcov inf）、#12（tau 回归）、#13（sigma=0 保护）
> - 明确 #5 和 #6 返回 `None`，#3/#4/#7 返回 `ArctanFitResult(reliable=False)`

### 4.2 新建 `test/unit/enhanced_sampling/constrained_ti/test_block_average_dense.py`

| # | 测试 | 断言 | 优先级 |
|---|------|------|--------|
| 1 | `_generate_block_sizes(2000)` | 包含 1..20，总数 ~28-38，严格递增 | P0 |
| 2 | `_generate_block_sizes(40)` | 正确处理小 N，max_b=10，返回 `[1..10]` | P0 |
| 3 | `_generate_block_sizes(4)` | max_b=1，返回 `[1]` | P0 |
| 4 | `_generate_block_sizes(0)` | 返回 `[1]`，不 crash | P0 |
| 5 | `_generate_block_sizes(3)` | max_b=0 → 返回 `[1]` | P0 |
| 6 | AR(1)（N=5000, φ=0.9, seed=42）密集采样 | `arctan.reliable=True`，`sem_asymptote` 与理论一致（rtol < 15%） | P0 |
| 7 | `dense_sampling=False` 传统模式 | `block_sizes` 为 2 的幂次，`arctan` 仍被计算 | P0 |
| 8 | `tau_corr_implied` 与 ACF τ\_corr 对比 | 偏差 < 25% | P1 |
| 9 | sigma ≈ 0 的常数序列 | `arctan` 为 `None` 或 `tau_corr_implied=nan`，不 crash | P1 |

> **Rev.2 变更**：
> - 新增 #4 和 #5（Round 2 发现的边界 bug 的测试）
> - tau_corr 对比容差放宽到 25%（间接计算，误差传播约 2x）
> - #7 明确 `dense_sampling=False` 下 arctan 仍然被计算（只是 block_sizes 不同）

### 4.3 修改 `test/unit/enhanced_sampling/constrained_ti/test_workflow.py`

| # | 测试 | 断言 | 构造方法 | 优先级 |
|---|------|------|---------|--------|
| 1 | `report.block_avg.arctan` 字段存在 | `arctan is not None`（AR(1) N=5000） | 真实 `analyze_standalone` | P0 |
| 2 | 三级路径 1：arctan 可靠 → `sem_final = sem_asymptote` | `sem_final ≈ arctan.sem_asymptote` | 真实 AR(1) N=5000 φ=0.9 | P0 |
| 3 | 三级路径 2：arctan 不可靠 + plateau → `sem_final = sem_plateau` | `sem_final ≈ sem_plateau` | **手动构造** mock `ConstraintPointReport`（`arctan=ArctanFitResult(reliable=False, ...)`，`plateau_reached=True`）| P0 |
| 4 | 三级路径 3：arctan 不可靠 + no plateau → `sem_final = sem_auto` | `sem_final ≈ sem_auto` | **手动构造** mock（`arctan=None`，`plateau_reached=False`）| P0 |
| 5 | failure_reasons 包含 "Arctan fit" | 路径 2/3 触发时检查 | mock 构造 | P0 |
| 6 | `test_clean_ar1_all_pass` 回归 | `all_passed` 不变 + SEM 量级不变（rtol < 20%） | 原测试 + 新增断言 | **P0** |
| 7 | CSV 新增 4 列存在且格式正确 | 解析输出 CSV 验证列名 | 真实运行 | P1 |
| 8 | `_BLOCK_KEYS` dispatch 测试 | `analyze_standalone(..., arctan_r2_min=0.99)` 确实传递到 block_average | 检查 `report.block_avg.arctan.r2` 阈值行为 | P1 |

> **Rev.2 变更**：
> - 路径 2 明确使用**手动构造 mock**（Round 2 指出该路径无法用真实数据可靠触发）
> - `test_clean_ar1_all_pass` 回归从 P1 升为 **P0**（密集采样默认开启，是行为变更）
> - 新增 #5 failure_reasons 内容检查（P0）
> - 新增 #8 `_BLOCK_KEYS` dispatch 验证

### 4.4 新建 `test/integration/enhanced_sampling/test_constrained_ti_regression.py`

| # | 测试 | 断言 | 优先级 |
|---|------|------|--------|
| 1 | `data_example/ti/ti_target_0.302356` 回归 | `arctan.reliable=True`，`sem_asymptote` vs `sem_auto` 偏差 < 5%，`sem_auto ≈ 2.84e-04`（rtol < 5%） | P0 |
| 2 | 诊断图生成（arctan 存在） | PNG 文件存在且不报错 | P1 |
| 3 | 诊断图生成（mock scipy 不可用，arctan=None） | PNG 文件存在且不报错（plot.py 不因 arctan=None crash） | P1 |

> **Rev.2 变更**：新增 #3 smoke test（Round 2 发现 plot.py 在 `arctan=None` 时的行为未被测试）。

### 4.5 数值精度约定

- 所有 `pytest.approx` 明确使用 `rel=` 参数（而非默认绝对容差）
- 随机种子固定（seed=42 等），确保可复现
- AR(1) 测试使用 N ≥ 5000 以减少统计波动
- 合成数据测试 rtol < 5%；真实/AR(1) 数据 rtol < 15%；间接量（tau）rtol < 25%

---

## 5. 实施顺序

| Phase | 文件 | 内容 |
|-------|------|------|
| 1 | `config.py` | 新增 2 个常量（R²阈值、最小点数） |
| 2 | `models.py` | 新增 `ArctanFitResult` dataclass；`BlockAverageResult` 新增 `arctan: ArctanFitResult \| None` |
| 3 | `_arctan_fit.py` | 全新文件：`fit_arctan_sem()` 纯函数（WLS + NaN 过滤 + pcov 检查 + tau 计算），返回 `ArctanFitResult \| None` |
| 4 | `block_average.py` | `_generate_block_sizes()` + `dense_sampling` 开关 + 调用 `fit_arctan_sem()` 传入 `sigma`/`n_total` |
| 5 | `workflow.py` | 三级 SEM 层级 + failure_reasons + CSV 列（末尾追加）+ `_BLOCK_KEYS` 更新 |
| 6 | `plot.py` | x 轴动态 base + 拟合曲线叠加 + 渐近线 + 标题更新 |
| 7 | `__init__.py` | import + `__all__` 新增 `ArctanFitResult` |
| 8 | `pyproject.toml` | `[project.optional-dependencies]` 新增 `ti = ["scipy>=1.10"]` |
| 9 | 测试文件 | test_arctan_fit + test_block_average_dense + test_workflow 修改 + 集成回归测试 |
| 10 | `context4agent/` | 同步更新 architecture 文档 + DESIGN.md |

---

## 6. 验证清单

1. `pip install .` → `pytest test/unit/enhanced_sampling/constrained_ti/ -v` 全部通过
2. 用 `data_example/ti/ti_target_0.302356` 跑单点诊断：
   - `arctan.reliable = True`
   - `arctan.sem_asymptote` 与 `autocorr.sem_auto` 偏差 < 5%
   - 诊断图显示 arctan 拟合曲线（红线）和渐近线（绿色虚线）
3. `dense_sampling=False` 模式下，现有测试行为不变（向后兼容验证）
4. 全量回归：`pytest test/unit/ -v` + `pytest test/integration/ -v` 无新增失败
5. 使用高 τ/N 比率数据验证 arctan 在困难场景下的行为（手动/半自动）

---

## 附录 A：Rev.2 修订记录

| 来源 | 问题 | 修正 |
|------|------|------|
| R2-架构-A | `dataclasses.replace` 无先例 | `fit_arctan_sem` 新增 `n_total`/`sigma_series` 参数，内部计算 `tau_corr_implied` |
| R2-架构-B | `dense_sampling=True` 隐式变更行为 | `test_clean_ar1_all_pass` 回归升为 P0 |
| R2-架构-C | `_generate_block_sizes(n<4)` 返回空数组 | 添加 `max_b < 1` 保底返回 `[1]` |
| R2-架构-D | scipy 不可用时返回值矛盾 | 统一：`fit_arctan_sem` 返回 `None`；`reliable=False` 仅用于已尝试但质量不达标 |
| R2-架构-E | WLS 权重文字不一致 | 统一为 `sigma_i = sqrt(block_sizes[i]/n_total)` + 注释说明 |
| R2-架构-F | `__init__.py` import 未更新 | 明确列出 import 语句 |
| R2-架构-G | `_BLOCK_KEYS` 未更新 | §2.5 新增第 5 点 |
| R2-架构-H | log base 固定 log10 影响旧模式 | 动态选择 base |
| R2-架构-I | `DEFAULT_ARCTAN_CROSS_VALID_RTOL` 死常量 | 删除 |
| R2-算法-2 | 16→20 间隙修复不充分 | 连续整数扩展到 1..20，几何从 21 开始 |
| R2-算法-6 | 近乎水平检测 1.1 是死代码 | 删除，依赖 R² 机制 |
| R2-算法-7 | `DEFAULT_ARCTAN_CROSS_VALID_RTOL` 未使用 | 删除 |
| R2-测试-1 | scipy mock 层级不匹配 | 明确 `fit_arctan_sem` 返回 `None`，block_average 透传 |
| R2-测试-2 | 路径 2 无法用真实数据触发 | 明确使用手动构造 mock |
| R2-测试-3 | `_BLOCK_KEYS` dispatch 无测试 | §4.3 新增 #8 |
| R2-测试-6 | AR(1) 容差 10% 偏紧 | 放宽到 15% |
| R2-测试-7 | plot.py arctan=None 无 smoke test | §4.4 新增 #3 |

## 附录 B：Round 3 实施建议（非阻塞）

以下 2 条建议来自 Round 3 测试审核，属于实施细节级别澄清，不影响计划结构：

**(a) §4.3 #3/#4 的 mock 层级**：路径 2/3 的测试应 patch `analyze_block_average` 的**返回值**（让它返回一个包含 `arctan=ArctanFitResult(reliable=False, ...)` 的 `BlockAverageResult`），然后通过 `analyze_standalone` 端到端调用验证 `report.sem_final` 的选择。不应直接构造 `ConstraintPointReport`（它是输出而非输入）。

**(b) §4.1 #7 的 mock target path**：由于 `curve_fit` 是在 `fit_arctan_sem()` 函数体内延迟导入的（`from scipy.optimize import curve_fit`），mock target 应为 `scipy.optimize.curve_fit`（而非 `_arctan_fit.curve_fit`），或通过 `builtins.__import__` 控制延迟导入行为。具体方案在实施时确定。
