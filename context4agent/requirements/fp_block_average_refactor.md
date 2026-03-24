# 重构：用 Flyvbjerg & Petersen 方法替换 block average 模块

## 背景

当前 block average 模块 (`analysis/block_average.py` + `analysis/_arctan_fit.py`) 过于复杂：

- dense 采样策略 (1..20 连续 + 1.25x 几何级数)
- 旧 plateau 检测（最后 window 个点的相对展度）
- arctan 曲线拟合 (`SEM(B) = A·arctan(B·x)`) 作为主 SEM 估计
- 3 层 fallback: arctan → plateau → ACF

实测问题：大 B 端 SEM(B) 噪声导致 arctan 拟合不稳定（R² 经常 < 0.75），5 个真实数据点中 3 个 arctan 判为 unreliable。

经测试验证，经典 F&P (1989) 方法简洁鲁棒：
- pow2 block sizes
- δSEM(B) = SEM(B) / √(2(n_b-1)) 作为每个点的不确定度
- 连续 level 增量 < δSEM 即为平台
- 无需曲线拟合

5 个真实 TI 点上 F&P SEM 与 ACF SEM_auto 一致（~10% 内）。

## 审核结论（2026-03-24）

### 假设一般性
- δSEM 公式在 B < τ 时不严格（blocks 间有相关），但不影响平台检测
- pow2 间隔对标准 TI 数据足够（τ ≈ 5-50）
- 极长 τ（> N/16）由 ACF fallback 兜底

### 数值稳定性
- 常数序列、NaN、极短序列均安全（上游 validate_and_trim 保护）
- n_blocks < 2 时应直接跳过（避免 ddof=1 除零）

### 出图 + 架构
- 高亮 plateau 窗口点以传达判定依据
- 存 `n_total` 而非 `n_blocks` 数组（减少冗余）
- CSV 直接删旧列 + 加新列（0.1.x 阶段 breaking change 可接受）

## 删除清单

| 文件 | 删除内容 |
|------|---------|
| `models.py` | `ArctanFitResult` dataclass |
| `models.py` | `BlockAverageResult` 旧字段：`sem_plateau`, `sem_at_max_B`, `plateau_rtol`, `cross_valid_ok`, `arctan` |
| `analysis/block_average.py` | `_generate_block_sizes()` (dense 策略)、旧 plateau 逻辑 |
| `analysis/_arctan_fit.py` | **整个文件** |
| `config.py` | `DEFAULT_PLATEAU_RTOL`, `DEFAULT_PLATEAU_WINDOW`, `DEFAULT_CROSS_VALID_RTOL`, `DEFAULT_ARCTAN_R2_MIN`, `DEFAULT_ARCTAN_MIN_POINTS`, `DEFAULT_MIN_BLOCKS` |
| `workflow.py` | arctan 相关 sem_final 分支、cross_valid 检查、CSV 列 (`sem_arctan`, `arctan_r2`, `arctan_reliable`, `tau_corr_implied`, `sem_block`, `sem_at_max_B`) |
| `plot.py` | arctan 拟合曲线、渐近线绘制 |
| `__init__.py` | `ArctanFitResult` 导出 |
| `test_arctan_fit.py` | **整个文件** |
| `test_block_average_dense.py` | **整个文件**（替换为 `test_block_average.py`） |

## 实施步骤

### Step 1: `config.py` — 常量替换

删除旧常量，新增：
```python
DEFAULT_FP_MIN_BLOCKS: int = 4       # pow2 最大 B 处的最小 block 数
DEFAULT_FP_CONSECUTIVE: int = 2      # 连续几个 level 不显著增长视为平台
DEFAULT_CROSS_CHECK_RTOL: float = 0.15  # SEM_block vs SEM_auto 交叉验证容差
```

### Step 2: `models.py` — 数据模型

删除 `ArctanFitResult`。`BlockAverageResult` 重写为：

```python
@dataclass(frozen=True)
class BlockAverageResult:
    """Output of Step 3: Flyvbjerg-Petersen block averaging."""
    block_sizes: np.ndarray      # pow2 block sizes
    sem_curve: np.ndarray        # SEM(B) at each level
    delta_sem: np.ndarray        # δSEM(B) = SEM(B) / √(2(n_b-1))
    n_total: int                 # series length (for downstream n_b = n_total // B)
    plateau_index: int | None    # 首次检测到平台的索引
    plateau_sem: float           # 平台处 SEM (称为 SEM_block)
    plateau_delta: float         # 平台处 δSEM
    plateau_block_size: int | None
    plateau_reached: bool
    passed: bool | None          # plateau_reached AND SEM <= SEM_max
```

### Step 3: `analysis/block_average.py` — 核心算法

删除 `_generate_block_sizes()`（dense）、旧 plateau 逻辑、arctan import。

保留 `_generate_block_sizes_pow2()` → 重命名为 `_generate_block_sizes()`。

`analyze_block_average()` 简化：

```python
def analyze_block_average(
    series: np.ndarray,
    *,
    sem_max: float | None = None,
    min_blocks: int = DEFAULT_FP_MIN_BLOCKS,
    n_consecutive: int = DEFAULT_FP_CONSECUTIVE,
) -> BlockAverageResult:
```

算法：
1. 生成 pow2 block sizes
2. 对每个 B：n_b = N // B；**n_b < 2 时跳过**（避免 ddof=1 除零）
3. SEM(B) = std(blocks, ddof=1) / √n_b
4. δSEM(B) = SEM(B) / √(2(n_b-1))
5. 平台检测：遍历相邻 level，若 SEM(B_{i+1}) - SEM(B_i) < δSEM(B_i) 连续 `n_consecutive` 次 → plateau 在首次满足处
6. plateau_sem = SEM(plateau_index)；未达平台则取最大 B 处

### Step 4: `analysis/_arctan_fit.py` — 删除

### Step 5: `workflow.py`

**`_BLOCK_KEYS`**：`{"min_blocks", "n_consecutive"}`

**sem_final** 简化为 2 层：
```python
if block_avg.plateau_reached:
    sem_final = block_avg.plateau_sem
else:
    sem_final = autocorr.sem_auto
    failure_reasons.append("F&P plateau not reached; falling back to SEM_auto.")
```

**交叉验证**：SEM_block vs SEM_auto 差距 > 15% 时发出 warning：
```python
if block_avg.plateau_reached:
    rel_diff = abs(block_avg.plateau_sem - autocorr.sem_auto) / max(block_avg.plateau_sem, autocorr.sem_auto)
    if rel_diff > DEFAULT_CROSS_CHECK_RTOL:
        failure_reasons.append(
            f"SEM_block ({block_avg.plateau_sem:.2e}) and SEM_auto ({autocorr.sem_auto:.2e}) disagree by {rel_diff:.0%}."
        )
```

**`analyze_block_average()` 调用**：不再传 `sem_auto`。

**CSV 列**（直接删旧列 + 加新列）：
```
xi, lambda_mean, sigma_lambda, tau_corr, n_eff, sem_auto,
sem_block, delta_sem_block, plateau_B, plateau_reached,
sem_final, sem_final_method, sem_max,
geweke_z, geweke_reliable, drift_D,
passed, failure_reasons
```

`sem_final_method` 列标明来源（`"plateau"` / `"acf"`）。

### Step 6: `plot.py`

Block average 面板：
- `errorbar()` pow2 点 + δSEM 误差棒
- x 轴 log2
- **高亮 plateau 窗口点**：plateau_index 及之后满足条件的点用绿色方形 marker，其余用蓝色圆形
- 平台处绿色水平线 + δSEM 色带
- 标题：`Block Average (plateau=yes/no)`

Summary text：`SEM_block` 替代旧 `SEM_block`/`SEM_arctan`，删除 `SEM_at_max_B`/`τ_implied`。

### Step 7: `__init__.py`

删除 `ArctanFitResult` 导出。

### Step 8: 测试

**删除**：`test_arctan_fit.py`、`test_block_average_dense.py`

**新增** `test_block_average.py`：
- pow2 生成正确性
- F&P plateau 在 AR(1) 数据上检测
- δSEM 公式验证
- n_b < 2 跳过
- 边界情况（常数序列、极短序列）

**修改** `test_workflow.py`：
- 3-tier → 2-tier 路径
- 新增 SEM_block vs SEM_auto 交叉验证 warning 测试

**修改** `test_constrained_ti_regression.py`：F&P plateau 断言替代 arctan 断言

### Step 9: context4agent 文档同步

## 验证

1. `pip install .`
2. `pytest test/unit/enhanced_sampling/constrained_ti/ -v`
3. `pytest test/integration/enhanced_sampling/ -v`
4. 对 5 个真实 TI 点跑 `standalone_diagnostics()` 出诊断图
