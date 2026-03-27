# constrained_ti — 开发备忘

## 定位

约束热力学积分 (constrained TI) 的收敛诊断模块。对每个约束点的 Lagrange 乘子时间序列做四步诊断，判断采样是否充分。不从 `md_analysis.__init__` re-export — 需直接 `from md_analysis.enhanced_sampling.constrained_ti import ...`。

## 四步诊断流程

1. **ACF**（`analysis/autocorrelation.py`）→ τ_corr, N_eff, SEM_auto = σ√(2τ/N)
2. **Block averaging**（`analysis/block_average.py`）→ F&P (1989) pow2 平台检测 → SEM_block
3. **Running average**（`analysis/running_average.py`）→ drift D < factor × SEM
4. **Geweke**（`analysis/geweke.py`）→ stationarity z-test

## sem_final 选择（2-tier）

```
F&P plateau reached → SEM_block (primary)
otherwise           → SEM_auto  (ACF fallback)
```

交叉验证：|SEM_block − SEM_auto| / max > 15% 时发出 warning。

## 模块结构

| 文件 | 用途 |
|---|---|
| `config.py` | 阈值常量（`DEFAULT_FP_MIN_BLOCKS=4`, `DEFAULT_FP_CONSECUTIVE=2`, `DEFAULT_CROSS_CHECK_RTOL=0.15` 等） |
| `models.py` | frozen dataclass: `BlockAverageResult`, `AutocorrResult`, `RunningAverageResult`, `GewekeResult`, `ConstraintPointReport`, `TIReport` |
| `workflow.py` | 编排器：`analyze_single_point`, `analyze_standalone`, `analyze_ti`, `standalone_diagnostics`, CSV 导出 |
| `plot.py` | 2×2 诊断图（running avg / ACF / block avg / summary） |
| `integration.py` | 梯形积分权重、SEM targets、自由能积分 |
| `io.py` | 自动发现约束点目录 |
| `analysis/` | 四步诊断引擎 → `analysis/CLAUDE.md` |

## BlockAverageResult 字段

```python
block_sizes: np.ndarray      # pow2 block sizes
sem_curve: np.ndarray        # SEM(B) at each level
delta_sem: np.ndarray        # δSEM(B) = SEM(B) / √(2(n_b−1))
n_total: int                 # series length
plateau_index: int | None    # first plateau index
plateau_sem: float           # SEM at plateau (= SEM_block)
plateau_delta: float         # δSEM at plateau
plateau_block_size: int | None
plateau_reached: bool
passed: bool | None
```

## CSV 输出列

λ = dA/dξ 不是纯能量（量纲为 Hartree/ξ_unit），保持 a.u. 输出。
仅积分后的自由能 ΔA 转换为 eV。

```
xi, lambda_mean, sigma_lambda, tau_corr, n_eff, sem_auto,
sem_block, delta_sem_block, plateau_B, plateau_reached,
sem_final, sem_final_method, sem_max,
geweke_z, geweke_reliable, drift_D,
passed, failure_reasons
```

## 陷阱

- `analyze_block_average` 不再接受 `sem_auto`, `dense_sampling`, `arctan_*` 参数（2026-03-24 F&P 重构已删除 arctan）
- `_BLOCK_KEYS` 仅含 `{"min_blocks", "n_consecutive"}`，通过 `engine_overrides` 传递
- `ArctanFitResult` 和 `analysis/_arctan_fit.py` 已删除
- **符号约定**：`ConstraintPointReport.lambda_mean` 存储 CP2K 输出的原始 Shake 乘子 ⟨λ⟩；积分时由 `workflow.analyze_ti` 取反（`dA/dξ = −⟨λ⟩`）后存入 `TIReport.forces`，CSV 列 `dA_dxi` 和自由能图均使用取反后的值
