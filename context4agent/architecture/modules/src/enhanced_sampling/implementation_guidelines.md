# md_analysis.enhanced_sampling — Implementation Guidelines

## 职责边界

增强抽样方法的数据处理与可视化。当前包含两个子包：慢增长（slow-growth）和约束热力学积分（constrained TI）。

## 目录组织

```
enhanced_sampling/
  __init__.py          # 空包（re-export 预留）
  slowgrowth/          # 慢增长热力学积分子包
    __init__.py        # re-export 数据类 + 分析/绘图/CSV 函数
    config.py          # 输出文件名常量
    SlowGrowth.py      # 数据类 + 积分逻辑
    SlowGrowthPlot.py  # 绘图 + CSV 导出 + 统一入口
  constrained_ti/      # 约束 TI 收敛诊断 + 自由能积分
    __init__.py        # re-export public API
    config.py          # 阈值常量 + 输出文件名
    models.py          # frozen dataclasses + 异常层次
    analysis/          # 四个独立诊断引擎
      _acf_core.py     #   FFT 加速 ACF 计算
      autocorrelation.py # τ_corr, N_eff, SEM_auto
      block_average.py #   F&P pow2 平台检测
      running_average.py # 累积均值漂移检查
      geweke.py        #   平稳性 z 检验
    integration.py     # 梯形权重、误差传播、最优分配
    workflow.py        # 编排器 + standalone_diagnostics + CSV 导出
    plot.py            # 2×2 诊断图 + 自由能曲线图
    io.py              # 约束点目录发现 + 批量解析
```

## 依赖方向

### slowgrowth
- `slowgrowth` → `utils.RestartParser.ColvarParser`（`ColvarMDInfo`）
- `slowgrowth` → `utils.config`（`HA_TO_EV`）
- `slowgrowth` → `utils._io_helpers`（`_write_csv`）

### constrained_ti
- `constrained_ti.io` → `utils.RestartParser.ColvarParser`（仅 I/O）
- `constrained_ti.workflow` → `constrained_ti.analysis.*`、`integration`、`plot`、`io`
- `constrained_ti.workflow` → `utils.config`（`HA_TO_EV`）
- `constrained_ti.workflow` → `utils._io_helpers`（`_write_csv`）

### 模块内单向约束
- `workflow.py` 可导入 `plot.py`；`plot.py` **不得**导入 `workflow.py`
- `io.py` **不得**导入 `integration.py`（权重计算在 `workflow.py` 中）
- 各分析引擎之间**不得**互相导入
- **无反向依赖**：其他模块不依赖 `enhanced_sampling`

## constrained_ti 四步诊断流程

1. **Running Average**（`running_average.py`）— 累积均值漂移 D < drift_factor × SEM
2. **ACF**（`autocorrelation.py`）— Sokal (1997) 自洽截断 → τ_corr → N_eff ≥ 50
3. **Block Average**（`block_average.py`）— F&P (1989) pow2 块 → δSEM 平台检测
4. **Geweke**（`geweke.py`）— 前 10% vs 后 50% z 检验（|z| < 1.96）

### SEM 选择（2-tier）

```
F&P plateau reached → SEM_block (primary)
otherwise           → SEM_auto  (ACF fallback)
```

交叉验证：|SEM_block − SEM_auto| / max > 15% 时发出 warning。

### 符号约定

- `ConstraintPointReport.lambda_mean`：CP2K 输出的原始 SHAKE 乘子 ⟨λ⟩
- 积分时取反：`dA/dξ = −⟨λ⟩`（在 `workflow.analyze_ti` 中执行）
- λ 量纲为 Hartree/ξ_unit（a.u.），仅积分后的 ΔA 转换为 eV

### 积分方法

- 非均匀梯形积分：`compute_trapezoid_weights(xi)` 支持升序和降序 ξ
- `reverse=True` 时 ξ 降序排列，权重为负，ΔA = A(ξ_min) − A(ξ_max)
- SEM 目标分配：`SEM_max,k = ε_tol / (|w_k| × √K)`
- 误差传播：`σ_A = √(Σ (w_k × SEM_k)²)`

### 单位约定

- 用户 API 接受 eV（`epsilon_tol_ev`）
- 内部转换为 Hartree：`epsilon_tol_au = epsilon_tol_ev × EV_TO_HARTREE`
- SEM 目标与 λ 同单位（a.u.）
- 仅最终 ΔA 输出为 eV

## 扩展性

未来可在 `enhanced_sampling/` 下新增同级子包（如 `metadynamics/`），结构与 `slowgrowth/` 平行。顶层 `__init__.py` 可在需要时添加 re-export。
