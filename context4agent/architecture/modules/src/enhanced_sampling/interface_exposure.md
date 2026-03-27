# md_analysis.enhanced_sampling — Interface Exposure

## 模块角色

增强抽样分析工作流的顶层包。包含 `slowgrowth`（慢增长）和 `constrained_ti`（约束 TI 收敛诊断）两个子包。

## Public API

`enhanced_sampling/__init__.py` 未 re-export 任何符号（空包），所有公开接口通过子包直接导入。

## 推荐导入方式

```python
# slowgrowth 子包
from md_analysis.enhanced_sampling.slowgrowth import (
    Slowgrowth,
    SlowgrowthFull,
    SlowgrowthSegment,
    slowgrowth_analysis,
    plot_slowgrowth_quick,
    plot_slowgrowth_publication,
    write_slowgrowth_csv,
)

# constrained_ti 子包
from md_analysis.enhanced_sampling.constrained_ti import (
    # 异常
    ConvergenceError,
    InsufficientSamplingError,
    # 数据模型
    ConstraintPointInput,
    ConstraintPointReport,
    TIReport,
    TIPointDefinition,
    RunningAverageResult,
    AutocorrResult,
    BlockAverageResult,
    GewekeResult,
    # 恒电位校正
    ConstantPotentialCorrection,
    ConstantPotentialResult,
)
from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    analyze_standalone,
    analyze_single_point,
    analyze_ti,
    standalone_diagnostics,
    write_convergence_csv,
    write_free_energy_csv,
    write_single_point_csv,
)
from md_analysis.enhanced_sampling.constrained_ti.io import (
    discover_ti_points,
    load_ti_series,
)
from md_analysis.enhanced_sampling.constrained_ti.plot import (
    plot_point_diagnostics,
    plot_free_energy_profile,
)
from md_analysis.enhanced_sampling.constrained_ti.correction import (
    compute_constant_potential_correction,
    write_corrected_free_energy_csv,
    plot_corrected_free_energy_profile,
)
```

## 关键函数签名

### constrained_ti.workflow

```python
analyze_standalone(series, *, dt=1.0, xi=0.0, sem_target=None,
                   equilibration=0, **engine_overrides) → ConstraintPointReport

analyze_ti(xi_values, lambda_series_list, dt, *,
           epsilon_tol_ev=0.05, equilibration=0,
           **engine_overrides) → TIReport

standalone_diagnostics(restart_path, log_path, *, equilibration=0,
                       sem_target=None, colvar_id=None,
                       output_dir=None) → dict[str, Path | ConstraintPointReport]
```

### constrained_ti.io

```python
discover_ti_points(root_dir, *, pattern="auto", reverse=False) → list[TIPointDefinition]
    # pattern: "ti_target" | "xi" | "auto"
    # reverse: True → ξ 降序（初态 = max ξ）

load_ti_series(point_defs) → list[tuple[float, np.ndarray, float]]
    # 返回 (xi, lambda_series, dt_fs)，未做均衡裁剪
```

## 顶层 re-export

`enhanced_sampling` **未**从 `md_analysis.__init__` re-export，需显式导入完整路径。
