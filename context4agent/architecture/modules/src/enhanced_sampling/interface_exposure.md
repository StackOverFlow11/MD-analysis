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
                   equilibration=0, time_start_fs=0.0,
                   auto_equilibration=False,
                   **engine_overrides) → ConstraintPointReport

analyze_ti(xi_values, lambda_series_list, dt, *,
           epsilon_tol_ev=0.05, equilibration=0,
           time_starts=None, auto_equilibration=False,
           **engine_overrides) → TIReport

standalone_diagnostics(restart_path, log_path, *, equilibration=0,
                       sem_target=None, colvar_id=None,
                       output_dir=None,
                       auto_equilibration=False) → dict[str, Path | ConstraintPointReport]
```

### engines

引擎抽象层（半私有，下划线前缀）。

```python
class ConstraintMDParser(Protocol):
    name: str
    def is_constraint_directory(self, directory: Path) -> bool: ...
    def parse_metadata(self, directory: Path) -> ColvarRestart: ...
    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog: ...

CP2KParser  # 内置实现（recognises *.restart + *.LagrangeMultLog）
infer_parser(directory) -> ConstraintMDParser  # sniff 注册的 parser
register_parser(name, factory)                  # 添加新 engine（如 VASPParser）
get_parser(name) -> ConstraintMDParser
ParserInferenceError                             # sniff 失败异常
```

### constrained_ti.io

```python
discover_ti_points(root_dir, *, parser="auto", dir_filter=None,
                   reverse=False, strict=False) → list[TIPointDefinition]
    # parser: "auto" sniffs，或 ConstraintMDParser 实例 / 注册名（如 "cp2k"）
    # dir_filter: None=parser.is_constraint_directory；str=glob；callable=自定义
    # strict=True: 候选目录 metadata 解析失败 → FileNotFoundError
    # 排序：ξ 主键 + 目录名 tiebreaker；reverse=True → 降序

load_ti_series(point_defs) → list[tuple[float, np.ndarray, float]]
    # 返回 (xi, lambda_series, dt_fs)，未做均衡裁剪
    # metadata 已在 discover 阶段缓存，此处仅读 λ(t)
```

### slowgrowth.SlowGrowth

```python
SlowgrowthFull.from_directory(directory, *, parser="auto", colvar_id=None)
    # parser-driven 入口；与 from_paths(restart, log) 数值结果一致

SlowgrowthFull.from_paths(restart_path, log_path, *, colvar_id=None)
    # 旧入口，CP2K 专用；保留以向后兼容
```

## 顶层 re-export

`enhanced_sampling` **未**从 `md_analysis.__init__` re-export，需显式导入完整路径。
