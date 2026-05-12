# constrained_ti — 开发备忘

## 定位

约束热力学积分 (constrained TI) 的收敛诊断模块。对每个约束点的 Lagrange 乘子时间序列做四步诊断，判断采样是否充分。不从 `md_analysis.__init__` re-export — 需直接 `from md_analysis.enhanced_sampling.constrained_ti import ...`。

Agent 入口：`dispatch("ti_full_analysis", {root_dir, output_dir, parser, dir_filter, reverse, equilibration, epsilon_tol_ev, auto_equilibration, point_slice})`。调用 workflow 层 wrapper `run_ti_full_from_root(...)`（签名与 `TaskContract` 一致）。`agent/_handlers.py` 里的 handler 是**薄层**：只把 `TIFullAnalysisReport` 的文件路径路由到 `TaskResult.outputs`、JSON 可序列化数值路由到 `TaskResult.summary`，不含业务逻辑。

`summary.per_point`（14 个字段：`point_index` / `xi` / `n_analyzed` / `time_start_fs` / `time_end_fs` / `time_total_fs` / `tau_corr` / `n_eff` / `sem_final_au` / `sem_max_au` / `geweke_z` / `drift_D` / `passed` / `failure_reasons`）是**未来 Resources 层**判定"非平衡漂移 / N_eff 太少 / 遍历性假阳"等业务失败类别的信号源。

**Strict discovery 差异**：`run_ti_full_from_root` 调 `discover_ti_points(..., strict=True)`，任何被 `dir_filter` 选中但 metadata 解析失败的目录直接抛 `FileNotFoundError`（agent 映射为 `file_not_found`）；CLI 312 菜单路径保持 `strict=False`（跳过失败目录并 WARN）。

**IO 层架构**：discover/load 完全 engine-agnostic，靠 `md_analysis.engines.protocols.ConstraintMDParser` Protocol 适配；CP2K 实现是 `md_analysis.engines.cp2k.CP2KParser`（engines 包 import 时自动注册）。早期的 `enhanced_sampling/_parsers.py` shim 已在 engines 重构中删除，不要再创建本地 parser 缓存层。`TIPointDefinition` 字段为 `directory + parser + metadata(ConstraintMetadata)`，`xi` 是 property（从 `metadata.colvars.primary.target_au` 推导，不再从目录名解析）。`discover_ti_points(parser="auto", dir_filter=None)`：parser="auto" 嗅探注册 parser，dir_filter=None 用 `parser.is_constraint_directory` 做内容过滤（目录命名自由）；用户传 `dir_filter="ti_target_*"` 这种 glob 时，跳过 sniff 直接用第一个注册 parser（默认 CP2KParser）。新增 engine 适配 = 实现 Protocol + `register_parser` 注册，io.py 不动。

## 四步诊断流程

1. **ACF**（`analysis/autocorrelation.py`）→ τ_corr, N_eff, SEM_auto = σ√(2τ/N)
2. **Block averaging**（`analysis/block_average.py`）→ F&P (1989) pow2 平台检测 → SEM_block
3. **Running average**（`analysis/running_average.py`）→ drift D < factor × SEM
4. **Geweke**（`analysis/geweke.py`）→ stationarity z-test

## sem_final 选择

始终使用 block-average SEM：

```
F&P plateau reached → SEM at plateau start
otherwise           → SEM at largest valid block size
```

交叉验证：|SEM_block − SEM_auto| / max > 15% 时发出 warning（无论平台是否达到）。

## 模块结构

| 文件 | 用途 |
|---|---|
| `config.py` | 阈值常量（`DEFAULT_FP_MIN_BLOCKS=4`, `DEFAULT_FP_CONSECUTIVE=2`, `DEFAULT_CROSS_CHECK_RTOL=0.15` 等）+ 输出文件名 |
| `models.py` | frozen dataclass: `BlockAverageResult`, `AutocorrResult`, `RunningAverageResult`, `GewekeResult`, `ConstraintPointReport`, `TIReport` |
| `workflow.py` | 编排器：`analyze_single_point`, `analyze_standalone`, `analyze_ti`, `standalone_diagnostics`, `run_ti_full_from_root` (agent wrapper), `TIFullAnalysisReport`, `_parse_point_slice`, CSV 导出 |
| `plot.py` | 2×2 诊断图（running avg / ACF / block avg / summary）+ 自由能曲线图 |
| `integration.py` | 梯形积分权重、SEM targets、自由能积分 |
| `io.py` | Engine-agnostic 约束点目录发现 + 批量加载（`discover_ti_points(parser="auto", dir_filter=None, reverse, strict)`、`load_ti_series`）；解析委托给 `_parsers.py` |
| `correction.py` | 恒电势自由能修正（Nørskov）：`ConstantPotentialCorrection`, `ConstantPotentialResult`, `compute_constant_potential_correction`。`plot_corrected_free_energy_profile` 已迁至 `plot.py`，`correction.py` 仅 re-export 以保持旧导入路径 |
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

## Auto-equilibration（可选）

`analyze_standalone` / `analyze_ti` / `standalone_diagnostics` 均支持 `auto_equilibration=True`：
- 二分砍前半：每轮取后半段数据重跑四步诊断
- 通过 → 返回；数据不足（< `DEFAULT_AUTO_EQUIL_MIN_FRAMES=100`）→ 返回不收敛
- `failure_reasons` 中记录使用帧数和迭代次数
- 默认关闭（`False`），不影响现有行为
- `_is_converged()` 内部函数：TI 上下文用 `passed`，standalone 检查 geweke + running_avg + N_eff
- `_auto_equilibrate()` 内部函数：循环体调 `analyze_single_point`，完全复用现有诊断引擎

## 陷阱

- `analyze_block_average` 不再接受 `sem_auto`, `dense_sampling`, `arctan_*` 参数（2026-03-24 F&P 重构已删除 arctan）
- `_BLOCK_KEYS` 仅含 `{"min_blocks", "n_consecutive"}`，通过 `engine_overrides` 传递
- `ArctanFitResult` 和 `analysis/_arctan_fit.py` 已删除
- **符号约定**：`ConstraintPointReport.lambda_mean` 存储 CP2K 输出的原始 Shake 乘子 ⟨λ⟩；积分时由 `workflow.analyze_ti` 取反（`dA/dξ = −⟨λ⟩`）后存入 `TIReport.forces`，CSV 列 `dA_dxi` 和自由能图均使用取反后的值
- **Bug 9028656**：`correction.py` 的 `_get_electrode_area` 调用 `_sorted_frame_dirs(bader_dir)` 时缺少必需的 `dir_pattern` 参数，导致恒电位修正（菜单 313）运行时 TypeError。修复：补上 `DEFAULT_DIR_PATTERN`。跨包调用 `_frame_utils` 的私有函数时务必检查签名是否同步

## 恒电势修正（correction.py）

Nørskov 修正公式：`ΔF_Φ(ξ) = ΔF_q(ξ) + [σ(ξ) − σ_ref] × [Φ(ξ) − Φ_ref] × A / 2`

- σ 从各 `ti_target_*/bader/` 的 Bader 帧系综平均得到（`trajectory_surface_charge`）
- Φ 由 calibration mapper 从 σ 外推（`mapper.predict(σ)`）
- A 为电极表面积（从 POSCAR 晶胞计算，`AREA_VECTOR_INDICES`）
- 基准 = IS/FS 中点：`σ_ref = (σ_IS + σ_FS) / 2`，`Φ_ref = (Φ_IS + Φ_FS) / 2`（最小化最大修正量）
- 不做误差分析（修正项视为精确），保留 TIReport 的 λ 误差
- 缺少 bader/ 目录时 WARN 并跳过修正
- 依赖：`electrochemical.charge`（σ 计算）、`electrochemical.calibration`（σ→Φ）
