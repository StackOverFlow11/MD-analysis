# `md_analysis.workflows` 接口暴露约定（当前实现）

> 对应代码：`src/md_analysis/workflows/`（`__init__.py` + `models.py` + 7 个域模块）
>
> 本文档定义编程入口 facade 包的公开接口与暴露边界。`run_*` 完整清单以 `workflows/__init__.py` 的 `__all__` 为权威。

## 1. 接口角色定义

- `workflows` 是 programmatic 入口的 **canonical** 模块：按域（water / potential / charge / calibration / enhanced_sampling / scripts / composite）拆分 `run_*` 函数，每个函数封装一种分析或工作目录生成，返回统一的 `WorkflowResult`。
- `md_analysis.main` 是同一批符号的 import-only 薄 re-export facade（`WorkflowResult` + 全部 30 个 `run_*`），仅为保留历史的 `from md_analysis.main import ...` 导入形态；新代码应使用 `md_analysis.workflows.<domain>` 显式路径。
- CLI（`md_analysis.cli`）的全部菜单命令（Phase 6.6 入口收敛后）在 `execute()` 内通过 `lazy_import("md_analysis.workflows.<domain>", "run_*")` 调度到本包，无底层业务函数直调。
- 本包只做 facade：组织参数、输出路径、artifact 收集与 metadata；科学计算保留在 `water/` / `electrochemical/` / `enhanced_sampling/` / `scripts/` 业务模块（各域的委托目标见第 3 节）。

## 2. `WorkflowResult` 契约（`models.py`）

```python
@dataclass(frozen=True)
class WorkflowResult:
    name: str                    # workflow 稳定标识（如 "water_three_panel"）
    output_dir: Path             # 所有 artifacts 所在目录
    artifacts: dict[str, Path]   # 短键 → 实际写盘文件/目录路径
    metadata: dict[str, Any]     # 轻量 scalars（帧数、方法、模式标志、单位）
    extra: Any | None = None     # 强类型 report 对象（可选）
```

- `to_dict()`：返回 JSON 友好字典；`extra` 若实现 `to_dict()` 则委托调用，否则原样透传。
- `require_artifacts_exist(result)`：校验 `artifacts` 中每个路径在磁盘上存在（相对路径按 `output_dir` 解析）；缺失时抛 `MissingArtifactError`（`MDAnalysisError` 子类），错误消息一次性列出全部 `(key, path)` 对。这是 workflow 测试与 CLI 的统一落盘断言点。
- 批量 workflow（`scripts.py`）的 artifacts 用 `workdir_<i>` 索引键（十进制枚举序号，非 zero-padding），artifacts 值也可以是目录（`require_artifacts_exist` 接受目录存在性）。

## 3. 公开接口清单（`__init__.py` 的 `__all__`，33 个符号）

### 3.1 模型（`models.py`）

- `WorkflowResult`、`MissingArtifactError`、`require_artifacts_exist`

### 3.2 水分析（`water.py`，5 个）

| 函数 | artifacts 键 | 委托目标 |
|---|---|---|
| `run_water_three_panel` | `density_csv` / `orientation_csv` / `adsorbed_profile_csv` / `adsorbed_range_txt` / `adsorbed_theta_csv` / `plot_png` | `water.plot_water_three_panel_analysis` |
| `run_water_density` | `density_csv` | `water.water_mass_density_z_distribution_analysis` |
| `run_water_orientation` | `orientation_csv` | `water.water_orientation_weighted_density_z_distribution_analysis` |
| `run_ad_water_orientation` | `adsorbed_profile_csv` / `adsorbed_range_txt` | `water.ad_water_orientation_analysis` |
| `run_ad_water_theta` | `theta_csv` | `water.compute_adsorbed_water_theta_distribution` |

单步 facade 直接写调用方给的 `output_dir`（不追加子目录）。

### 3.3 电势分析（`potential.py`，6 个）

| 函数 | artifacts 键 | 说明 |
|---|---|---|
| `run_potential_full` | `electrode_csv` 或 `center_csv`(+`fermi_csv`) / `phi_z_png` / `thickness_sensitivity_csv` | 复合入口，按输入在 `output_dir` 下追加 `electrode/`、`center/`、`fermi/`、`phi_z/`、`thickness_sensitivity/` 子目录；实际跑了哪些子分析记录在 `metadata["sub_analyses"]` |
| `run_center_potential` | `center_csv` | 单步 |
| `run_fermi_energy` | `fermi_csv` | 单步 |
| `run_electrode_potential` | `electrode_csv` | 单步（U vs SHE） |
| `run_phi_z_profile` | `phi_z_png` | 单步（唯一返回 PNG 的 potential 单步 facade） |
| `run_thickness_sensitivity` | `thickness_sensitivity_csv` | 单步（厚度扫描） |

全部 facade 同时支持 continuous（cube glob + 可选 `md_out_path` + xyz）与 distributed（`input_mode="distributed"` + `sp_*`）两种输入模式，由底层业务函数按 `input_mode` 分发。

### 3.4 电荷分析（`charge.py`，3 个）

| 函数 | 输出子目录 | artifacts 键 |
|---|---|---|
| `run_surface_charge` | `output_dir/<method>/`（双侧）或 `output_dir/<method>_<side>/`（单侧） | `charge_csv` / `charge_png` |
| `run_tracked_charge` | `output_dir/tracked/` | `tracked_charge_csv` / `tracked_charge_png` |
| `run_counterion_charge` | `output_dir/counterion_tracking/` | `counterion_charge_csv` / `counterion_summary_csv`（`counterion_charge_png` 仅在底层实际产出时列出） |

### 3.5 标定（`calibration.py`，2 个）

- `run_calibration_fit`：σ → φ 拟合并持久化 JSON。artifacts：`calibration_json`（必有）、`calibration_csv` / `calibration_png`（仅当传 `output_dir`）。`csv_path` 与 `data_points` 互斥。拟合质量指标（`r_squared` / `rmse` / `equation`）在 `metadata`，report 对象在 `extra`。
- `run_calibration_predict`：对已有 calibration JSON 的 **只读查询**，`artifacts` 恒为空；预测值在 `extra`（`CalibrationPredictionReport`）并镜像到 `metadata`。

### 3.6 增强抽样（`enhanced_sampling.py`，5 个）

| 函数 | artifacts 键 | `extra` |
|---|---|---|
| `run_slowgrowth_quick_plot` | `csv` / `quick_png` | `SlowgrowthAnalysisReport` |
| `run_slowgrowth_publication_plot` | `csv` / `publication_png` | 同上（两入口共享 `_run_slowgrowth_with_style`） |
| `run_ti_single_diagnostics` | `csv` / `diagnostics_png`（仅在底层产出时列出） | `ConstraintPointReport` |
| `run_ti_full_analysis` | `convergence_csv` / `free_energy_csv` / `free_energy_png` / `diagnostics_png_{i}`（逐点 PNG 按整数索引展平） | `TIFullAnalysisReport` |
| `run_ti_constant_potential_correction` | 上述 TI 基线全部键 + `corrected_free_energy_csv` / `corrected_free_energy_png` | `ConstantPotentialResult` |

`run_ti_full_analysis` / `run_ti_constant_potential_correction` 的 `strict: bool = True`（损坏约束点目录抛 `FileNotFoundError`；`False` 则 WARN+skip）。CLI 312/313 显式传 `strict=False` 保留菜单路径的宽松历史行为。

### 3.7 工作目录生成（`scripts.py`，8 个）

`run_bader_single` / `run_bader_batch` / `run_ti_single` / `run_ti_batch` / `run_potential_single` / `run_potential_batch` / `run_sp_single` / `run_sp_batch`：

- 只生成输入文件 / 工作目录，**不提交** PBS / SLURM / VASP / CP2K 作业。
- Single 入口 artifacts：`{"workdir": <dir>}`；Batch 入口 artifacts：`{"workdir_<i>": <dir_i>, ...}`。
- Batch `metadata` 固定包含 `n_successful` / `n_skipped` / `n_failed` 三计数 + `workdir_paths` 列表（当前底层 succeed-all-or-raise，后两计数恒为 0，键为契约占位）。
- `extra`：Bader / TI / SP 批量入口携带底层 `*BatchReport`；Potential 批量无 report 类型，`extra=None`。
- `run_ti_batch` 的 `colvar_id: int | None = None`、`overwrite: bool = False`、`verbose: bool = False`（`verbose` 仅控制 tqdm 进度条，不进 metadata；CLI 422 传 `overwrite=True` + `verbose=True`）。

### 3.8 复合入口（`composite.py`，1 个）

- `run_interface_analysis`：water 三联图 + 完整电势流水线，按惯例布局写入 `<output_dir>/water/` 与 `<output_dir>/electrochemical/potential/`。合并两个 leaf 的 artifacts；键冲突时抛 `RuntimeError`（不允许静默覆盖）。范围仅限 water + potential，取代已删除的 legacy `run_all`。

## 4. 推荐导入方式

- `from md_analysis.workflows import run_water_three_panel, WorkflowResult`
- `from md_analysis.workflows.water import run_water_density`（按域显式导入，新代码首选）
- `from md_analysis.workflows import require_artifacts_exist`（运行后落盘断言）
- `from md_analysis.main import run_potential_full`（兼容历史的顶层形态）

## 5. 非公开边界

- 各域模块内的私有 helper（`_run_slowgrowth_with_style` / `_common_metadata` / `_potential_distributed_kwargs` / `_build_batch_metadata` / `_flatten_workdirs_to_artifacts` / `_resolve_single_frame` 等）以下划线开头，不进 `__all__`，不属于公开契约。
- `__all__` 之外的任何名字（含模块内直接定义的临时对象）不承诺稳定性。

## 6. 稳定性级别与兼容承诺

- **Stable**：`run_*` 函数签名、`WorkflowResult` 字段形状、各函数的 artifacts 键（`test/unit/cli/test_workflow_wiring.py` 对 CLI 读取的 artifact 键做了契约锁定）。
- `metadata` 新增键视为非破坏性变更；删除或改义现有键属于破坏性变更。
- 入口重构期间删除的 legacy 名（`run_water_analysis` / `run_potential_analysis` / `run_charge_analysis` / `run_tracked_charge_analysis` / `run_counterion_charge_analysis` / `run_all`）未保留 alias；新名见第 3 节。

## 7. 接口变更触发条件

- `__all__` 新增 / 删除 / 重命名 `run_*` 或模型符号
- 任一 `run_*` 的签名、返回的 `WorkflowResult` 字段或 artifacts 键变化
- `WorkflowResult` / `MissingArtifactError` / `require_artifacts_exist` 行为变化
- 输出目录布局约定变化（如 charge 子目录、composite 惯例布局）

## 8. 接口变更流程（必须执行）

1. 更新 `src/md_analysis/workflows/` 对应模块与 `__init__.py.__all__`
2. 同步更新 `md_analysis/main.py` 的 re-export 清单
3. 更新本文档与同目录 `implementation_guidelines.md`
4. 运行 `test/unit/workflows/` 与 `test/unit/cli/test_workflow_wiring.py`（CLI 契约锁定）
5. 运行 `test/integration/test_main.py` 端到端验收
