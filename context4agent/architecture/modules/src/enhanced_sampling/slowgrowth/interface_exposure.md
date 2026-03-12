# md_analysis.enhanced_sampling.slowgrowth — Interface Exposure

## Public API

| Symbol | Module | Description |
|---|---|---|
| `Slowgrowth` | SlowGrowth.py | 基础数据类（frozen dataclass），包含分析就绪的数组 |
| `SlowgrowthFull` | SlowGrowth.py | 完整慢增长轨迹，从 restart + log 文件解析 |
| `SlowgrowthSegment` | SlowGrowth.py | `SlowgrowthFull` 的切片，自由能重新归零 |
| `slowgrowth_analysis` | SlowGrowthPlot.py | 统一入口：解析、切片、绘图 + CSV |
| `plot_slowgrowth_quick` | SlowGrowthPlot.py | 快速双轴图（CV 底轴 + MD step 顶轴） |
| `plot_slowgrowth_publication` | SlowGrowthPlot.py | 出版级双轴图（legend 标注能量值） |
| `write_slowgrowth_csv` | SlowGrowthPlot.py | 表格 CSV 导出 |

## Config Constants

| Constant | Value | Source |
|---|---|---|
| `DEFAULT_SG_QUICK_PNG_NAME` | `"slowgrowth_quick.png"` | config.py |
| `DEFAULT_SG_PUBLICATION_PNG_NAME` | `"slowgrowth_publication.png"` | config.py |
| `DEFAULT_SG_CSV_NAME` | `"slowgrowth_data.csv"` | config.py |

## 数据类层次

```
Slowgrowth (frozen dataclass)
├── steps: np.ndarray           # 绝对 MD 步编号
├── times_fs: np.ndarray        # 模拟时间 (fs)
├── target_au: np.ndarray       # 目标 CV 值 (a.u.)
├── lagrange_shake: np.ndarray  # Shake Lagrange 乘子
├── free_energy_au: np.ndarray  # 累积自由能变 (Hartree)
├── timestep_fs: float          # MD 时间步 (fs)
├── target_growth_au: float     # 每步 CV 增量 (a.u./step)，已从 ConstraintInfo 的 per-a.u.-time 值转换
├── n_steps -> int (property)
└── reversed() -> Slowgrowth    # 反转初末态

SlowgrowthFull(Slowgrowth)
├── md_info: ColvarMDInfo       # 原始解析元数据
├── from_md_info(cls, md_info, *, colvar_id) -> SlowgrowthFull
├── from_paths(cls, restart_path, log_path, *, colvar_id) -> SlowgrowthFull
└── segment(start, end, *, ref_step) -> SlowgrowthSegment

SlowgrowthSegment(Slowgrowth)
└── parent: SlowgrowthFull      # 源 Full 对象引用
```

## `slowgrowth_analysis` Signature

```python
def slowgrowth_analysis(
    restart_path: str,
    log_path: str,
    *,
    initial_step: int = 0,
    final_step: int | None = None,
    output_dir: Path | None = None,
    plot_style: str = "both",        # "quick" | "publication" | "both"
    colvar_id: int | None = None,
) -> dict[str, Path]
```

返回键：`"csv"`、`"quick_png"`（如适用）、`"publication_png"`（如适用）。

## `plot_slowgrowth_quick` Signature

```python
def plot_slowgrowth_quick(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    png_name: str = DEFAULT_SG_QUICK_PNG_NAME,
    absolute_steps: np.ndarray | None = None,
) -> Path
```

## `plot_slowgrowth_publication` Signature

```python
def plot_slowgrowth_publication(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    png_name: str = DEFAULT_SG_PUBLICATION_PNG_NAME,
) -> Path
```

## `write_slowgrowth_csv` Signature

```python
def write_slowgrowth_csv(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    csv_name: str = DEFAULT_SG_CSV_NAME,
) -> Path
```
