# utils — 开发备忘

## 定位

共享底层工具包，被所有分析模块和 `engines/` 依赖。**不**反向 import `engines/` 或业务模块。

Phase 2-7 重构后已划分为三个子层（`formats` / `structure` / `io`），物理常量集中在 `constants.py`。

## 三层结构

| 子层 | 职责 | 不应该做的事 |
|---|---|---|
| `utils/formats` | 解析单个文件 / 单类文件格式 | 不判断 workflow，不发现目录，不持有 engine registry |
| `utils/structure` | 几何 / 化学语义（金属层、水拓扑、周期聚类） | 不读 CP2K/VASP 文件 |
| `utils/io` | 路径发现、通用 CSV 写出、cell 文件选择调度 | 不解析具体文件内容 |

`engines/` 是 utils 之上的薄门面层，负责把 CP2K 的具体文件组合（`*.restart` + `*.LagrangeMultLog`、`md.out` + cube 等）映射成 engine-neutral 的 dataclass（`ConstraintMetadata` / `LambdaSeries` / `PotentialFrame` / `FermiRecord`）。

## 约定

- **`__init__.py` 不 re-export 任何符号**（`__all__ = []`）。所有调用方按子模块直接路径导入：
  - `from md_analysis.utils.constants import HA_TO_EV, DEFAULT_LAYER_TOL_A`
  - `from md_analysis.utils.formats.cube import read_cube_header_and_values, slab_average_potential_ev`
  - `from md_analysis.utils.formats.bader import load_bader_atoms`
  - `from md_analysis.utils.formats.cp2k_cell import parse_abc_from_md_inp`
  - `from md_analysis.utils.formats.cp2k_colvar import parse_colvar_restart, parse_lagrange_mult_log`
  - `from md_analysis.utils.formats.cp2k_stdout import parse_md_out_fermi, parse_sp_out_fermi`
  - `from md_analysis.utils.formats.cp2k_xyz import read_xyz_atoms_for_steps`
  - `from md_analysis.utils.structure.layer import detect_interface_layers`
  - `from md_analysis.utils.structure.water import detect_water_molecule_indices`
  - `from md_analysis.utils.structure.cluster import cluster_1d_periodic`
  - `from md_analysis.utils.io._frame_discovery import extract_step_time_from_dirname, discover_frame_dirs`
  - `from md_analysis.utils.io._io_helpers import _cumulative_average, _write_csv, _write_csv_from_arrays`
  - `from md_analysis.utils.io.cell_resolver import resolve_cell_abc`
- **新业务模块应优先调 `engines.cp2k` 门面（如 `read_constraint_metadata` / `read_fermi_series` / `read_continuous_potential_frames`）**，而不是直接走 `utils/formats/cp2k_*`；后者保留是给 engines 层用的实现细节。
- 下划线前缀函数（如 `_compute_bisector_cos_theta_vec`、`_write_csv_from_arrays`、`_parse_md_out_fermi`）被多模块 cross-layer 使用，路径仍为子模块直接路径，视为不稳定的内部依赖
- **两个 config.py**（重要！）：
  - `utils/constants.py`：物理常量（`AU_TIME_TO_FS`、`HA_TO_EV`、`BOHR_TO_ANG`）、cSHE 常量、默认参数、轴映射
  - `md_analysis/config.py`（上级目录）：用户持久化配置
- **单位约定**：距离 Å、能量 eV（内部 Hartree→eV 转换）、分数坐标 [0,1)、时间 fs

## 子目录

| 目录 | 用途 |
|---|---|
| `formats/` | 单文件解析器：cube / bader / cp2k_{cell,colvar,stdout,xyz} / vasp_{report,outcar,locpot}（VASP 三个为 Phase 8 占位） |
| `structure/` | 几何 / 拓扑：layer / water / cluster |
| `io/` | 路径与 IO 调度：`_frame_discovery` / `_io_helpers` / `cell_resolver` |

## 陷阱与历史 Bug

- **`formats/cube.py`**：CP2K 输出的 cube 文件使用 Fortran `D` 指数格式（如 `1.23D-04`），`_float()` 辅助函数将 `D` 替换为 `E`
- **`formats/cube.py`**：`read_cube_atoms(path, header)` 为公开函数，从 cube 文件头解析原子坐标 + cell → `ase.Atoms`
- **`formats/bader.py`**：POTCAR 中元素符号可能带 `_pv`/`_sv` 后缀，解析时需去除
- **`formats/cp2k_colvar.py`**：保留 `ColvarRestart = ConstraintMetadata` 和 `LagrangeMultLog = LambdaSeries` 两个 alias（Phase 5b rename 兼容层），canonical 名通过 `engines.models` 暴露
- **`formats/cp2k_stdout.py`** 与 **`formats/cp2k_xyz.py`**：Phase 7a 从 `electrochemical/potential/_frame_source.py` 抽出；后者现已退化为 thin wrapper（见 `electrochemical/potential/CLAUDE.md`）
- **`utils/formats/vasp_*.py`**（Phase 8 占位）：runtime **不**反向 import `engines/`，类型注解通过 `TYPE_CHECKING` + 字符串引用
- **常量精度**：`AU_TIME_TO_FS = 0.02418884326585`（CODATA 值），不要随意修改
