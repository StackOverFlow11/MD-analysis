# utils — 开发备忘

## 定位

共享底层工具包，被所有分析模块和 `engines/` 依赖。**不**反向 import `engines/` 或业务模块。

Phase 1（破坏性重构）后 `formats/` 进一步按 engine 家族切分为 `cp2k/ vasp/ common/ bader/` 子包；`io/` 瘦身为纯通用 I/O helper，不再承载 CP2K `.restart` / `md.inp` 选择语义；物理常量集中在 `constants.py`。

## 三层结构

| 子层 | 职责 | 不应该做的事 |
|---|---|---|
| `utils/formats` | 解析单个文件 / 单类文件格式（按 engine 家族分子包） | 不判断 workflow，不发现目录，不持有 engine registry |
| `utils/structure` | 几何 / 化学语义（金属层、水拓扑、周期聚类） | 不读 CP2K/VASP 文件 |
| `utils/io` | 通用路径发现、通用 CSV 写出 | 不解析具体文件内容；不持有 engine 文件选择/调度语义 |

`engines/` 是 utils 之上的薄门面层，负责把 CP2K 的具体文件组合（`*.restart` + `*.LagrangeMultLog`、`md.out` + cube 等）映射成 engine-neutral 的 dataclass（`ConstraintMetadata` / `LambdaSeries` / `PotentialFrame` / `FermiRecord`）。

## 约定

- **`__init__.py` 不 re-export 任何符号**（`__all__ = []`）。所有调用方按子模块直接路径导入：
  - `from md_analysis.utils.constants import HA_TO_EV, DEFAULT_LAYER_TOL_A`
  - `from md_analysis.utils.formats.common.cube import read_cube_header_and_values, slab_average_potential_ev`
  - `from md_analysis.utils.formats.bader.acf import load_bader_atoms`
  - `from md_analysis.utils.formats.bader.potcar import _read_potcar_zval`
  - `from md_analysis.utils.formats.cp2k.cell import parse_abc_from_md_inp`
  - `from md_analysis.utils.formats.cp2k.colvar import parse_colvar_restart, parse_lagrange_mult_log`
  - `from md_analysis.utils.formats.cp2k.stdout import parse_md_out_fermi, parse_sp_out_fermi`
  - `from md_analysis.utils.formats.cp2k.xyz import read_xyz_atoms_for_steps`
  - `from md_analysis.utils.structure.layer import detect_interface_layers`
  - `from md_analysis.utils.structure.water import detect_water_molecule_indices`
  - `from md_analysis.utils.structure.cluster import cluster_1d_periodic`
  - `from md_analysis.utils.io._frame_discovery import extract_step_time_from_dirname, discover_frame_dirs`
  - `from md_analysis.utils.io._io_helpers import _cumulative_average, _write_csv, _write_csv_from_arrays`
- **新业务模块应优先调 `engines.cp2k` 门面（如 `read_constraint_metadata` / `read_fermi_series` / `read_continuous_potential_frames`）**，而不是直接走 `utils/formats/cp2k/*`；后者保留是给 engines 层用的实现细节。
- 下划线前缀函数（如 `_compute_bisector_cos_theta_vec`、`_write_csv_from_arrays`、`_parse_md_out_fermi`）被多模块 cross-layer 使用，路径仍为子模块直接路径，视为不稳定的内部依赖
- **两个 config.py**（重要！）：
  - `utils/constants.py`：物理常量（`AU_TIME_TO_FS`、`HA_TO_EV`、`BOHR_TO_ANG`）、cSHE 常量、默认参数、轴映射
  - `md_analysis/config.py`（上级目录）：用户持久化配置
- **单位约定**：距离 Å、能量 eV（内部 Hartree→eV 转换）、分数坐标 [0,1)、时间 fs

## 子目录

| 目录 | 用途 |
|---|---|
| `formats/common/` | engine-neutral 文件格式：`cube`（CP2K cube + 后续 LOCPOT-as-cube 共享） |
| `formats/cp2k/` | CP2K 单文件解析：`cell` / `colvar` / `stdout` / `xyz` |
| `formats/vasp/` | VASP 占位：`report` / `outcar` / `locpot`（hard-error，未注册到 parser registry） |
| `formats/bader/` | Bader 解析：`acf`（`_read_acf` + `load_bader_atoms`） / `potcar`（`_read_potcar_zval`） / `_errors`（`BaderParseError`） |
| `structure/` | 几何 / 拓扑：`layer` / `water` / `cluster` |
| `io/` | 路径与 IO 调度：`_frame_discovery` / `_io_helpers`（cell 文件选择已下沉到上层，不再放 `utils/io/`） |

## 陷阱与历史 Bug

- **`formats/common/cube.py`**：CP2K 输出的 cube 文件使用 Fortran `D` 指数格式（如 `1.23D-04`），`_float()` 辅助函数将 `D` 替换为 `E`
- **`formats/common/cube.py`**：`read_cube_atoms(path, header)` 为公开函数，从 cube 文件头解析原子坐标 + cell → `ase.Atoms`
- **`formats/bader/potcar.py`**：POTCAR 中元素符号可能带 `_pv`/`_sv` 后缀，解析时需去除
- **`formats/bader/_errors.py`**：`BaderParseError` 单独存放，避免 `acf.py` 与 `potcar.py` 循环 import（`acf.py` 同时依赖 `_errors` 和 `potcar`）
- **`formats/cp2k/colvar.py`**：解析器返回 CP2K 专有的 raw 类型（`Cp2kConstraintMetadataRaw` / `Cp2kLambdaSeriesRaw` 等）；engine-neutral canonical 类型（`ConstraintMetadata` / `LambdaSeries` / `ConstraintRun`）由 `engines.models` 持有，raw→canonical 转换在 `engines.cp2k`
- **`formats/cp2k/stdout.py`** 与 **`formats/cp2k/xyz.py`**：formats 抽取期间从 `electrochemical/potential/_frame_source.py` 抽出；后者现已退化为 thin wrapper（见 `electrochemical/potential/CLAUDE.md`）
- **`formats/vasp/*.py`** 占位：runtime **不**反向 import `engines/`，类型注解通过 `TYPE_CHECKING` + 字符串引用
- **常量精度**：`AU_TIME_TO_FS = 0.02418884326585`（CODATA 值），不要随意修改
