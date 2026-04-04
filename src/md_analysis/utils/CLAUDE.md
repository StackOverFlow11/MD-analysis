# utils — 开发备忘

## 定位

共享底层工具包，被所有分析模块依赖。包含物理常量、文件解析器、结构分析工具。

## 约定

- **`__init__.py` 是集中 re-export hub**：54 项在 `__all__` 中（均为公开符号）。下划线前缀函数（如 `_compute_bisector_cos_theta_vec`）在 `__init__.py` 中导入但不在 `__all__` 中，供 water 层内部 cross-layer 使用
- **两个 config.py**（重要！）：
  - `utils/config.py`：物理常量（`AU_TIME_TO_FS`、`HA_TO_EV`、`BOHR_TO_ANG`）、cSHE 常量、默认参数、轴映射
  - `md_analysis/config.py`（上级目录）：用户持久化配置
- **`_io_helpers.py`**：带下划线的私有模块，提供 `_cumulative_average()`、`_write_csv()`（dict rows）和 `_write_csv_from_arrays()`（numpy arrays），被全模块共享。所有 CSV 输出统一通过这两个函数
- **单位约定**：距离 Å、能量 eV（内部 Hartree→eV 转换）、分数坐标 [0,1)、时间 fs

## 陷阱与历史 Bug

- **CubeParser**：CP2K 输出的 cube 文件使用 Fortran `D` 指数格式（如 `1.23D-04`），`_float()` 辅助函数将 `D` 替换为 `E`
- **CubeParser**：`read_cube_atoms(path, header)` 为公开函数，从 cube 文件头解析原子坐标 + cell → `ase.Atoms`。被分布式电势分析（`_frame_source.py`）和 `PhiZProfile.py` 使用
- **BaderParser**：POTCAR 中元素符号可能带 `_pv`/`_sv` 后缀，解析时需去除
- **常量精度**：`AU_TIME_TO_FS = 0.02418884326585`（CODATA 值），不要随意修改

## 子目录

| 目录 | 用途 |
|---|---|
| `StructureParser/` | 金属层检测、水分子拓扑、周期聚类 → `StructureParser/CLAUDE.md` |
| `RestartParser/` | CP2K restart 和 LagrangeMultLog 解析 → `RestartParser/CLAUDE.md` |
