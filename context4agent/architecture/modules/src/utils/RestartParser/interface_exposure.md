# utils.RestartParser — 接口暴露

> 对应代码：`src/md_analysis/utils/RestartParser/__init__.py`
>
> 本文档定义 `RestartParser` 子包的符号级公开接口。

## 1. 接口角色定义

- `RestartParser` 负责解析 CP2K restart 文件和 LagrangeMultLog 文件。
- 公开范围由 `__init__.py` 的 `__all__` 严格定义（11 个符号）。
- 主要消费者：`enhanced_sampling.slowgrowth`、`enhanced_sampling.constrained_ti`、`scripts.TIGen`

## 2. 当前公开接口清单

### 2.1 `CellParser.py` 导出

异常：

- `CellParseError(MDAnalysisError)`
  - CP2K 文件格式错误时抛出（缺少 `&CELL` 块、缺少 A/B/C 向量、非正交 cell、缺少 ABC 行）

函数：

- `parse_abc_from_restart(restart_path: str | Path) -> tuple[float, float, float]`
  - 输入：CP2K `.restart` 文件路径
  - 输出：正交 cell 的 `(a, b, c)` 长度（Angstrom）
  - 语义：解析 `&CELL ... &END CELL` 块中 A/B/C 向量的对角元素

- `parse_abc_from_md_inp(md_inp_path: str | Path) -> tuple[float, float, float]`
  - 输入：CP2K 输入文件路径（如 `md.inp`）
  - 输出：正交 cell 的 `(a, b, c)` 长度（Angstrom）
  - 语义：匹配 `ABC [angstrom] a b c` 行

### 2.2 `ColvarParser.py` 导出

异常：

- `ColvarParseError(MDAnalysisError)`
  - 解析 COLVAR restart 或 LagrangeMultLog 文件失败时抛出

数据结构（frozen dataclass）：

- `ConstraintInfo`
  - COLLECTIVE 约束参数
  - 字段：
    - `colvar_id: int` — 约束序号
    - `target_au: float` — step_start 时刻的目标值（a.u.）
    - `target_growth_au: float` — 目标增长速率（per a.u. time，非 per step）
    - `intermolecular: bool` — 是否为分子间约束

- `ColvarInfo`
  - 约束集合容器
  - 字段：
    - `constraints: tuple[ConstraintInfo, ...]` — 有序约束元组
  - 方法：
    - `__len__() -> int` — 约束数量
    - `__getitem__(colvar_id: int) -> ConstraintInfo` — 按 colvar_id 查找（KeyError 如未找到）
    - `__iter__() -> Iterator[ConstraintInfo]` — 迭代约束
  - 属性：
    - `primary: ConstraintInfo` — 第一个约束

- `ColvarRestart`
  - restart 文件完整元数据
  - 字段：
    - `project_name: str` — 项目名称
    - `step_start: int` — 起始步数
    - `time_start_fs: float` — 起始时间（fs）
    - `timestep_fs: float` — MD 时间步（fs）
    - `total_steps: int` — 总步数
    - `colvars: ColvarInfo` — 约束信息
    - `lagrange_filename: str | None` — LagrangeMultLog 文件名（CP2K 不一定写出，可选）
    - `cell_abc_ang: tuple[float, float, float]` — 正交 cell 长度（Angstrom）
    - `fixed_atom_indices: tuple[int, ...] | None` — 固定原子索引（0-based，排序后）

- `LagrangeMultLog`
  - 拉格朗日乘子时序数据
  - 字段：
    - `shake: np.ndarray` — SHAKE 乘子（shape 取决于约束数/格式）
    - `rattle: np.ndarray` — RATTLE 乘子
    - `n_steps: int` — 步数
    - `n_constraints: int` — 约束数
  - 属性：
    - `collective_shake: np.ndarray` — shape `(n_steps,)`，集体变量 SHAKE 乘子
    - `collective_rattle: np.ndarray` — shape `(n_steps,)`，集体变量 RATTLE 乘子

- `ColvarMDInfo`
  - 完整慢增长 MD 会话：组合 restart 配置与 Lagrange 输出数据
  - 字段：
    - `restart: ColvarRestart` — 输入配置
    - `lagrange: LagrangeMultLog` — 输出乘子数据
  - 属性：
    - `n_steps: int` — 步数
    - `steps: np.ndarray` — 绝对步数 `[0, 1, ..., n_steps-1]`
    - `times_fs: np.ndarray` — 绝对时间（fs）
  - 方法：
    - `target_series_au(colvar_id: int | None = None) -> np.ndarray` — 返回正确对齐的 ξ(k) 序列（a.u.）
  - 工厂方法：
    - `ColvarMDInfo.from_paths(restart_path: str | Path, log_path: str | Path) -> ColvarMDInfo` — 一步解析两个文件

函数：

- `parse_colvar_restart(restart_path: str | Path) -> ColvarRestart`
  - 解析 CP2K COLVAR restart 文件
  - 使用 `finditer` 解析所有 `&COLLECTIVE` 块
  - 内部复用 `CellParser.parse_abc_from_restart()` 提取 cell 参数

- `parse_lagrange_mult_log(log_path: str | Path) -> LagrangeMultLog`
  - 解析 LagrangeMultLog 文件
  - 自动检测单约束/多约束格式
  - CP2K 溢出值 `***` 自动转换为 `nan`

- `compute_target_series(restart: ColvarRestart, n_steps: int, *, colvar_id: int | None = None) -> np.ndarray`
  - 重建 ξ(t) 序列（原子单位）
  - 公式：`xi(k) = target_au + (k - step_start) * target_growth_au * dt_au`
  - `dt_au = timestep_fs / AU_TIME_TO_FS`
  - `colvar_id` 可选参数：指定使用哪个约束，默认使用 primary

## 3. Stability

Stable — 被 `slowgrowth`、`constrained_ti` 和 `scripts.TIGen` 使用。
