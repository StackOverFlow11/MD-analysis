# Phase 3 — engines 数据接口设计

> ## 实施状态(2026-05-14 更新;Phase 5 业务迁移收尾)
>
> 本文件为 **Phase 4 实施前的设计快照**;Phase 4 / Phase 5 实施过程中部分设计步骤
> (如 CP2K-名兼容 alias 的引入)是过渡态,不代表当前代码状态。
>
> | Phase / Commit | 状态 | 落地 commit |
> |---|---|---|
> | Phase 4 Commit 1(D10 + D8 constraint models migration) | ✅ 完成 | `63c4021` |
> | Phase 4 Commit 2(D7 CellSpec + read_cell) | ✅ 完成 | `e9f7e01`(+ `a84145f` doc-only) |
> | Phase 4 Commit 3(D11 CenterPotentialScalarFrame) | ✅ 完成 | `f4f4f77` |
> | Phase 5 Commit 1(water/cli cell 消费 → read_cell) | ✅ 完成 | `37107a5` |
> | Phase 5 Commit 2(SG/TI/cli 命名清理 + 删 CP2K-名兼容 alias) | ✅ 完成 | `3397a32` |
> | Phase 5 Commit 3(potential Fermi-only dict → typed) | ✅ 完成 | `6101873` |
>
> **Phase 5 业务迁移代码层收尾**;后续:
> - charge / correction 业务迁移留 Phase 8(VASP 阶段一并审视)
> - 入口层(cli / agent / scripts)迁移按 `overall_reconstruction_plan.md` Phase 6 推进
>
> **过渡态描述(已不再现存)**,继续读时请注意这些不是当前代码:
> - §3.3 ConstraintRun 中 `.restart` / `.lagrange` legacy alias property 段:
>   Phase 5 Commit 2 移除
> - §6.1 / §6.3 / §6.4 / §7 / §9 commit-1 步骤里提及的 `ColvarRestart` /
>   `LagrangeMultLog` / `ColvarMDInfo` 兼容 alias:Phase 5 Commit 2 移除
> - §11 Phase 5 表 + commit 边界设想已落地到下表
> - §11 `CenterPotential.parse_md_out_fermi` 选项 a/b 已拍板执行选项 b 轻量变体
>   (业务消费 typed `read_fermi_series`,不强连 `read_center_potential_scalar_frame`):
>   Phase 5 Commit 3 落地
>
> 之后做 Phase 6 入口层迁移计划时,以**代码现状 + 本 status banner** 为准,
> §3-§9 的设计契约文字保留为历史快照。
>
> 用途:作为 Phase 4 实施的**单一输入合约**;本文件 doc-only。
> 输入链路:
> - Phase 2 业务层数据需求 catalog(`phase2_business_data_requirements.md`)
> - Round 1/2/3 局部计划讨论(归档在 `temp/engines_design.md` /
>   `engine_design_round2.md` / `engine_design_round3.md`)
> - codex Round 4 通过 + 3 小约束
>
> 设计原则:Round 2 4 条红线
> 1. **不夹带科学行为变化** — 所有 facade / dataclass 输出值与现有 utils 层 1:1
> 2. **utils 不反向 import engines** — utils.formats.* 不出现 `import engines`
>    (TYPE_CHECKING 块除外)
> 3. **不过度抽象** — 没有业务消费的接口不进 Phase 4 实施;Protocol 不为单一
>    facade 而引入
> 4. **每个模型显式标注 Phase 4 实施 vs 仅设计暂不落地**

---

## 0. 设计目标 / 红线

### 0.1 设计目标

按 `overall_reconstruction_plan.md` Phase 3:

> 先设计,给用户审阅,确认后再实现。
> - 确定 `models`。
> - 确定 `protocols`。
> - 判断是否需要泛型 `ParserRegistry[T]`。
> - 判断哪些数据接口暂不抽象,只保留 facade。

本文件给出 6 个候选模型的完整设计 + 每个模型的 Phase 4 实施门槛。

### 0.2 红线(贯穿全文)

- **R1 — 不夹带科学行为变化**:任何新 dataclass / facade 的输出值必须与现有 utils
  层实现 numerically equal(byte-equal 优先);所有派生 property 公式与 utils 层 1:1
- **R2 — utils 不反向 import engines**:`utils.formats.cp2k.colvar` 等模块在 runtime
  不出现 `import engines`;TYPE_CHECKING 块允许字符串注解
- **R3 — engines 不算业务科学量**:engines 提供"原始/中间物理量"(Hartree potential、
  Fermi level、constraint metadata),业务公式(cSHE U_vs_SHE 转换、自由能积分、
  σ→φ 标定)**仍在业务层**
- **R4 — 实施 vs 设计分离**:本文件每个模型必须显式标注"Phase 4 实施"或"仅设计
  暂不落地",并给出门槛条件

### 0.3 Phase 4 实施 / 仅设计 一览(本文件结论摘要)

| 模型 | Phase 4 实施? | 原因 |
|---|---|---|
| `CellSpec`(§1) | **实施** | water + cli 主流程消费,2 处生产 caller |
| `StructureSnapshot`(§2) | **仅设计** | 当前业务无消费,作为未来扩展接口;Phase 4 不落地 |
| `ConstraintRun` composite(§3) | **实施** | SG + TI 主流程消费,5 处 caller |
| `ChargeFrame` / `ChargeTrajectory`(§4) | **仅设计** | 业务自有 `BaderTrajectoryData` 满足需求;Phase 4 不迁移 charge 业务 |
| `CenterPotentialScalarFrame`(§5) | **实施 dataclass + facade** | 引入 engines 入口,**但不**改 CenterPotential 业务代码(留 Phase 5) |
| `ConstraintMetadata`/`LambdaSeries` 物理迁移(§6) | **实施** | codex 红线 #2,Round 2 D10=A |
| `ParserRegistry[T]`(§8) | **仅设计** | 触发门槛(第二类 Protocol)未达成 |

---

## 1. `CellSpec`(D7 Layer 1) — Phase 4 实施

### 1.1 解决的问题

water `_common.py` + cli `_params.py` 直读 `utils.formats.cp2k.cell.parse_abc_from_md_inp`
拿 `(a, b, c)`;此 facade 不在 engines 层,违反 Phase 5 业务层"通过 engines 消费"
方向。同时 `ConstraintMetadata.cell_abc_ang` 已经是 tuple,但只能从 `.restart` 拿,
不从 `md.inp` 拿,且为正交特化形态。

### 1.2 字段契约

```python
@dataclass(frozen=True)
class CellSpec:
    """Engine-neutral cell descriptor."""

    cell_matrix_ang: np.ndarray
    pbc: tuple[bool, bool, bool] = (True, True, True)
```

| 字段 | shape / 类型 | 单位 | 约定 |
|---|---|---|---|
| `cell_matrix_ang` | `np.ndarray (3, 3)` | Å | **行向量约定**:每行 = 一条 lattice vector,与 `ase.Atoms.cell.array` / CP2K `&CELL A/B/C` 行序一致 |
| `pbc` | `tuple[bool, bool, bool]` | — | 默认 `(True, True, True)`,当前周期界面体系全 True |

### 1.3 派生属性(read-only `@property`)

```python
@property
def abc_ang(self) -> tuple[float, float, float]:
    """Norms of the three lattice vectors, in Å.

    NOTE: For non-orthorhombic cells this is NOT a box-edge length;
    callers that need "正交盒子长度" semantics MUST first check
    `is_orthorhombic` and fall back to a full-matrix branch otherwise.
    """

@property
def is_orthorhombic(self) -> bool:
    """True iff the off-diagonal elements of cell_matrix_ang are zero
    within a small numerical tolerance."""
```

**关键约束(codex Round 4 补充 #1)**:对非正交 cell,`abc_ang` 是**三条 lattice
vector 的范数**,不代表正交盒子长度;**只有 `is_orthorhombic=True` 时**才能被 water
当前逻辑当作 `(a, b, c)` 使用。docstring **必须**显式写明这一点。

### 1.4 Facade 签名

```python
# engines/cp2k.py
def read_cell(path: Path | str) -> CellSpec: ...
    """Read cell from a CP2K *.restart or md.inp file.

    Suffix-based dispatch (no content sniffing): any path whose
    ``Path.suffixes`` contains ``.restart`` (including bak variants
    such as ``*.restart.bak-1``) is sent to ``parse_abc_from_restart``;
    every other suffix (``md.inp`` / ``*.inp`` / etc.) is sent to
    ``parse_abc_from_md_inp``.
    Always returns CellSpec with cell_matrix_ang shape (3, 3); for
    md.inp's ABC-only format, the matrix is constructed as
    diag(a, b, c) and is_orthorhombic is True.
    """
```

### 1.5 不动什么(R3 边界)

- engines 不读 ase.Atoms;`read_cell` 只返回 cell,不返回 structure
- 业务层把 `CellSpec.cell_matrix_ang` 喂给 ase.Atoms 的责任仍在业务层(`_iter_trajectory`
  等)
- `parse_abc_from_md_inp` / `parse_abc_from_restart` 在 utils.formats.cp2k.cell
  保留(facade 内部委托给它们)

### 1.6 Phase 4 实施门槛

- ✅ 新增 `engines.models.CellSpec`
- ✅ 新增 `engines.cp2k.read_cell(path)` facade
- ✅ 在 `engines/__init__.py` re-export `CellSpec`
- ⚠️ Phase 4 **不**改业务 caller(留 Phase 5 业务迁移)

---

## 2. `StructureSnapshot`(D7 Layer 2) — 仅设计暂不落地

### 2.1 解决的问题

user 在 Round 2 提出"Cell + 原子序号 + 原子元素 + 分数坐标 + 原子受力"的复合载体。
当前业务**无强需求**(Bader 业务自有 ase.Atoms;water 业务用 ase.Atoms;potential
业务用 ase.Atoms);但 Phase 8 VASP placeholder + 未来力学 / 形变分析 可能消费此类
载体。

**Phase 4 不落地;本节只锁字段契约,作为未来扩展接口设计。**

### 2.2 字段契约(完整,作为未来锁定)

```python
@dataclass(frozen=True)
class StructureSnapshot:
    """Engine-neutral structural snapshot (single frame).

    Phase 3 status: DESIGNED but NOT implemented in Phase 4.
    Reserved for future force-aware analysis and VASP rollout.
    """

    cell: CellSpec
    atom_indices: np.ndarray
    atomic_numbers: np.ndarray
    symbols: tuple[str, ...]
    frac_coords: np.ndarray
    forces_ev_per_ang: np.ndarray | None = None
```

| 字段 | shape | 单位 | 约定 |
|---|---|---|---|
| `cell` | `CellSpec` | — | 见 §1 |
| `atom_indices` | `(n_atoms,) int` | — | **源文件原始顺序**(0-indexed);例如 CP2K xyz 中原子在文件里出现的位置。Bader 链路里 IndexMapper 用于在 POSCAR↔XYZ 之间映射时,POSCAR snapshot 的 `atom_indices` 反映回原始 XYZ 文件位置 |
| `atomic_numbers` | `(n_atoms,) int` | — | Z 数 |
| `symbols` | `tuple[str, ...]` 长 n_atoms | — | 化学符号(如 `("Cu", "Cu", "Ag", ...)`) |
| `frac_coords` | `(n_atoms, 3) float` | — | **wrap 到 `[0, 1)`**;消除整数周期表示歧义 |
| `forces_ev_per_ang` | `(n_atoms, 3) float \| None` | eV/Å | 默认 None;非 None 时 shape 必须 `(n_atoms, 3)` |

### 2.3 Facade 签名(草案,**Phase 4 不实施**)

```python
# engines/cp2k.py — placeholder docstring only, NOT implemented in Phase 4
def read_structure_snapshot(...) -> StructureSnapshot: ...
```

### 2.4 触发实施的条件

任一发生时,本节升格为 Phase 4+ 实施:
- 出现需要 force 的业务流程(elastic / phonon / NEB-like)
- VASP placeholder 实施需要 frame-level 结构载体
- Bader IndexMapper 决定走 engines-level 抽象(D9 触发)

---

## 3. `ConstraintRun` composite view(D8) — Phase 4 实施

### 3.1 解决的问题

`utils.formats.cp2k.colvar.ColvarMDInfo` 是 CP2K-specific composite,SG / TI / cli
共 5 处直接消费(`SlowGrowth.py:18,179` + `workflow.py:681` + `cli/_enhanced_sampling.py:44`
等)。需要在 engines 层提供同形 composite,业务层不再走 utils.formats.

### 3.2 字段契约

```python
@dataclass(frozen=True)
class ConstraintRun:
    """Engine-neutral composite view of a constraint-MD run.

    Combines metadata (restart-time inputs) with the resulting
    Lagrange-multiplier time series.  Derived series are computed
    via @property (not stored as fields) to keep this object an
    "inputs-only" snapshot.
    """

    metadata: ConstraintMetadata
    lambda_series: LambdaSeries
```

| 字段 | 类型 | 说明 |
|---|---|---|
| `metadata` | `ConstraintMetadata` | 见 §6;CP2K 端从 `.restart` 解析 |
| `lambda_series` | `LambdaSeries` | 见 §6;CP2K 端从 `*.LagrangeMultLog` 解析 |

### 3.3 派生属性(read-only `@property`)

完整保留 `ColvarMDInfo` 现有 4 个派生项,**加 2 个 Phase 5b legacy alias property**
(Phase 4 实施时引入,Phase 5 业务命名清理时删除):

```python
# Phase 5b legacy aliases — let business code (17 sites) using
# .restart / .lagrange keep working without changes; Phase 5
# naming cleanup removes these along with the 17 caller sites.
@property
def restart(self) -> ConstraintMetadata:
    """Phase 5b legacy alias for `metadata`."""
    return self.metadata

@property
def lagrange(self) -> LambdaSeries:
    """Phase 5b legacy alias for `lambda_series`."""
    return self.lambda_series

# Original 4 derived properties:
@property
def n_steps(self) -> int: ...

@property
def steps(self) -> np.ndarray: ...
    # shape (n_steps,) int, [0, 1, ..., n_steps-1]

@property
def times_fs(self) -> np.ndarray: ...
    # shape (n_steps,) float, steps * timestep_fs

def target_series_au(self, colvar_id: int | None = None) -> np.ndarray:
    # shape (n_steps,) float;CP2K 公式:
    #   xi(k) = target_au + (k - step_start) * target_growth_au * dt_au
    # where dt_au = timestep_fs / AU_TIME_TO_FS
```

**关于 `from_paths` classmethod**:历史 `ColvarMDInfo.from_paths(restart_path,
log_path)` classmethod **不**在 ConstraintRun 上重实现。理由:`from_paths` 需要调
parser + raw→canonical 转换 helper,转换 helper 在 `engines.cp2k`;若在
`engines.models.ConstraintRun` 上加 classmethod 调 `engines.cp2k.read_constraint_run_from_files`,
会形成 `engines.models → engines.cp2k` import cycle。**解决方案**:`from_paths`
作为 standalone facade `engines.cp2k.read_constraint_run_from_files` 提供(见
§3.4);3 处业务调用切到该 facade。

**关键约束(codex Round 3 §2.3)**:

- 派生 property 公式 **目前**与 `ColvarMDInfo` 完全一致(CP2K SHAKE/RATTLE 口径)
- **VASP adapter** 只有在确认 VASP REPORT 语义能映射到同一公式时才能构造
  `ConstraintRun`;否则:
  - **不允许** silent fallback(不能返回字段全 None / NaN 的"半 valid"对象)
  - **必须** raise `NotImplementedError` 或 engine-specific exception
- docstring 必须显式说明:"派生 property 公式按 CP2K SHAKE/RATTLE 约定;其他
  engine 实现 facade 时需在子类 docstring 里 override 公式说明"

### 3.4 Facade 签名

```python
# engines/cp2k.py — directory-level (the historical facade)
def read_constraint_run(directory: Path | str) -> ConstraintRun: ...
    """Read a CP2K constraint-MD point as a composite view.

    Internally:
      metadata = read_constraint_metadata(directory)
      lambda_series = read_lambda_series(directory)
      return ConstraintRun(metadata, lambda_series)

    Does NOT extend ConstraintMDParser Protocol; the engine-neutral
    layer simply composes two existing parser calls.
    """

# engines/cp2k.py — file-level facades (Phase 4 Commit 1 additions)
def read_constraint_metadata_from_restart(
    restart_path: Path | str,
) -> ConstraintMetadata: ...
    """Read a single CP2K *.restart file -> canonical ConstraintMetadata.

    File-level analogue of read_constraint_metadata(directory). Used by
    callers that already located the .restart file (e.g. scripts/TIGen
    and cli/_scripts SG preview).
    """

def read_lambda_series_from_log(
    log_path: Path | str,
) -> LambdaSeries: ...
    """Read a single CP2K *.LagrangeMultLog file -> canonical LambdaSeries.

    File-level analogue of read_lambda_series(directory).
    """

def read_constraint_run_from_files(
    restart_path: Path | str,
    log_path: Path | str,
) -> ConstraintRun: ...
    """Read a CP2K constraint-MD run from explicit (restart, log) paths.

    File-level analogue of read_constraint_run(directory); replaces the
    historical ColvarMDInfo.from_paths(...) class-method (which is NOT
    re-implemented on ConstraintRun to avoid an engines.models ->
    engines.cp2k import cycle; see §3.3).

    Internally:
      metadata = read_constraint_metadata_from_restart(restart_path)
      lambda_series = read_lambda_series_from_log(log_path)
      return ConstraintRun(metadata=metadata, lambda_series=lambda_series)
    """
```

### 3.5 不动什么

- **不**扩展 `ConstraintMDParser` Protocol(Round 3 §1 已锁:facade 内部组合已有
  parser)
- 现有 `read_constraint_metadata` / `read_lambda_series` 不动

`ColvarMDInfo` 在 Phase 4 Commit 1(constraint models migration,D10 + D8 原子
提交)物理迁移到 `engines.models`(作为 `ColvarMDInfo = ConstraintRun` 别名),
`utils.formats.cp2k.colvar` 删除定义、**不**保留任何 alias re-export(R2 硬
约束,见 §6.4 同 commit 切 caller + §9 Commit 1 改动范围)。

### 3.6 Phase 4 实施门槛

- ✅ 新增 `engines.models.ConstraintRun`(物理定义)
- ✅ 新增 `engines.cp2k.read_constraint_run(directory)` facade
- ✅ 在 `engines/__init__.py` re-export `ConstraintRun`
- ✅ `engines.models` 加 `ColvarMDInfo = ConstraintRun` 别名(物理定义在 engines)
- ✅ `utils.formats.cp2k.colvar` **删除** `ColvarMDInfo` 定义,**不**做 utils 层
  re-export(R2 硬约束)
- ✅ 同 commit 把 §6.4 表中所有 `ColvarMDInfo` import 站点(SG×2 + TI×1 + CLI×1
  lazy_import + 测试)切到 `engines.models`
- ⚠️ 与 §6(D10 dataclass 物理迁移)**强耦合**:`ConstraintRun` 字段引用 canonical
  `ConstraintMetadata` / `LambdaSeries`,这两个必须先完成物理迁移到 engines.models

---

## 4. `ChargeFrame` / `ChargeTrajectory`(D9) — 仅设计暂不落地

### 4.1 解决的问题

业务自有 `BaderTrajectoryData`(`electrochemical/charge/Bader/BaderData.py`)已满足
当前 σ 计算 + counterion 追踪流程。但 Phase 8 VASP / 未来 CP2K Mulliken /
Gaussian Hirshfeld 接入时,**需要统一的 charge 数据契约**。本节锁定该契约,**Phase 4
不实施**(charge 业务直读保留)。

### 4.2 `ChargeFrame` 字段契约

```python
@dataclass(frozen=True)
class ChargeFrame:
    """Single-frame engine-neutral charge data.

    Phase 3 status: DESIGNED but NOT implemented in Phase 4.
    Business code continues using BaderTrajectoryData.
    """

    step: int
    time_fs: float | None
    atom_indices: np.ndarray
    charges: np.ndarray
    atom_order: str
    is_remapped: bool
    source_engine: str
    charge_type: str
```

| 字段 | shape / 类型 | 约定 |
|---|---|---|
| `step` | `int` | MD step 号 |
| `time_fs` | `float \| None` | 仿真时间,fs;子目录名能提供时填,否则 None |
| `atom_indices` | `(n_atoms,) int` | 索引数组,与 `charges` 同长 |
| `charges` | `(n_atoms,) float` | 电荷数组;**对应顺序见 `atom_order`** |
| `atom_order` | `str` enum | `"poscar"` / `"xyz"` / `"snapshot"` / `"custom"`;明确 `charges[i]` 对应哪个原子序列 |
| `is_remapped` | `bool` | **该实例**是否已经做过 POSCAR→XYZ remap |
| `source_engine` | `str` enum | `"vasp"` / `"cp2k"` / `"unknown"` / `"custom"` (codex Round 4 补充 #2:**不**用 `"future"`;`"future"` 不是数据来源) |
| `charge_type` | `str` enum | `"bader_net"` / `"bader_electron_count"` / `"mulliken"` / `"hirshfeld"` |

### 4.3 单位 + 符号约定(docstring 强制要求)

- **`charges` 单位**:电子数 / e(无量纲;对 Bader 是 ZVAL - 原始 Bader 电子数;对
  Mulliken/Hirshfeld 是分配电子数)
- **Bader 符号约定**:`bader_net_charge = ZVAL - bader_charge`,**正号 = 失电子**
- **当前 Bader 数据源**:VASP 后处理链路(POSCAR + ACF.dat + POTCAR)

设计文档要求:future implementer 必须把以上 3 条**逐字**抄进 `ChargeFrame` 的
docstring,避免概念混淆。

### 4.4 `ChargeTrajectory` 字段契约

```python
@dataclass(frozen=True)
class ChargeTrajectory:
    """Trajectory-level charge data.

    Phase 3 status: DESIGNED but NOT implemented in Phase 4.
    """

    steps: np.ndarray
    times_fs: np.ndarray
    atom_indices: np.ndarray
    charges: np.ndarray
    atom_order: str
    is_remapped: bool
    source_engine: str
    charge_type: str
```

| 字段 | shape | 约定 |
|---|---|---|
| `steps` | `(n_frames,) int` | |
| `times_fs` | `(n_frames,) float` | |
| `atom_indices` | `(n_atoms,) int` | 假定跨帧不变;若帧间原子集合变化(counterion 检测场景),应用 `ChargeFrame` 列表而非此模型 |
| `charges` | `(n_frames, n_atoms) float` | 跨帧矩阵;**对应顺序见 `atom_order`** |
| 其他字段 | 同 `ChargeFrame` |

### 4.5 不动什么

- 业务层 `BaderTrajectoryData` 保留;Phase 4 不迁
- `electrochemical/charge/Bader/{BaderData,AtomCharges,SurfaceCharge}.py` 7 处
  `load_bader_atoms` 调用点保留
- `engines.charge` 子门面 Phase 4 **不**创建

### 4.6 触发实施的条件

任一发生时,本节升格为 Phase 4+ 实施:
- Phase 8 VASP placeholder 落地 + 业务要用统一 engine charge facade
- CP2K Mulliken 接入(`long_term.md` 已 flag)
- Bader 业务跨引擎复用(如 Gaussian Hirshfeld 在同一 σ 计算流程内)

---

## 5. `CenterPotentialScalarFrame`(D11) — Phase 4 实施(dataclass + facade only)

### 5.1 解决的问题

`electrochemical/potential/CenterPotential.py:447` 在 Fermi-only 分支直读
`parse_md_out_fermi`,绕过 engines;同时业务做 slab-averaged Hartree potential 计算
后,需要一个**标量级别**的帧载体(不带 cube values 巨数组,只带聚合标量),用于厚度
灵敏度 / 跨帧统计。

### 5.2 字段契约

```python
@dataclass(frozen=True)
class CenterPotentialScalarFrame:
    """Engine-neutral scalar-level frame for slab-averaged potential analysis.

    Distinct from PotentialFrame: this carries scalar aggregates
    (slab-centered Hartree potential, Fermi level, slab geometry)
    instead of the heavy (cube_path, values, header) raw payload.
    """

    step: int
    time_fs: float | None
    center_source: str
    center_z_ang: float | None
    slab_thickness_ang: float
    phi_center_ev: float
    fermi_level_ev: float | None
    phi_z_std_ev: float | None = None
    n_slices: int | None = None
```

| 字段 | 类型 | 单位 | 约定 |
|---|---|---|---|
| `step` | `int` | — | MD step 号 |
| `time_fs` | `float \| None` | fs | 同 `PotentialFrame.time_fs` |
| `center_source` | `str` enum | — | `"interface"`(界面检测)/ `"cell"`(几何中心)/ `"manual"`(用户给定) |
| `center_z_ang` | `float \| None` | Å | slab 中心 z 坐标(分数 → Å);若 `center_source="cell"` 可填 `cell_c / 2`,`None` 表示业务下游再决 |
| `slab_thickness_ang` | `float` | Å | slab 取的厚度(`thickness_sensitivity_analysis` 用) |
| `phi_center_ev` | `float` | **eV** | slab-averaged Hartree potential(`slab_average_potential_ev` 已 eV) |
| `fermi_level_ev` | `float \| None` | **eV** | Fermi level(业务从 Hartree 转 eV);连续模式拿 `md.out`,分布式模式拿 `sp.out` |
| `phi_z_std_ev` | `float \| None` | eV | slab 内 φ(z) 空间标准差(`thickness_sensitivity_analysis` 用) |
| `n_slices` | `int \| None` | — | slab 内 grid 切片数(诊断用) |

**单位后缀拍板(codex Round 3 §2.5)**:

- 中心势字段:**`phi_center_ev`**(eV) — 与 cSHE 公式中 `φ_center` 符号一致
- Fermi 字段:**`fermi_level_ev`**(eV) — 与 `phi_center_ev` 单位一致,业务下游
  公式 `U = -E_Fermi + φ_center + ...` 直接代入,不再二次单位转换
- **禁止**含糊命名:不用 `fermi_energy_level` / `fermi_value` / `potential`
- 业务层 `Hartree → eV` 转换的责任仍在 engines facade 内(facade 拿 `fermi_raw` Hartree
  乘 `HA_TO_EV`)

### 5.3 engines 不算什么(R3 + codex Round 2 D11 红线)

- ❌ **不**算 `U_vs_SHE` / `U_vs_RHE` / `U_vs_PZC`
- ❌ **不**使用 `DP_A_H3O_W_EV` / `MU_HPLUS_G0_EV` / `DELTA_E_ZP_EV` 等 cSHE
  reference 常数
- ❌ **不**做 reference scale 转换

以上仍在 `electrochemical.potential` 业务层。

### 5.4 Facade 签名(Phase 3 敲定)

```python
# engines/cp2k.py
def read_center_potential_scalar_frame(
    frame: PotentialFrame,
    *,
    center_z_ang: float,
    slab_thickness_ang: float,
    center_source: str = "manual",
) -> CenterPotentialScalarFrame: ...
```

**契约**:

| 参数 | 类型 | 说明 |
|---|---|---|
| `frame` | `PotentialFrame` | 已解析帧;facade 不重新读 cube / md.out。`frame.header` + `frame.values` 提供 cube 数据;`frame.fermi_raw` 提供 Hartree Fermi(facade 内部 `× HA_TO_EV` 得 `fermi_level_ev`);`frame.step` / `frame.time_fs` 直接透传 |
| `center_z_ang` | `float` | slab 中心 z 坐标(Å);**必须由 caller 提供**,facade **不**做 interface 检测 |
| `slab_thickness_ang` | `float` | slab 取的厚度(Å) |
| `center_source` | `str` | 仅作为 metadata 透传到 `CenterPotentialScalarFrame.center_source`;默认 `"manual"`,业务层若先做 interface 检测可传 `"interface"`(`"cell"` 也是合法值,语义见 §5.2) |

**facade 内部行为**(纯聚合,不引入新业务公式):

> **注意**:实际 `slab_average_potential_ev` 签名是
> `(header, values, thickness_ang, *, z_center_ang=None)` — 三个 positional
> 参数 + 一个 keyword-only 参数(`z_center_ang`,注意是 `z_center_ang` 不是
> `center_z_ang`)。下面 pseudo-code 严格按实际签名调用。

```
phi_center_ev, info = slab_average_potential_ev(
    frame.header,
    frame.values,
    slab_thickness_ang,
    z_center_ang=center_z_ang,
)
return CenterPotentialScalarFrame(
    step=frame.step,
    time_fs=frame.time_fs,
    center_source=center_source,
    center_z_ang=center_z_ang,
    slab_thickness_ang=slab_thickness_ang,
    phi_center_ev=phi_center_ev,
    fermi_level_ev=(
        frame.fermi_raw * HA_TO_EV if frame.fermi_raw is not None else None
    ),
    phi_z_std_ev=info.get("phi_z_std_ev"),
    n_slices=info.get("n_slices"),
)
```

**红线**(codex Round 5 MEDIUM):

- engines **不**做 interface 检测;业务层先调 `detect_interface_layers(atoms)` →
  `_extract_interface_geometry(...)` 得到 `center_z_ang`,再传 facade
- engines **不**读 `md.out` / `sp.out`;Fermi 走 `frame.fermi_raw`(已由
  `read_continuous_potential_frames` / `read_distributed_potential_frames` 在帧构造
  时解析)
- 这把 engines 的责任锁死在"已知 frame + 已知 center + 厚度 → 标量 reduce";业务
  层 interface 检测 + Fermi 选源都在 frame 构造前完成

### 5.5 Phase 4 实施门槛(**核心约束**)

- ✅ 新增 `engines.models.CenterPotentialScalarFrame`(dataclass)
- ✅ 新增 `engines.cp2k.read_center_potential_scalar_frame(...)` facade
- ✅ 在 `engines/__init__.py` re-export
- ⚠️ **Phase 4 不改 `electrochemical/potential/CenterPotential.py`**:
  - 业务保留 dict-style `parse_md_out_fermi` 访问(`engines/CLAUDE.md:73-78` 历史
    决议)
  - 业务保留 `slab_average_potential_ev` 直接调用
  - 新 facade **仅作为 Phase 5 业务迁移入口**;Phase 4 只交付接口
- ⚠️ 与 dict 接口迁移**解耦**:本节落地 dataclass + facade 不强制要求 dict 接口
  下线;Phase 5/6 各自独立立项

---

## 6. `ConstraintMetadata` / `LambdaSeries` 物理迁移(D10) — Phase 4 实施(**最高风险**)

### 6.1 解决的问题

当前 `ConstraintMetadata` / `LambdaSeries` **物理定义**在
`utils/formats/cp2k/colvar.py`,`engines.models` 从 utils 反 import 后 re-export。
现状**不**违反 import 边界(engines→utils 方向允许),但违反 user/codex 对齐时锁定的:

> `utils.formats.cp2k.*` 只做 CP2K 单文件解析 → 返回 CP2K-specific 原始结果或局部
> dataclass;`engines.cp2k` 负责把 CP2K-specific 结果组装/转换成 `engines.models`
> 中的 canonical 数据结构

### 6.2 目标终态

```
engines/models.py
    # ── Nested neutral types (physical definition here) ──
    @dataclass(frozen=True)
    class ConstraintInfo:        # neutral: collective-variable constraint
        colvar_id: int
        target_au: float
        target_growth_au: float
        intermolecular: bool

    @dataclass(frozen=True)
    class ColvarInfo:            # neutral: collection of constraints
        constraints: tuple[ConstraintInfo, ...]
        # property primary + __len__ / __getitem__ / __iter__ 保留

    # ── Canonical models (physical definition here) ──
    @dataclass(frozen=True)
    class ConstraintMetadata:
        project_name: str
        step_start: int
        time_start_fs: float
        timestep_fs: float
        total_steps: int
        colvars: ColvarInfo      # nested type is now neutral, lives in
                                 # engines.models — see §6.2.1
        lagrange_filename: str | None
        cell_abc_ang: tuple[float, float, float]
        fixed_atom_indices: tuple[int, ...] | None

    @dataclass(frozen=True)
    class LambdaSeries:
        shake: np.ndarray
        rattle: np.ndarray
        n_steps: int
        n_constraints: int
        # property collective_shake / collective_rattle 保留

utils/formats/cp2k/colvar.py
    # CP2K-specific raw types — field set identical to canonical ones
    # but flagged distinct so the conversion boundary is explicit.
    # Parsers return these; engines.cp2k converts them to canonical
    # engines.models types. NO runtime import of engines (R2).
    @dataclass(frozen=True)
    class Cp2kConstraintInfoRaw:
        colvar_id: int
        target_au: float
        target_growth_au: float
        intermolecular: bool

    @dataclass(frozen=True)
    class Cp2kColvarInfoRaw:
        constraints: tuple[Cp2kConstraintInfoRaw, ...]
        # No accessor methods on raw side; engines.cp2k builds the
        # accessor-bearing canonical ColvarInfo during conversion.

    @dataclass(frozen=True)
    class Cp2kConstraintMetadataRaw:
        project_name: str
        step_start: int
        time_start_fs: float
        timestep_fs: float
        total_steps: int
        colvars: Cp2kColvarInfoRaw
        lagrange_filename: str | None
        cell_abc_ang: tuple[float, float, float]
        fixed_atom_indices: tuple[int, ...] | None

    @dataclass(frozen=True)
    class Cp2kLambdaSeriesRaw:
        shake: np.ndarray
        rattle: np.ndarray
        n_steps: int
        n_constraints: int

    def parse_colvar_restart(path) -> Cp2kConstraintMetadataRaw: ...
    def parse_lagrange_mult_log(path) -> Cp2kLambdaSeriesRaw: ...

engines/cp2k.py
    def _cp2k_raw_to_constraint_info(raw) -> ConstraintInfo: ...
    def _cp2k_raw_to_colvar_info(raw) -> ColvarInfo: ...
    def _cp2k_raw_to_constraint_metadata(raw) -> ConstraintMetadata: ...
    def _cp2k_raw_to_lambda_series(raw) -> LambdaSeries: ...

    def read_constraint_metadata(directory) -> ConstraintMetadata:
        raw = parse_colvar_restart(...)
        return _cp2k_raw_to_constraint_metadata(raw)
    def read_lambda_series(directory) -> LambdaSeries:
        raw = parse_lagrange_mult_log(...)
        return _cp2k_raw_to_lambda_series(raw)
```

#### 6.2.1 嵌套类型也物理迁移(codex Round 5 HIGH 2)

`ConstraintMetadata.colvars: ColvarInfo` 的嵌套类型 `ColvarInfo` / `ConstraintInfo`
**字段本身已经是 engine-neutral 物理量**(`colvar_id` / `target_au` /
`target_growth_au` / `intermolecular`),但**物理定义位置**留在
`utils.formats.cp2k.colvar` 会让 canonical model 携带 CP2K-specific 语义,与 D10
"canonical model 完全 neutral"目标冲突(TYPE_CHECKING 解决 import 不解决语义泄漏)。

**修订决议**:`ColvarInfo` 与 `ConstraintInfo` **物理迁移到** `engines.models`,
成为 neutral nested contract;CP2K 端用同字段的 `Cp2kColvarInfoRaw` /
`Cp2kConstraintInfoRaw` raw 类型(平行结构、不同 import 路径),由 engines.cp2k
转换函数 `_cp2k_raw_to_colvar_info` / `_cp2k_raw_to_constraint_info` 一对一字段
映射成 canonical 实例。

VASP 端将来实现时:
- 用各自的 `VaspConstraintInfoRaw` / `VaspColvarInfoRaw`(字段集可能不同,VASP
  REPORT 表达若不同则 raw 端各自扩展)
- 用 `_vasp_raw_to_colvar_info(...)` 转换;**若 VASP 语义无法无损映射到 neutral
  `ConstraintInfo` 的 4 字段**,**禁止** silent fallback,**必须** raise
  `NotImplementedError`(与 §3.3 `ConstraintRun` VASP fallback 禁令一致)

**Raw 类型命名(codex Round 4 补充 #3)**:**不**用下划线开头。

- `Cp2kConstraintMetadataRaw` — 跨层转换契约,需要被 engines.cp2k 和测试 import,
  命名必须明确
- `Cp2kLambdaSeriesRaw` — 同上
- 下划线适合**私有 helper**(如 `_cp2k_raw_to_constraint_metadata`),不适合**跨
  层契约类型**

### 6.3 Phase 4 实施步骤(分阶段,降低回归风险)

#### Step 1 — 引入 raw 类型(纯新增,不破坏)

- 在 `utils/formats/cp2k/colvar.py` 新增 `Cp2kConstraintInfoRaw` /
  `Cp2kColvarInfoRaw` / `Cp2kConstraintMetadataRaw` / `Cp2kLambdaSeriesRaw` 四个
  dataclass(字段集与当前 `ConstraintInfo` / `ColvarInfo` / `ConstraintMetadata` /
  `LambdaSeries` 1:1)
- 在该文件**保留**现有 `ConstraintInfo` / `ColvarInfo` / `ConstraintMetadata` /
  `LambdaSeries` 定义(暂时)
- `utils.formats.cp2k.colvar` **不**从 engines import 任何东西(R2)

#### Step 2 — 物理迁移 canonical 定义

- 在 `engines/models.py` 新增 `ConstraintInfo` / `ColvarInfo` /
  `ConstraintMetadata` / `LambdaSeries` 的**物理定义**(含 `ColvarInfo` accessor
  methods + `LambdaSeries` `collective_shake` / `collective_rattle` property)
- 移除 `engines.models:41` 的 `from ..utils.formats.cp2k.colvar import ...` re-export
- 在 `utils/formats/cp2k/colvar.py` **删除**现有 `ConstraintInfo` / `ColvarInfo` /
  `ConstraintMetadata` / `LambdaSeries` 类定义
- `ColvarRestart` / `LagrangeMultLog` Phase 5b 历史 alias 也**物理搬到** engines.models
  (`ColvarRestart = ConstraintMetadata; LagrangeMultLog = LambdaSeries`);
  `utils.formats.cp2k.colvar` 不持有它们,不 import 它们(R2 保持)
- 旧 test 调用点 `from utils.formats.cp2k.colvar import ColvarRestart` 等
  **同 commit 切到** `from engines.models import ColvarRestart`(详见 §6.4)

#### Step 3 — parser 函数签名变更

- `parse_colvar_restart(path) -> ConstraintMetadata` → `-> Cp2kConstraintMetadataRaw`
- `parse_lagrange_mult_log(path) -> LambdaSeries` → `-> Cp2kLambdaSeriesRaw`

#### Step 4 — engines.cp2k 转换函数

- 新增 `_cp2k_raw_to_constraint_info(raw) -> ConstraintInfo`
- 新增 `_cp2k_raw_to_colvar_info(raw) -> ColvarInfo`
- 新增 `_cp2k_raw_to_constraint_metadata(raw) -> ConstraintMetadata`
- 新增 `_cp2k_raw_to_lambda_series(raw) -> LambdaSeries`
- 改造 `read_constraint_metadata` / `read_lambda_series` 的内部:从直接返回 parser
  输出 → 调用 raw → 转换 → 返回 canonical

#### Step 5 — `ColvarMDInfo` 物理迁移 + 同步切调用点(codex Round 5 HIGH 1)

> 破坏性重构原则,**不**保留 utils 层 alias。utils.formats.cp2k.colvar 在
> Phase 4 结束时**完全不持有** engines 类型,**不** runtime import engines
> (R2 全程满足)。

5a. **删除** `utils/formats/cp2k/colvar.py` 中的 `ColvarMDInfo` 定义。
5b. **不**在 utils 层做 `ColvarMDInfo = ConstraintRun` alias re-export(那会
    引入 utils → engines runtime import,违反 R2)。
5c. `ColvarMDInfo` 作为**别名**物理搬到 `engines.models`:`ColvarMDInfo =
    ConstraintRun`(让短期旧消费者 `from engines.models import ColvarMDInfo` 不破;
    本身**不**违反 R2,因为 alias 住在 engines)。
5d. **同 commit 把 §6.4 表中所有 `ColvarMDInfo` import 站点从
    `utils.formats.cp2k.colvar` 机械切到 `engines.models`**(详见 §6.4 新表):
    业务 5 处 + 测试 3 处 + CLI 2 处 + scripts 0 处(scripts 用 raw type,不切)。
    cli 的两处 lazy_import 字符串同步从 `"md_analysis.utils.formats.cp2k.colvar"`
    改成 `"md_analysis.engines.models"` /
    `"md_analysis.engines.cp2k"`(parse_colvar_restart 走 engines.cp2k facade)。
5e. agent 的异常 FQN(`ColvarParseError`)**不**切——它是 utils CP2K-specific
    exception,Phase 5b 之前的 alias **不**搬到 engines。它的 FQN 字符串保留
    `"md_analysis.utils.formats.cp2k.colvar.ColvarParseError"`。

**预估工作量**:Phase 1 同等规模(机械 import path rename;Phase 1 已成功完成
17 处)。Step 5d 涉及 10 处源码切 + 全套测试改 import,与 Phase 1 utils 重命名
同复杂度。

### 6.4 旧调用点迁移清单(**Phase 4 同 commit 机械切,不保留 utils 层 alias**)

破坏性重构原则:Phase 4 结束时 utils.formats.cp2k.colvar **不**持有任何 engines
canonical 类型 / engines alias,**不** runtime import engines。所有
`ColvarRestart` / `LagrangeMultLog` / `ColvarMDInfo` / `ColvarInfo` /
`ConstraintInfo` import 站点必须 **同 commit** 从 utils 路径切到 engines.models 路径。

`parse_colvar_restart` / `parse_lagrange_mult_log` 函数本身留在 utils(只签名换 raw),
直接调用 parser 的 site **可**保留 utils import,但**只**用 raw type;若它们之前
用 dataclass 字段访问(`.colvars`、`.shake` 等),需要换到 engines.cp2k facade
(`read_constraint_metadata` / `read_lambda_series`)拿 canonical 实例。

| 调用点 | 当前 import / 用法 | Phase 4 改造(同 commit 切) |
|---|---|---|
| `scripts/TIGen.py:19` | `from ..utils.formats.cp2k.colvar import parse_colvar_restart` + 用 `restart.colvars.primary.target_au` 等 | 切到 `from ..engines.cp2k import read_constraint_metadata_from_restart`;调用方式 `parse_colvar_restart(restart_path)` → `read_constraint_metadata_from_restart(restart_path)`;返 canonical,业务 `.colvars.primary` 等 accessor 不变(5 处 caller:line 19 + 284 / 486 / 567 / 674) |
| `test/unit/utils/test_slowgrowth_parser.py:11` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarInfo, ColvarMDInfo, ColvarParseError, ColvarRestart, ConstraintInfo, LagrangeMultLog, compute_target_series, parse_colvar_restart, parse_lagrange_mult_log)` | 拆分三组:① `ColvarInfo` / `ColvarMDInfo` / `ColvarRestart` / `ConstraintInfo` / `LagrangeMultLog` 切到 `md_analysis.engines.models`;② `compute_target_series` 切到 `md_analysis.engines.cp2k`(Phase 4 Commit 1 物理迁移);③ `ColvarParseError` / `parse_colvar_restart` / `parse_lagrange_mult_log` 保留 utils path |
| `test/unit/utils/test_slowgrowth_parser.py:367` | `from md_analysis.utils.formats.cp2k.colvar import _parse_fixed_atoms_list` | 不切(private utils helper) |
| `test/unit/scripts/test_ti_gen.py:22` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarRestart, ColvarInfo, ConstraintInfo)` | 全部切到 `from md_analysis.engines.models import (...)` |
| `test/unit/engines/test_facade.py:67` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarRestart, LagrangeMultLog)` | 切到 `from md_analysis.engines.models import (ColvarRestart, LagrangeMultLog)`;`test_facade.py:72-73` 的 `assert ColvarRestart is ConstraintMetadata` 仍成立(canonical alias 物理在 engines.models) |
| `engines/models.py:41`(re-export) | `from ..utils.formats.cp2k.colvar import ConstraintMetadata, LambdaSeries` | **删除该 re-export**;改为本地物理定义 |
| `engines/protocols.py:18` | `from .models import ConstraintMetadata, LambdaSeries` | 不动 |
| `engines/cp2k.py:33` | `from ..utils.formats.cp2k.colvar import (parse_colvar_restart, parse_lagrange_mult_log)` | 不动 import,但调用方式改:parser 返 raw → 走 `_cp2k_raw_to_constraint_metadata` / `_cp2k_raw_to_lambda_series` 转换 helper |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:18` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | 同 commit 切到 `from ...engines.models import ColvarMDInfo`(`ColvarMDInfo = ConstraintRun` alias 物理在 engines.models) |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:150` | `md_info = ColvarMDInfo.from_paths(restart_path, log_path)` | 切到 `md_info = read_constraint_run_from_files(restart_path, log_path)`;同 commit 加 `from ...engines.cp2k import read_constraint_run_from_files` |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:179` | 同 `:18`(`from_directory` 内 lazy import) | 同上 |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:186` | `md_info = ColvarMDInfo(restart=parser_obj.parse_metadata(directory), lagrange=parser_obj.parse_lambda_series(directory))` | kwarg rename:`metadata=parser_obj.parse_metadata(directory), lambda_series=parser_obj.parse_lambda_series(directory)`(ColvarMDInfo alias 不变,只换 canonical kwarg 名) |
| `enhanced_sampling/constrained_ti/workflow.py:681` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | 同 SlowGrowth:18 |
| `enhanced_sampling/constrained_ti/workflow.py:683` | `md_info = ColvarMDInfo.from_paths(restart_path, log_path)` | 同 SlowGrowth:150 |
| `cli/_enhanced_sampling.py:44` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo")` | 同 commit 切到 `lazy_import("md_analysis.engines.models", "ColvarMDInfo")` |
| `cli/_enhanced_sampling.py:47` | `info = ColvarMDInfo.from_paths(restart_path, log_path)` | 改 lazy_import 字符串 `lazy_import("md_analysis.engines.cp2k", "read_constraint_run_from_files")`,调用切到 `info = read_constraint_run_from_files(restart_path, log_path)` |
| `cli/_scripts.py:267` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart")` | 同 commit 切到 `lazy_import("md_analysis.engines.cp2k", "read_constraint_metadata_from_restart")`;局部变量同步 rename;业务 `restart.colvars.primary` 等 canonical accessor 不变 |
| `agent/_tasks_legacy.py:677` | FQN 字符串 `"md_analysis.utils.formats.cp2k.colvar.ColvarParseError"` | **不**切(`ColvarParseError` 是 utils CP2K-specific exception,保留在 utils) |

**Phase 4 测试基线门槛**:同 commit 切完上述站点后,全量 unit 729 / integration 44
**必须**不变;`rg "engines" src/md_analysis/utils --type py` 仅在 docstring 和
TYPE_CHECKING 块出现,**无** runtime import。

### 6.5 风险与缓解

| 风险 | 缓解 |
|---|---|
| parser 签名变化破坏测试 | Step 1-5 分阶段;每步跑分层测试;Step 3 之前 raw 与 canonical 并存 |
| `is ConstraintMetadata` 类型断言破裂 | `test/unit/engines/test_facade.py:72-73` 断言 `ColvarRestart is ConstraintMetadata`;Phase 4 把 `ColvarRestart` 别名物理放在 engines.models(`ColvarRestart = ConstraintMetadata`),test import 同 commit 切到 engines.models,断言仍成立 |
| **R2 红线**:utils 反向 import engines | **硬约束**(codex Round 5 HIGH 1):utils.formats.cp2k.colvar **不** runtime import engines;不持有 engines alias;raw 类型 + 转换 helper 双重隔离。任何引入 utils → engines 的过渡方案在 Phase 4 实施时都**禁止**采用 |
| §6.4 表中 11 处调用点迁移工作量 | Phase 1 同等机械 import path rename 规模(Phase 1 实际改了 17 处);commit 拆分后每笔可控;§6.4 表给出每处的精确切换目标 |

### 6.6 Phase 4 实施门槛

- ✅ Step 1-5 全部完成
- ✅ 全量 unit 729 / integration 44 基线不变
- ✅ R2 architecture scan 通过(`rg "engines" src/md_analysis/utils` 仅在 TYPE_CHECKING
  和 docstring)
- ✅ raw → canonical 转换函数有单测覆盖

---

## 7. 模型边界(Round 2 必备章节)

按 Round 2 codex 要求,本节单独列每个模型解决什么问题 + Phase 4 实施 vs 仅设计。

| 模型 | 解决的问题 | Phase 4 实施 | 仅设计原因(若适用) |
|---|---|---|---|
| `CellSpec` | water/cli 需要 engine-neutral cell 读取入口 | ✅ | — |
| `StructureSnapshot` | 未来 force-aware 业务 + VASP 阶段需要复合结构载体 | ❌ | 当前业务无消费;Phase 4 不引入未触发的抽象 |
| `ConstraintRun` | SG/TI/cli 5 处消费 `ColvarMDInfo` 需要 engine-neutral 入口 | ✅ | — |
| `ChargeFrame` / `ChargeTrajectory` | 未来跨引擎电荷数据统一契约 | ❌ | 业务自有 `BaderTrajectoryData` 满足;Phase 8 VASP 真启动后再实施 |
| `CenterPotentialScalarFrame` | engines 提供 slab-averaged Hartree + Fermi 标量帧 | ✅ (dataclass + facade) | 业务 caller 改造留 Phase 5 |
| `ConstraintMetadata` / `LambdaSeries` 物理迁移 | codex R2 + user 对齐红线 | ✅ | — |
| `ParserRegistry[T]` | 未来多类 Protocol 共享 registry | ❌ | 触发门槛(第二类 Protocol)未达成 |

---

## 8. Protocol / Registry(D12) — 仅设计

### 8.1 当前状态

```python
# engines/protocols.py
_REGISTRY: dict[str, Callable[[], ConstraintMDParser]] = {}

def register_parser(name, factory) -> None: ...
def get_parser(name) -> ConstraintMDParser: ...
def infer_parser(directory) -> ConstraintMDParser: ...
def resolve_parser(parser) -> ConstraintMDParser: ...
```

单一 Protocol(`ConstraintMDParser`)+ 单一 registry。

### 8.2 `ParserRegistry[T]` 设计草案

```python
# engines/protocols.py — Phase 3 design, NOT Phase 4 implementation
from typing import Generic, TypeVar

T = TypeVar("T")

class ParserRegistry(Generic[T]):
    """Generic engine-parser registry for one data contract."""

    def __init__(self, contract_name: str) -> None:
        self._contract_name = contract_name
        self._reg: dict[str, Callable[[], T]] = {}

    def register(self, name: str, factory: Callable[[], T]) -> None: ...
    def get(self, name: str) -> T: ...
    def infer(self, directory: Path) -> T: ...

# Concrete registries (future)
constraint_registry: ParserRegistry[ConstraintMDParser] = (
    ParserRegistry("constraint_md")
)
# When charge Protocol exists in Phase 8+:
# charge_registry: ParserRegistry[ChargeParser] = ParserRegistry("charge")
```

### 8.3 Phase 4 实施门槛(codex Round 3 §2.7)

**Phase 4 是否真的实现取决于:是否同时落地第二类稳定 Protocol。**

- 如果 Phase 4 只新增 `engines.cp2k.read_cell` / `read_constraint_run` /
  `read_center_potential_scalar_frame` 等**模块级 facade**(非 Protocol) → **不**
  重写现有 registry,继续用 `_REGISTRY` 单体
- 如果 Phase 4 同时引入新 Protocol(如 `ChargeParser`) → 才升级到 `ParserRegistry[T]`
- **不**为了 registry 重写本身破坏现有可工作的 `ConstraintMDParser` 注册路径

按 §4(charge 仅设计)+ §1(cell facade 是模块级)+ §3(composite 在 facade 内组合,
不扩展 Protocol)+ §5(scalar frame facade 是模块级),Phase 4 **不**引入第二类
Protocol。

**结论**:Phase 4 **不**实施 `ParserRegistry[T]`;保持当前 `_REGISTRY` 单体。

### 8.4 约束(codex Round 3 §2.7)

- VASP placeholder 仍**不得**注册成可用 parser
- registry 按"数据契约"分,不按 workflow 分(草案已遵循)
- 保留现有 `ConstraintMDParser` 兼容路径

---

## 9. Phase 4 实施门槛清单(commit 拆分)

Phase 4 总共 **3 笔 commit**(每笔后测试基线必须绿;**禁止**积压回归):

### Commit 1 — Constraint models migration(D10 + D8 原子提交,最高风险)

> codex Round 7 HIGH:D10 物理迁移 + D8 ConstraintRun + ColvarMDInfo 必须**单一
> 原子 commit**。理由:`ColvarMDInfo` 依赖 `ConstraintMetadata` / `LambdaSeries`;
> 若 D10 commit 后 ColvarMDInfo 仍留 utils 层,它要么 import engines(违反 R2),
> 要么改成依赖 raw 类型(但业务消费需要 canonical)。单一 commit 跳过这个矛盾态。

**改动范围**(同一 commit 完整完成,无中间不一致态):

1. **utils.formats.cp2k.colvar 新增 raw 类型**:`Cp2kConstraintInfoRaw` /
   `Cp2kColvarInfoRaw` / `Cp2kConstraintMetadataRaw` / `Cp2kLambdaSeriesRaw`
   (字段集与 canonical 1:1;R2:utils 不 import engines)
2. **engines.models 物理拥有 canonical 类型 + 别名**:
   - `ConstraintInfo` / `ColvarInfo`(嵌套 neutral,含 accessor methods,codex
     Round 5 HIGH 2)
   - `ConstraintMetadata` / `LambdaSeries`(canonical,含派生 property)
   - `ConstraintRun`(D8 composite)
   - 别名:`ColvarRestart = ConstraintMetadata` / `LagrangeMultLog = LambdaSeries` /
     `ColvarMDInfo = ConstraintRun`
3. **utils.formats.cp2k.colvar 删除原 canonical 定义**:`ConstraintInfo` /
   `ColvarInfo` / `ConstraintMetadata` / `LambdaSeries` / `ColvarMDInfo` /
   `ColvarRestart` / `LagrangeMultLog`(全部移除;**不**做 re-export)
4. **parser 函数签名变更**:
   - `parse_colvar_restart(path) -> Cp2kConstraintMetadataRaw`
   - `parse_lagrange_mult_log(path) -> Cp2kLambdaSeriesRaw`
5. **engines.cp2k 新增转换 helper + composite facade**:
   - `_cp2k_raw_to_constraint_info` / `_cp2k_raw_to_colvar_info` /
     `_cp2k_raw_to_constraint_metadata` / `_cp2k_raw_to_lambda_series`
   - 改造 `read_constraint_metadata` / `read_lambda_series` 内部:走 parser →
     raw → 转换 → canonical
   - 新增 `read_constraint_run(directory)` facade
6. **engines/__init__.py re-export**:`ConstraintRun` 加入 public API
7. **§6.4 表中所有 import 站点同 commit 切**(11 处机械 import path rename):
   - `ColvarRestart` / `LagrangeMultLog` / `ColvarInfo` / `ConstraintInfo` /
     `ColvarMDInfo` 从 `utils.formats.cp2k.colvar` 切到 `engines.models`
   - `cli/_enhanced_sampling.py:44` lazy_import 字符串切到 `engines.models`
   - `cli/_scripts.py:267` lazy_import 字符串按 §6.4 选定方向切
   - `scripts/TIGen.py:19` 按 §6.4 选定方向切
   - `agent/_tasks_legacy.py:677` 异常 FQN 字符串**不**切(保留 utils path)

**门槛**:
- utils.formats.cp2k.colvar **无** runtime engines import(R2 全程满足)
- utils 层**不**持有 engines canonical 或 engines alias
- 全量 unit + integration(729 + 44 基线)
- `test/unit/engines/test_facade.py:72-73` 类型断言专项验证
  (`ColvarRestart is ConstraintMetadata` 仍成立)
- `ConstraintRun.target_series_au` property 公式专项测试与 utils 旧
  `ColvarMDInfo.target_series_au` numerically equal
- import 边界:`rg "from .*engines|^import.*engines" src/md_analysis/utils` 命中
  全部必须落在 TYPE_CHECKING / docstring / comment 内(人工 / AST 复核;不依赖
  `grep -v TYPE_CHECKING` 行级过滤,见 §10)

### Commit 2 — D7 CellSpec

- 新增 `engines.models.CellSpec`
- 新增 `engines.cp2k.read_cell(path)` facade(内部委托 `parse_abc_from_md_inp` /
  `parse_abc_from_restart`)

**测试**:`CellSpec.abc_ang` / `is_orthorhombic` 派生 property 专项测试

### Commit 3 — D11 CenterPotentialScalarFrame

- 新增 `engines.models.CenterPotentialScalarFrame`
- 新增 `engines.cp2k.read_center_potential_scalar_frame(...)` facade
- **不**改 `CenterPotential.py` 业务代码

**测试**:facade 输出与现有 `slab_average_potential_ev` + Fermi 提取 numerically
equal

---

## 10. Phase 4 测试策略(codex Round 3 §2.8)

| 测试类型 | 范围 | 通过门槛 |
|---|---|---|
| **model property 小单测** | 每个新 dataclass 的派生 property:`CellSpec.abc_ang` / `CellSpec.is_orthorhombic` / `ConstraintRun.target_series_au` / `ConstraintRun.times_fs` 等 | 派生值正确,与 utils 旧实现 numerically equal(`np.testing.assert_allclose` atol=1e-12) |
| **CP2K raw → engines canonical 转换测试** | `_cp2k_raw_to_constraint_metadata` / `_cp2k_raw_to_lambda_series` | 转换前后字段值 byte-equal;字段集完整(用 `dataclasses.fields()` 对照) |
| **import 边界扫描** | `rg "from .*engines\|^import.*engines" src/md_analysis/utils --type py` 然后**人工 / AST 复核** | 通过门槛:**无 runtime engines import**(即所有命中必须在 `if TYPE_CHECKING:` 块内、docstring 内或注释内)。`grep -v TYPE_CHECKING` 行级过滤**不充分**——`utils/formats/vasp/{report,outcar}.py` 的 `if TYPE_CHECKING:` 块**下一行** import 不会被行级过滤排除;必须人工或 AST 确认每个命中位于 TYPE_CHECKING block / docstring / comment 之中(参考实现:解析 `ast.If` 节点 `test == Name("TYPE_CHECKING")`,只把 `body` 内 `Import` / `ImportFrom` 节点视为合规) |
| **data_example 回归** | `data_example/potential/` + `data_example/bader/` + `data_example/sg/` | integration 全过(44 passed) + CSV 列值 byte-equal(diff fixture 输出) |
| **业务 caller 不破** | full unit + integration | 729 + 44 基线不变 |

每笔 commit 单独跑测试,任一回归即修复;不允许"积压后统一处理"。

---

## 11. Phase 5 业务迁移影响面(基于 Phase 2 catalog 重算)

| 业务节点 | 当前状态(2026-05-14) | Phase 5 任务 |
|---|---|---|
| `water._common._parse_abc_from_md_inp` | ✅ 已切到 `engines.cp2k.read_cell(...).abc_ang`(Phase 5 Commit 1) | 完成 |
| `electrochemical.potential.CenterPotential.parse_md_out_fermi` (line 447) | ✅ 已切到 `engines.cp2k.read_fermi_series(...)` typed `FermiRecord`(Phase 5 Commit 3,选项 b 轻量变体)| 完成(`read_center_potential_scalar_frame` 业务消费推迟到 future Phase)|
| `electrochemical.charge.Bader.*` (7 处) | 保持 | 留 Phase 8 charge engines facade 落地 |
| `enhanced_sampling.slowgrowth.SlowGrowth.{from_paths,from_directory}` (2 处) | ✅ canonical `ConstraintRun`(Phase 5 Commit 2 命名清理完成,兼容 alias 已删) | 完成 |
| `enhanced_sampling.constrained_ti.workflow.standalone_diagnostics` (1 处) | ✅ 同上 | 完成 |
| `enhanced_sampling.constrained_ti.correction._get_electrode_area` (1 处) | 仍直读 `load_bader_atoms` | Phase 8 cleanup:同 charge facade 一起迁,或独立提"从 POSCAR 读 cell"小 helper |
| `cli/_params.py` `CellAbcParam.collect` | ✅ 已切到 `engines.cp2k.read_cell(...)` + source-vs-suffix gate(Phase 5 Commit 1)| 完成 |

Phase 5 commit 状态:
- ✅ **Commit 1**(`37107a5`):水/cli cell 迁移(D7 CellSpec 消费方激活)
- ✅ **Commit 2**(`3397a32`):SG/TI/cli 命名清理(`ColvarMDInfo` → `ConstraintRun`,
  删 CP2K-名兼容 alias 及 `.restart`/`.lagrange` legacy property)
- ✅ **Commit 3**(`6101873`):potential Fermi-only dict → typed
  (`fermi_energy_analysis` 切到 `read_fermi_series`;选项 b 轻量变体,
  `read_center_potential_scalar_frame` 业务消费推迟到 future Phase)
- ⬜ charge / correction 业务迁移留 Phase 8(VASP 阶段)

入口层(cli / agent / scripts)迁移**不在 Phase 5**;按 `overall_reconstruction_plan.md`
Phase 6 统一处理。

---

## 设计文档完结

已落地状态见**文件顶部 status banner**;后续 Phase 5 Commit 3 起步时,以
status banner + 代码现状为准。

注:本节末尾原"待办"流程已完成(本文件经过 codex Round 4-11 全部审阅 +
Phase 4 全 commit 实施 + Phase 5 Commit 1+2 实施)。
