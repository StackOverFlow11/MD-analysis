# Phase 3 — engines 数据接口设计

> 用途：作为 Phase 4(实施 engines)的**单一输入合约**。本文件 doc-only，不动代码。
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

    Auto-detects file type by extension and content sniffing.
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

完整保留 `ColvarMDInfo` 现有 4 个派生项:

```python
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
# engines/cp2k.py
def read_constraint_run(directory: Path | str) -> ConstraintRun: ...
    """Read a CP2K constraint-MD point as a composite view.

    Internally:
      metadata = read_constraint_metadata(directory)
      lambda_series = read_lambda_series(directory)
      return ConstraintRun(metadata, lambda_series)

    Does NOT extend ConstraintMDParser Protocol; the engine-neutral
    layer simply composes two existing parser calls.
    """
```

### 3.5 不动什么

- **不**扩展 `ConstraintMDParser` Protocol(Round 3 §1 已锁:facade 内部组合已有
  parser)
- `ColvarMDInfo` 在 `utils.formats.cp2k.colvar` 保留为**过渡期 alias**(`ColvarMDInfo
  = ConstraintRun`),让旧测试不破;Phase 5 业务迁移后再删
- 现有 `read_constraint_metadata` / `read_lambda_series` 不动

### 3.6 Phase 4 实施门槛

- ✅ 新增 `engines.models.ConstraintRun`(物理定义)
- ✅ 新增 `engines.cp2k.read_constraint_run(directory)` facade
- ✅ 在 `engines/__init__.py` re-export `ConstraintRun`
- ✅ `utils.formats.cp2k.colvar` 加 `ColvarMDInfo = ConstraintRun` alias(过渡期)
- ⚠️ Phase 4 **不**改业务 caller(留 Phase 5)
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

### 5.4 Facade 签名(草案)

```python
# engines/cp2k.py
def read_center_potential_scalar_frame(
    *,
    cube_path: Path,
    md_out_path: Path | None = None,
    sp_out_path: Path | None = None,
    center_source: str = "interface",
    slab_thickness_ang: float,
    ...
) -> CenterPotentialScalarFrame: ...
```

具体签名细节(参数取自现有 `slab_average_potential_ev` + `parse_md_out_fermi` /
`parse_sp_out_fermi` 的输入集合)由 Phase 4 实施时按调用点最小化原则拍板。

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
    @dataclass(frozen=True)
    class ConstraintMetadata:    # physical definition here
        project_name: str
        step_start: int
        time_start_fs: float
        timestep_fs: float
        total_steps: int
        colvars: ColvarInfo      # NOTE: ColvarInfo / ConstraintInfo
                                 # are CP2K-specific sub-types — see §6.4
        lagrange_filename: str | None
        cell_abc_ang: tuple[float, float, float]
        fixed_atom_indices: tuple[int, ...] | None

    @dataclass(frozen=True)
    class LambdaSeries:           # physical definition here
        shake: np.ndarray
        rattle: np.ndarray
        n_steps: int
        n_constraints: int
        # property collective_shake / collective_rattle 保留

utils/formats/cp2k/colvar.py
    @dataclass(frozen=True)
    class Cp2kConstraintMetadataRaw:   # CP2K-specific raw type
        # fields equivalent to current ConstraintMetadata
        ...
    @dataclass(frozen=True)
    class Cp2kLambdaSeriesRaw:         # CP2K-specific raw type
        ...
    def parse_colvar_restart(path) -> Cp2kConstraintMetadataRaw: ...
    def parse_lagrange_mult_log(path) -> Cp2kLambdaSeriesRaw: ...
    # NO runtime import of engines (R2)

engines/cp2k.py
    def _cp2k_raw_to_constraint_metadata(raw) -> ConstraintMetadata: ...
    def _cp2k_raw_to_lambda_series(raw) -> LambdaSeries: ...
    def read_constraint_metadata(directory) -> ConstraintMetadata:
        raw = parse_colvar_restart(...)
        return _cp2k_raw_to_constraint_metadata(raw)
    def read_lambda_series(directory) -> LambdaSeries:
        raw = parse_lagrange_mult_log(...)
        return _cp2k_raw_to_lambda_series(raw)
```

**Raw 类型命名(codex Round 4 补充 #3)**:**不**用下划线开头。

- `Cp2kConstraintMetadataRaw` — 跨层转换契约,需要被 engines.cp2k 和测试 import,
  命名必须明确
- `Cp2kLambdaSeriesRaw` — 同上
- 下划线适合**私有 helper**(如 `_cp2k_raw_to_constraint_metadata`),不适合**跨
  层契约类型**

### 6.3 Phase 4 实施步骤(分阶段,降低回归风险)

#### Step 1 — 引入 raw 类型(纯新增,不破坏)

- 在 `utils/formats/cp2k/colvar.py` 新增 `Cp2kConstraintMetadataRaw` /
  `Cp2kLambdaSeriesRaw` 两个 dataclass(字段集与当前 `ConstraintMetadata` /
  `LambdaSeries` 1:1)
- 在该文件**保留**现有 `ConstraintMetadata` / `LambdaSeries` 定义(暂时)

#### Step 2 — 物理迁移 canonical 定义

- 在 `engines/models.py` 新增 `ConstraintMetadata` / `LambdaSeries` 的**物理定义**
- 移除 `engines.models:41` 的 `from ..utils.formats.cp2k.colvar import ...` re-export
- 在 `utils/formats/cp2k/colvar.py` 删除现有 `ConstraintMetadata` / `LambdaSeries`
  类定义
- `ColvarRestart = ConstraintMetadata` / `LagrangeMultLog = LambdaSeries` 两个
  Phase 5b 历史 alias **更新**:从 `engines.models` re-export(让旧测试不破)
  — 这意味着 `utils.formats.cp2k.colvar` **需要**从 engines re-export 来满足
  alias。**这违反 R2**(utils 不反向 import engines)。
- **解决方案**:这两个 alias **物理搬到** `engines.models`(`ColvarRestart =
  ConstraintMetadata; LagrangeMultLog = LambdaSeries`),`utils.formats.cp2k.colvar`
  不再持有它们。旧测试如果 import `from utils.formats.cp2k.colvar import ColvarRestart`
  会破——见 §6.4 旧调用点清单

#### Step 3 — parser 函数签名变更

- `parse_colvar_restart(path) -> ConstraintMetadata` → `-> Cp2kConstraintMetadataRaw`
- `parse_lagrange_mult_log(path) -> LambdaSeries` → `-> Cp2kLambdaSeriesRaw`

#### Step 4 — engines.cp2k 转换函数

- 新增 `_cp2k_raw_to_constraint_metadata(raw) -> ConstraintMetadata`(私有 helper,
  下划线 OK)
- 新增 `_cp2k_raw_to_lambda_series(raw) -> LambdaSeries`
- 改造 `read_constraint_metadata` / `read_lambda_series` 的内部:从直接返回 parser
  输出 → 调用 raw → 转换 → 返回 canonical

#### Step 5 — `ColvarMDInfo` 处理

- `ColvarMDInfo` 当前在 `utils/formats/cp2k/colvar.py`,且消费 `ConstraintMetadata`
  / `LambdaSeries`
- 选项 A:`ColvarMDInfo` 也搬到 `engines.models`,作为 `ConstraintRun` 的 alias
  (`ColvarMDInfo = ConstraintRun`)
- 选项 B:删 `ColvarMDInfo`,业务全切到 `ConstraintRun`
- **Phase 4 选 A**:保留 alias 让 SG/TI/cli 调用点暂时不破,Phase 5 业务迁移时切

### 6.4 旧调用点迁移清单(**Phase 4 实施时必须保证测试通过**)

| 调用点 | 类型 | Phase 4 改造方式 |
|---|---|---|
| `scripts/TIGen.py:19` | `from ..utils.formats.cp2k.colvar import parse_colvar_restart` + 用 `restart.colvars` | parser 返 `Cp2kConstraintMetadataRaw` 后,`raw.colvars` 字段名需要保持 — 见 Note 1 |
| `test/unit/utils/test_slowgrowth_parser.py` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarInfo, ColvarMDInfo, ColvarParseError, ColvarRestart, ConstraintInfo, LagrangeMultLog, compute_target_series, parse_colvar_restart, parse_lagrange_mult_log)` | `ColvarRestart` / `LagrangeMultLog` 改为 `from md_analysis.engines.models import (...)` 别名;`ColvarMDInfo` 改为 `from md_analysis.engines.models import ConstraintRun as ColvarMDInfo` |
| `test/unit/scripts/test_ti_gen.py:22` | 同上 import 模式 | 同上 |
| `test/unit/engines/test_facade.py:67` | `from md_analysis.utils.formats.cp2k.colvar import ColvarRestart, LagrangeMultLog` | 同上 |
| `engines/models.py:41`(re-export) | `from ..utils.formats.cp2k.colvar import ConstraintMetadata, LambdaSeries` | 删除该 import;改为本地定义 |
| `engines/protocols.py:18` | `from .models import ConstraintMetadata, LambdaSeries` | 不动(已经经过 engines.models) |
| `engines/cp2k.py:33` | `from ..utils.formats.cp2k.colvar import (parse_colvar_restart, parse_lagrange_mult_log)` | 不动 import,但调用方式变(从直接返回到走 raw 转换 helper) |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:18,179` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | Phase 4 保留 alias;Phase 5 切到 `from ...engines.models import ConstraintRun` |
| `enhanced_sampling/constrained_ti/workflow.py:681` | 同上 | 同上 |
| `cli/_enhanced_sampling.py:44` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo")` | Phase 4 保留;Phase 6 入口层迁移再切 |
| `cli/_scripts.py:267` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart")` | Phase 4 保留;Phase 6 切 |
| `agent/_tasks_legacy.py:677` | FQN 字符串 `"md_analysis.utils.formats.cp2k.colvar.ColvarParseError"` | `ColvarParseError` **不**迁(它是 utils CP2K-specific exception),保留 |

**Note 1**:`ColvarInfo` / `ConstraintInfo` 是 `ConstraintMetadata.colvars` 字段的
子类型,**Phase 4 决议保留在 `utils.formats.cp2k.colvar` 作为 CP2K-specific 嵌套
类型**——它们不是 engine-neutral 的(VASP 端的 collective variable 表达可能不同)。
`ConstraintMetadata.colvars: ColvarInfo` 字段类型注解在 `engines.models` 里走
TYPE_CHECKING + 字符串引用,运行时不导致 utils → engines 反向 import。

### 6.5 风险与缓解

| 风险 | 缓解 |
|---|---|
| parser 签名变化破坏测试 | Step 1-5 分阶段;每步跑分层测试;Step 3 之前 raw 与 canonical 并存 |
| `is ConstraintMetadata` 类型断言破裂 | `test/unit/engines/test_facade.py:72-73` 断言 `ColvarRestart is ConstraintMetadata`,Phase 4 后 `ColvarRestart` 物理搬到 engines.models 仍然成立 |
| utils 反向 import engines(R2 违反) | parser 函数 **不** import engines;只用 raw 类型;转换在 engines.cp2k 完成 |
| 业务/CLI/agent 调用点回归 | `ColvarMDInfo` / `ColvarRestart` / `LagrangeMultLog` 保留为 engines.models 别名;Phase 4 业务零改动 |

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

Phase 4 总共 **4 笔 commit**(每笔不超过当前已知最大 commit 规模):

### Commit 1 — D10 物理迁移(最高风险,先做)

- 引入 `Cp2kConstraintMetadataRaw` / `Cp2kLambdaSeriesRaw`
- 物理迁移 `ConstraintMetadata` / `LambdaSeries` 到 `engines.models`
- 物理迁移 `ColvarRestart` / `LagrangeMultLog` alias 到 `engines.models`
- parser 函数签名变更 + raw → canonical 转换 helper
- `ColvarMDInfo` 暂留 `utils/formats/cp2k/colvar.py` 等下一 commit

**测试**:全量 unit + integration;`test/unit/engines/test_facade.py:72-73` 类型
断言专项验证

### Commit 2 — D8 ConstraintRun

- 新增 `engines.models.ConstraintRun`(纯新增)
- 新增 `engines.cp2k.read_constraint_run(directory)` facade
- 把 `ColvarMDInfo` 改为 `engines.models.ColvarMDInfo = ConstraintRun` alias
- `utils.formats.cp2k.colvar.ColvarMDInfo` 改为 `from engines.models import
  ConstraintRun as ColvarMDInfo`(过渡期 utils → engines re-export 仅在 alias
  层,**docstring 必须**显式说明这是 Phase 5 业务迁移完成前的过渡)

**注意**:Commit 2 在 utils/formats/cp2k/colvar.py 引入 utils → engines runtime
import。这**短期**违反 R2,但是设计文档接受的过渡期成本。Phase 5 业务迁移结束后
**必须**把这条 re-export 删除。

**测试**:全量 unit + integration;`ConstraintRun.target_series_au` property 公式
专项测试

### Commit 3 — D7 CellSpec

- 新增 `engines.models.CellSpec`
- 新增 `engines.cp2k.read_cell(path)` facade(内部委托 `parse_abc_from_md_inp` /
  `parse_abc_from_restart`)

**测试**:`CellSpec.abc_ang` / `is_orthorhombic` 派生 property 专项测试

### Commit 4 — D11 CenterPotentialScalarFrame

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
| **import 边界扫描** | `rg "from .*engines\|^import.*engines" src/md_analysis/utils --type py \| grep -v TYPE_CHECKING` | 仅 `utils/formats/cp2k/colvar.py` 的 ColvarMDInfo re-export 在过渡期合法;其他全空 |
| **data_example 回归** | `data_example/potential/` + `data_example/bader/` + `data_example/sg/` | integration 全过(44 passed) + CSV 列值 byte-equal(diff fixture 输出) |
| **业务 caller 不破** | full unit + integration | 729 + 44 基线不变 |

每笔 commit 单独跑测试,任一回归即修复;不允许"积压后统一处理"。

---

## 11. Phase 5 业务迁移影响面(基于 Phase 2 catalog 重算)

Phase 4 完成后,业务层 import 形态:

| 业务节点 | Phase 4 后状态 | Phase 5 任务 |
|---|---|---|
| `water._common._parse_abc_from_md_inp` | 仍直读 `utils.formats.cp2k.cell` | 切到 `engines.cp2k.read_cell(...)` → `CellSpec.abc_ang` |
| `electrochemical.potential.CenterPotential.parse_md_out_fermi` (line 447) | 仍直读 | 选项 a:切到 `engines.cp2k.read_center_potential_scalar_frame(...)`(同时迁离 dict 接口);选项 b:保留 dict 接口直读,只在新增分支切 facade(Phase 5 拍板) |
| `electrochemical.charge.Bader.*` (7 处) | 保持 | Phase 4 不动,留待 Phase 8 charge engines facade 落地 |
| `enhanced_sampling.slowgrowth.SlowGrowth.{from_paths,from_directory}` (2 处) | 仍直读 `ColvarMDInfo`(已是 `ConstraintRun` alias) | 切到 `from engines.models import ConstraintRun` |
| `enhanced_sampling.constrained_ti.workflow.standalone_diagnostics` (1 处) | 同上 | 同上 |
| `enhanced_sampling.constrained_ti.correction._get_electrode_area` (1 处) | 仍直读 `load_bader_atoms` | Phase 5 cleanup:可独立提"从 POSCAR 读 cell"小 helper(charge 仍留 utils);或同 Phase 8 charge facade 一起迁 |

Phase 5 commit 边界初步设想:
- 1 笔:水/cli cell 迁移(D7 CellSpec 消费方)
- 1 笔:SG/TI/cli composite 迁移(D8 ConstraintRun 消费方)
- 1 笔:potential Fermi-only / dict 接口迁移(D11 + dict 历史决议同步处理)
- charge / correction 业务迁移留 Phase 8(VASP 阶段)

入口层(cli / agent / scripts)迁移**不在 Phase 5**;按 `overall_reconstruction_plan.md`
Phase 6 统一处理。

---

## 设计文档完结

待办:
1. user + codex 三次审本文件
2. 通过后进 Phase 4(实施),按 §9 4 笔 commit 推进
3. 任何 §6.3 / §9 子步骤回归,先停下找根因再继续;不允许在回归未消失时进下一 commit
