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
- 现有 `read_constraint_metadata` / `read_lambda_series` 不动

`ColvarMDInfo` 在 Phase 4 Commit 2 物理迁移到 `engines.models`(作为
`ColvarMDInfo = ConstraintRun` 别名),`utils.formats.cp2k.colvar` 删除定义、**不**
保留任何 alias re-export(R2 硬约束,见 §6.4 同 commit 切 caller)。

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
| `scripts/TIGen.py:19` | `from ..utils.formats.cp2k.colvar import parse_colvar_restart` + 用 `restart.colvars.primary.target_au` 等 | 切到 `from ..engines.cp2k import read_constraint_metadata`;调用方式从 `parse_colvar_restart(path)` 改为 `read_constraint_metadata(directory)`(canonical) — 或保留 parser 直调 + 通过 `_cp2k_raw_to_constraint_metadata(raw)` 转 canonical 后再访问 `.colvars.primary` |
| `test/unit/utils/test_slowgrowth_parser.py:11` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarInfo, ColvarMDInfo, ColvarParseError, ColvarRestart, ConstraintInfo, LagrangeMultLog, compute_target_series, parse_colvar_restart, parse_lagrange_mult_log)` | 拆分:`ColvarInfo` / `ColvarMDInfo` / `ColvarRestart` / `ConstraintInfo` / `LagrangeMultLog` 从 `md_analysis.engines.models` import;`ColvarParseError` / `compute_target_series` / `parse_colvar_restart` / `parse_lagrange_mult_log` 保留 utils path |
| `test/unit/utils/test_slowgrowth_parser.py:367` | `from md_analysis.utils.formats.cp2k.colvar import _parse_fixed_atoms_list` | 不切(private utils helper) |
| `test/unit/scripts/test_ti_gen.py:22` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarRestart, ColvarInfo, ConstraintInfo)` | 全部切到 `from md_analysis.engines.models import (...)` |
| `test/unit/engines/test_facade.py:67` | `from md_analysis.utils.formats.cp2k.colvar import (ColvarRestart, LagrangeMultLog)` | 切到 `from md_analysis.engines.models import (ColvarRestart, LagrangeMultLog)`;`test_facade.py:72-73` 的 `assert ColvarRestart is ConstraintMetadata` 仍成立(canonical alias 物理在 engines.models) |
| `engines/models.py:41`(re-export) | `from ..utils.formats.cp2k.colvar import ConstraintMetadata, LambdaSeries` | **删除该 re-export**;改为本地物理定义 |
| `engines/protocols.py:18` | `from .models import ConstraintMetadata, LambdaSeries` | 不动 |
| `engines/cp2k.py:33` | `from ..utils.formats.cp2k.colvar import (parse_colvar_restart, parse_lagrange_mult_log)` | 不动 import,但调用方式改:parser 返 raw → 走 `_cp2k_raw_to_constraint_metadata` / `_cp2k_raw_to_lambda_series` 转换 helper |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:18` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | 同 commit 切到 `from ...engines.models import ColvarMDInfo`(`ColvarMDInfo = ConstraintRun` alias 物理在 engines.models) |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:179` | 同上 | 同上 |
| `enhanced_sampling/constrained_ti/workflow.py:681` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | 同上 |
| `cli/_enhanced_sampling.py:44` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo")` | 同 commit 切到 `lazy_import("md_analysis.engines.models", "ColvarMDInfo")` |
| `cli/_scripts.py:267` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart")` | 同 commit 切到 `lazy_import("md_analysis.engines.cp2k", "read_constraint_metadata")`(走 facade 拿 canonical)或保留 utils path(仅消费 raw)— Phase 4 实施时选定后写在 commit message |
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

Phase 4 总共 **4 笔 commit**(每笔不超过当前已知最大 commit 规模):

### Commit 1 — D10 物理迁移(最高风险,先做)

- 引入 raw 类型(`Cp2kConstraintInfoRaw` / `Cp2kColvarInfoRaw` /
  `Cp2kConstraintMetadataRaw` / `Cp2kLambdaSeriesRaw`)在 utils.formats.cp2k.colvar
- 物理迁移 canonical 类型(`ConstraintInfo` / `ColvarInfo` / `ConstraintMetadata` /
  `LambdaSeries`)到 `engines.models`,**含嵌套 neutral 类型**(codex Round 5
  HIGH 2)
- 物理迁移 Phase 5b 历史 alias(`ColvarRestart` / `LagrangeMultLog`)到
  `engines.models`,utils.formats.cp2k.colvar **不**持有它们
- parser 函数签名变更(返 raw)+ raw → canonical 转换 helper(嵌套四级)
- §6.4 表中 `ColvarRestart` / `LagrangeMultLog` / `ColvarInfo` / `ConstraintInfo`
  import 站点(test/unit/engines/test_facade.py + test/unit/utils/test_slowgrowth_parser.py
  + test/unit/scripts/test_ti_gen.py)**同 commit 切**到 engines.models
- `ColvarMDInfo` 留 utils.formats.cp2k.colvar 等下一 commit

**门槛**:utils.formats.cp2k.colvar **无** runtime engines import(R2 全程满足)。

**测试**:全量 unit + integration;`test/unit/engines/test_facade.py:72-73` 类型
断言专项验证(`ColvarRestart is ConstraintMetadata` 仍成立);`rg "engines"
src/md_analysis/utils --type py | grep -v TYPE_CHECKING | grep -v "^.*:#" |
grep -v "^.*:\".*\""` 仅在 docstring/注释出现

### Commit 2 — D8 ConstraintRun + ColvarMDInfo 同步切

- 新增 `engines.models.ConstraintRun`(纯新增)
- 新增 `engines.cp2k.read_constraint_run(directory)` facade
- 物理迁移 `ColvarMDInfo` **别名**到 `engines.models`(`ColvarMDInfo =
  ConstraintRun`)
- **删除** `utils/formats/cp2k/colvar.py` 中的 `ColvarMDInfo` 定义;**不**在 utils
  层做 re-export
- §6.4 表中 5 处 `ColvarMDInfo` import 站点(SG×2 + TI×1 + CLI×1 lazy_import +
  agent 是字符串保留)**同 commit 机械切**到 engines.models
- cli `_scripts.py:267` 的 lazy_import 字符串按 §6.4 选定切换方向

**门槛**:utils.formats.cp2k.colvar **无** runtime engines import(R2 全程满足);
utils 层**不**持有 engines alias。

**测试**:全量 unit + integration(729 + 44 基线);`ConstraintRun.target_series_au`
property 公式专项测试与 utils 旧 `ColvarMDInfo.target_series_au` numerically equal

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
| **import 边界扫描** | `rg "from .*engines\|^import.*engines" src/md_analysis/utils --type py \| grep -v TYPE_CHECKING` | **全空**(R2 硬约束;过渡期 re-export 已删除,见 §6.3 Step 5 + §6.4 同 commit 切 caller 策略) |
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
| `enhanced_sampling.slowgrowth.SlowGrowth.{from_paths,from_directory}` (2 处) | **已 import** `from ...engines.models import ColvarMDInfo`(Phase 4 Commit 2 同步切完;`ColvarMDInfo` 是 `ConstraintRun` 别名,物理定义在 engines.models) | **命名清理**:把 `ColvarMDInfo` import 名换成 canonical `ConstraintRun`,删除别名依赖 |
| `enhanced_sampling.constrained_ti.workflow.standalone_diagnostics` (1 处) | 同上 | 同上 |
| `enhanced_sampling.constrained_ti.correction._get_electrode_area` (1 处) | 仍直读 `load_bader_atoms` | Phase 5 cleanup:可独立提"从 POSCAR 读 cell"小 helper(charge 仍留 utils);或同 Phase 8 charge facade 一起迁 |

Phase 5 commit 边界初步设想:
- 1 笔:水/cli cell 迁移(D7 CellSpec 消费方)
- 1 笔:SG/TI/cli 命名清理(`ColvarMDInfo` → `ConstraintRun`;Phase 4 已完成路径
  切换,Phase 5 只做 canonical 名字统一)
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
