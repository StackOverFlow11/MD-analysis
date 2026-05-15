# Phase 2 — 业务层数据需求 Catalog

> 用途：从顶层业务反推所需 engine-neutral 数据形态，作为 Phase 3 engines 数据接口设
> 计的**单一输入源**。本文件不指定 Phase 3 的类名、字段名、facade 名，也不假设 VASP
> 端能与 CP2K 单位无损对齐——只列出当前业务的实际消费证据与判断依据。

## 文档约定

- **当前调用点证据**：grep + Read 实证的 import 行 + 该 import 在业务函数体内被消费
  的字段/方法（精确到文件:行）。
- **业务实际需要的数据**：剥离文件格式细节后，业务工作流真正消费的物理量、形状、单
  位。
- **engines 当前覆盖现状**：对照 `engines.models` 的 4 个 frozen dataclass
  (`ConstraintMetadata` / `LambdaSeries` / `PotentialFrame` / `FermiRecord`) 与
  `engines.cp2k` 的 5 个 facade，写明覆盖到哪、缺到哪、为什么不够。
- **Phase 3 候选标注**（不预判类名/字段/facade，只记决策点）：
  - `NEW-FACADE`：业务此处需要一个 engine-neutral facade，engines 目前没有任何对应
    入口。
  - `EXTEND-MODEL`：业务已经走 engines，但需要扩展现有 dataclass 字段或派生方法。
  - `ENGINE-SCOPE`：业务依赖的数据是否归 engines 抽象本身就有歧义（如 Bader 输入来
    源 vs 业务语义），Phase 3 需要先拍板归属。
  - `STAY-AS-IS`：业务直读 utils.formats 在 Phase 2 范围内未发现明确升级理由，留作
    Phase 5 业务迁移决策点。

## engines 现有清单（事实陈述，非设计建议）

| dataclass / facade | 物理语义 | 当前 CP2K-only 实现 |
|---|---|---|
| `ConstraintMetadata` | 约束 MD 元数据：target / growth / timestep / fixed atoms / cell | `utils.formats.cp2k.colvar.parse_colvar_restart` |
| `LambdaSeries` | λ(t) 时间序列：shake / rattle / n_steps / n_constraints | `utils.formats.cp2k.colvar.parse_lagrange_mult_log` |
| `PotentialFrame` | 单帧 Hartree 势：step / time_fs / cube_path / header / values / fermi_raw / atoms | `engines.cp2k.read_continuous_potential_frames` / `engines.cp2k.read_distributed_potential_frames` |
| `FermiRecord` | 单条 Fermi 记录：step / time_fs / fermi_raw（+ legacy dict 桥接） | `utils.formats.cp2k.stdout.parse_md_out_fermi` |
| `read_constraint_metadata(directory)` | 目录 → ConstraintMetadata | CP2K |
| `read_lambda_series(directory)` | 目录 → LambdaSeries | CP2K |
| `read_fermi_series(md_out_path)` | md.out → list[FermiRecord] | CP2K |
| `read_continuous_potential_frames(...)` | 单目录 cube + md.out → list[PotentialFrame]（mode A） | CP2K |
| `read_distributed_potential_frames(...)` | `potential_t*_i*/` 子目录 → list[PotentialFrame]（mode B） | CP2K |

`ConstraintMDParser` Protocol（registry: `_REGISTRY`）+ CP2K 实现已注册；VASP 占
位**未**注册（确保 `infer_parser` 不会误派发到 stub）。

---

## 业务模块 1：`water/`

### 1.1 当前数据来源

| 文件:行 | import | 经 engines？ |
|---|---|---|
| `water/WaterAnalysis/_common.py:44` | `from ...utils.formats.cp2k.cell import parse_abc_from_md_inp as _parse_abc_from_md_inp` | ⚠️ 绕过 |

调用点（`_common.py` 内 `_iter_trajectory` / `_resolve_cell_abc` 调用链）：

```python
elif md_inp_path is not None:
    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)
```

`WaterDensity.py` / `AdWaterOrientation.py` 通过 `from ._common import _parse_abc_from_md_inp`
跟着用同一个 helper。

### 1.2 业务实际需要的数据

- 单一三元组 `(a, b, c)`：体系**正交**晶胞的轴长（Å），用于把 xyz 轨迹的笛卡尔坐标
  映射回分数坐标 + 给 ase.Atoms 套 cell + PBC。
- **不**需要完整 cell 矩阵；**不**需要 cell 随时间变化（NVT/NVE，cell 静态）。
- 单一可选 fallback：调用方可以直接传 `cell_abc=(a,b,c)`，没传时才从文件解析。

### 1.3 engines 当前覆盖现状

`ConstraintMetadata` 内部其实**已经携带 cell ABC**（`parse_abc_from_restart` 是
`utils.formats.cp2k.colvar` 内部 import），但：

- 该 cell 字段只能从 `.restart` 文件获得，不能从 `md.inp` 获得。
- `engines.cp2k` 没有任何"只读 cell 的"facade；`read_constraint_metadata` 返回的对
  象虽然含 cell，但要求传入"约束 MD 目录"——纯水分析的运行目录可能没有约束 MD 文
  件，强行走它会引入无关错误路径。
- water 当前的需求是"给我这个体系的 (a,b,c)"，按现有 engines 抽象，唯一 fit 的位
  置是 `read_constraint_metadata`，但语义不匹配。

### 1.4 Phase 3 候选

- **NEW-FACADE**（候选）：engines 是否需要一个"仅读 cell"的 facade（来源可以是
  `.restart`、`md.inp`、未来 VASP 的 POSCAR），返回某个 engine-neutral 形态的 cell
  描述。
- **单位/语义风险**（Phase 3 拍板时需要考虑）：
  - 当前 water 只支持正交 cell（`utils/structure/__init__.py:9` + `target_structure.md`
    "三基矢正交"前提），但 ASE Atoms 的 `.cell` 是 3×3 矩阵；engines facade 是否应
    该返回 3×3 还是 `(a,b,c)` 元组，影响下游业务签名。
  - VASP 的 POSCAR 同样可能是非正交，下游一旦扩展非正交也要求 engines facade 表达
    完整 lattice，而不仅 `(a,b,c)`。

---

## 业务模块 2：`electrochemical/potential/`

### 2.1 当前数据来源

| 文件:行 | import | 经 engines？ |
|---|---|---|
| `electrochemical/potential/_frame_source.py:31-35` | `from ...engines.cp2k import read_continuous_potential_frames, read_distributed_potential_frames` + `from ...engines.models import PotentialFrame` | ✅ thin wrapper，全经 engines |
| `electrochemical/potential/CenterPotential.py:34` | `from ...utils.formats.common.cube import (_float, discover_cube_files, extract_step_from_cube_filename, read_cube_header_and_values, slab_average_potential_ev)` | ⚠️ 绕过（直读 cube 工具集） |
| `electrochemical/potential/CenterPotential.py:447` | `from ...utils.formats.cp2k.stdout import parse_md_out_fermi`（函数体内 lazy） | ⚠️ 绕过（直读 md.out） |
| `electrochemical/potential/PhiZProfile.py:21` | `from ...utils.formats.common.cube import (CubeHeader, _float, discover_cube_files, extract_step_from_cube_filename, plane_avg_phi_z_ev, read_cube_atoms, read_cube_header_and_values, z_coords_ang)` | ⚠️ 绕过（直读 cube 工具集） |

`_frame_source.py` 本身已经是 thin forwarding wrapper（Phase 7b2 退化），不含业务
逻辑。

### 2.2 业务实际需要的数据

按 `electrochemical/potential/CLAUDE.md` 与函数体证据：

**整帧序列消费**（已 engines 化）：
- `list[PotentialFrame]`，每帧含 `step / time_fs / cube_path / header / values /
  fermi_raw / atoms`，覆盖整个 cSHE 计算流程（slab-averaged Hartree → φ_center → U
  vs SHE）。两种发现模式（continuous / distributed）均已统一在 `PotentialFrame` 后
  端。

**单工具消费**（仍走 utils.formats）：
- `slab_average_potential_ev(...)`：给定 `header + values + slab range`，返回
  `(phi_center_ev, info)`，业务用 `info["phi_z_std_ev"]` 做厚度灵敏度。
- `plane_avg_phi_z_ev(...)`：返回沿 c 轴的 φ(z) 1D profile，用于 `PhiZProfile.py`
  做跨帧 overlay。
- `z_coords_ang(...)`：cube grid → z 坐标。
- `read_cube_atoms(path, header)`：从 cube 头解析 atoms（独立于 xyz）。
- `discover_cube_files(...)`、`extract_step_from_cube_filename(...)`：cube 文件发
  现 + 文件名解析。**但这些已经被 `read_continuous_potential_frames` 在内部用过**；
  `CenterPotential.py` 在 frame 发现之外另有"fermi-only / distributed"分支会重新
  调用。
- `parse_md_out_fermi(...)`：连续模式下纯 Fermi 时序（不附 cube）；用于 213
  Fermi-only 路径。

### 2.3 engines 当前覆盖现状

- 整帧路径 ✅ 已经走 `PotentialFrame` + 两个 facade，没有破口。
- **工具集消费**是 gap：
  - `slab_average_potential_ev / plane_avg_phi_z_ev / z_coords_ang / read_cube_atoms`
    几个 cube 后处理函数纯算法（输入 header + values 数组），跟 engine 无关，**留
    在 utils.formats.common.cube 才是符合边界的**——业务直接 import 它们不算"违反
    engines 边界"，是数据-后处理工具的正常调用。
  - `parse_md_out_fermi` 是 CP2K 文件解析层，engines 已经在 `read_fermi_series`
    facade 里包装了它（`FermiRecord.from_legacy_dict`）。但 `CenterPotential.py:447`
    仍直读，理由：`CenterPotential.py:455-457` 还用 dict-style `r["step"]` 访问
    （`engines/CLAUDE.md` 行 73-78 已记录这个**有意保留的 dict 接口**）。

### 2.4 Phase 3 候选

- **STAY-AS-IS**：`common.cube` 的纯算法工具（`slab_average_potential_ev` /
  `plane_avg_phi_z_ev` / `z_coords_ang` / `read_cube_atoms`）按 Phase 1 边界本来就
  应该住在 `utils.formats.common.cube`，业务直接 import 是合法路径，不需要 engines
  facade。
- **EXTEND-MODEL**（候选）：`PotentialFrame` 目前只承载"已发现帧"。如果 Phase 3 决
  定 213 Fermi-only / distributed 路径也用统一 frame 模型，可以考虑引入"轻量帧"
  （只有 step/time/fermi，无 cube），让 `CenterPotential.py` 不再走 dict-style
  parse_md_out_fermi。
  - **决策注意**：`engines/CLAUDE.md:73-78` 明确"dict 接口有意保留，强迫迁移会牵
    涉 caller"。Phase 3 拍板前需要先决定"是否在 Phase 4 一起迁移 CenterPotential
    的 dict access"，否则 EXTEND-MODEL 不能落地。

---

## 业务模块 3：`electrochemical/charge/`

> ⚠️ Phase 2 **不**预判 charge 数据应该归 `engines.cp2k`、`engines.vasp` 还是新
> 的 `engines.charge`。只记录：当前输入文件链路是 VASP 后处理产物 (POSCAR / ACF /
> POTCAR)，业务消费的语义是 charge frame / trajectory。归属问题留 Phase 3 拍板。

### 3.1 当前数据来源

| 文件:行 | import | 经 engines？ |
|---|---|---|
| `electrochemical/charge/Bader/BaderData.py:12` | `from ....utils.formats.bader.acf import load_bader_atoms` | ⚠️ 绕过 |
| `electrochemical/charge/Bader/AtomCharges.py:15` | 同上 | ⚠️ 绕过 |
| `electrochemical/charge/Bader/SurfaceCharge.py:14` | 同上 | ⚠️ 绕过 |

`load_bader_atoms` 在该子包内被反复调用 7 次（`BaderData.py:115` /
`AtomCharges.py:221, 477` / `SurfaceCharge.py:345, 467` 等），均是"逐帧加载一个目
录的 POSCAR + ACF + POTCAR → 带电荷 Atoms"。

### 3.2 业务实际需要的数据

**单帧**：
- `ase.Atoms`（结构 + cell），并在 `atoms.arrays` 上挂两份数组：
  - `bader_charge` — 原始 Bader 电子数（POSCAR 序）
  - `bader_net_charge` — `ZVAL - bader_charge`，正号 = 失电子（POSCAR 序）

**轨迹**（`BaderTrajectoryData`，已是该子包的 frozen dataclass）：
- `steps: (n_frames,) int` / `times: (n_frames,) int(fs)` / `atom_indices_xyz:
  (n_atoms,) int` / `net_charges: (n_frames, n_atoms) float`，**已 remap 到 XYZ
  序**。

业务逻辑层（`SurfaceCharge.py` 的 σ 计算、`AtomCharges.py` 的指定原子追踪、
`counterion_charge_analysis` 的 counterion 自动检测）全部基于这两个对象工作。

### 3.3 engines 当前覆盖现状

- engines.models 目前**完全没有** charge frame / charge trajectory 任何对应
  dataclass；4 个现有 dataclass 都不覆盖这类数据。
- `engines.cp2k` 没有任何 charge facade（CP2K 本身的 Mulliken 也未实现）。
- `engines.vasp` 是 placeholder，不含任何 ACF/POTCAR 解析。
- 业务自己在 `electrochemical/charge/Bader/BaderData.py` 里维护了一个**业务层**的
  trajectory dataclass `BaderTrajectoryData`，自己处理 POSCAR↔XYZ remap，自己用
  `_sorted_frame_dirs` + `_frame_utils._extract_step_and_time` 做帧发现。

### 3.4 Phase 3 候选

- **ENGINE-SCOPE**（必须先拍板）：
  - 输入文件链路：POSCAR + ACF.dat + POTCAR——全部是 VASP 后处理产物。
  - 业务语义：charge frame / charge trajectory，与 CP2K/VASP 选哪个引擎无关。
  - 候选方案（**不在 Phase 2 拍板**）：
    1. 归 `engines.vasp`：尊重数据物理来源，与未来 `engines.vasp.read_locpot` 同
       门面。
    2. 新增 `engines.charge` 子门面：承认 Bader 是"按引擎切分电荷输出"模式的第一
       例，未来 CP2K 的 Mulliken / Gaussian 的 Hirshfeld 也按此扩展。
    3. 保留业务直读：业务自己持有 `BaderTrajectoryData`，engines 不抽象 charge。
       理由：Bader frame 几何严格绑定 VASP POSCAR 索引（IndexMapper 已是业务层的
       工程产物），强行 engine-neutral 化收益小于成本。
  - Phase 2 不选；Phase 3 设计阶段对照三方案的 `utils.formats.vasp.*` 占位实现进
    度 + `long_term.md` "多引擎适配策略"再决定。

- **横向影响**（不论 Phase 3 选哪条）：
  - charge 的 frame discovery (`_frame_utils._extract_step_and_time` +
    `_sorted_frame_dirs`) 与 `engines.cp2k.read_distributed_potential_frames` 在
    `potential_t*_i*` 上用的发现器**都委托到** `utils/io/_frame_discovery.py` 这
    个底层 frame-dir helper（step/time 名称解析共享）；缺的是 engines-level 的
    **charge/potential frame abstraction**——potential 走 `PotentialFrame`，
    charge 走业务层 `BaderTrajectoryData`，两者之间没有公共"engine frame"模型。
    Phase 3 若收敛 "engine frame" 概念，可与上面 ENGINE-SCOPE 一起拍板。

---

## 业务模块 4：`electrochemical/calibration/`

### 4.1 当前数据来源

`rg "utils\.formats|engines" src/md_analysis/electrochemical/calibration` → **空**。

`calibration/` 完全不依赖 `utils.formats.*` 或 `engines.*`：
- 输入：业务上游传入的 `(φ, σ)` 数据点（CSV 或手动）
- 内部：纯数值拟合（线性、多项式、样条、微分电容）
- 输出：mapper 对象 + JSON 配置

### 4.2 业务实际需要的数据

`CalibrationData` / `Mapper` 完全 self-contained，不消费任何文件解析层产物。

### 4.3 engines 当前覆盖现状

不适用。

### 4.4 Phase 3 候选

- **STAY-AS-IS**（无变化）：calibration 在 Phase 1-7 全程不应被 engines 重构影响。
  唯一与 engines 的交叉是 `enhanced_sampling/constrained_ti/correction.py` 把
  calibration `mapper.predict(σ)` 接到 Bader σ→Φ 链路，这部分歧义已在 charge 节
  和 constrained_ti 节单独讨论。

---

## 业务模块 5：`enhanced_sampling/slowgrowth/`

### 5.1 当前数据来源

| 文件:行 | import | 经 engines？ |
|---|---|---|
| `enhanced_sampling/slowgrowth/SlowGrowth.py:18` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo` | ⚠️ 绕过（顶层） |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:174` | `from ...engines import (ConstraintMDParser, infer_parser, resolve_parser)` | ✅ 经 engines（Protocol） |
| `enhanced_sampling/slowgrowth/SlowGrowth.py:179` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo`（`from_directory` 内 lazy） | ⚠️ 绕过 |

### 5.2 业务实际需要的数据

`SlowgrowthFull.from_md_info(md_info: ColvarMDInfo, ...)`（`SlowGrowth.py:99-139`）
完整列出业务消费字段：

- `md_info.steps`（绝对 step 数组）
- `md_info.times_fs`（绝对时间数组）
- `md_info.target_series_au(colvar_id)`（重建的 CV 目标序列）
- `md_info.lagrange.collective_shake`（λ(t) shake 序列）
- `md_info.restart.timestep_fs`（标量）
- `md_info.restart.colvars[colvar_id].target_growth_au`（标量，per-a.u.-time）

业务还保留 `md_info` 引用挂在 `SlowgrowthFull.md_info` 上，方便下游回溯 restart 元
数据（cell、fixed atoms 等）。

### 5.3 engines 当前覆盖现状

- `ConstraintMetadata` 已经包含 `colvars / timestep_fs / step_start /
  time_start_fs / fixed_atoms / cell`。
- `LambdaSeries` 已经包含 `shake / rattle / n_steps / n_constraints`。
- engines.protocols 已经定义 `ConstraintMDParser.parse_metadata(directory)` 与
  `.parse_lambda_series(directory)`。
- **缺口**：业务消费的不是 "ConstraintMetadata + LambdaSeries" 两个独立对象，而是
  它们组合后**派生**的 `target_series_au(...)`、`times_fs`、`steps` 这些**对齐
  到 absolute step axis 的派生序列**。当前这套派生在 `ColvarMDInfo` 上以方法/
  property 形式提供（`cp2k/colvar.py:124-176`），定义在 utils 层，engines.models
  没有对应物。

### 5.4 Phase 3 候选

- **EXTEND-MODEL**（候选，**不**在 Phase 2 定形态）：业务需要"约束 MD 运行数据的
  时间对齐组合视图"，当前的 `ColvarMDInfo` 是它的 CP2K-specific 实现。
  - **是 Phase 3 候选的理由**：
    - 业务签名（`from_md_info(md_info)`）就在告诉读者"我需要一个统一对象"，单独传
      `(metadata, lambda_series)` 两元组对消费方不方便。
    - `from_directory(directory, parser="auto")` 已经做出 engine-neutral 接入尝
      试，但内部仍 `ColvarMDInfo(restart=parser.parse_metadata(...),
      lagrange=parser.parse_lambda_series(...))`，**绑死了 CP2K-specific
      dataclass 名**。如果 VASP 接入时构造一个"看起来等价"的 ColvarMDInfo，会让
      utils.formats.cp2k 类型在 VASP 路径上被引用——违反 engine 边界。
  - **单位/语义风险**（Phase 3 拍板时务必考虑，**Phase 2 不预设**）：
    - `target_series_au` 当前用 `dt_au = timestep_fs / AU_TIME_TO_FS` 做 fs→a.u.
      换算 + `target_au + (k - step_start) * target_growth_au * dt_au` 重建。
      VASP 的 SHAKE/RATTLE 是否使用相同的 "growth per a.u.-time" 概念，需要 VASP
      原始格式调研后才能确认；不能默认 VASP 侧能无损映射到同一公式。
    - `colvars.primary` / `ConstraintMetadata.fixed_atoms` 是 CP2K block 解析直接
      得到的，VASP 端的 constraint metadata 表达形式不同。Phase 3 设计这层组合视
      图时，需要先看 VASP REPORT 的实际信息密度。
    - "派生序列"的所有权：放在 dataclass 上（property/方法） vs 由 facade 在读取
      时一次性算好放回 dataclass 字段——这是 Phase 3 设计取舍，会影响业务 caller
      代码风格。
  - **当前 engines 不够的根本原因**：现有 facade 只提供"输入态"（metadata / lambda
    series），没有提供"对齐到 absolute step axis 的派生量"；业务为了不重复实现该
    派生，只能 import utils.formats 层的 `ColvarMDInfo`。

---

## 业务模块 6：`enhanced_sampling/constrained_ti/`

### 6.1 当前数据来源

| 文件:行 | import | 经 engines？ |
|---|---|---|
| `enhanced_sampling/constrained_ti/io.py:16-21` | `from ...engines import (ConstraintMDParser, infer_parser, resolve_parser)` + `from ...engines.protocols import _REGISTRY` | ✅ 经 engines |
| `enhanced_sampling/constrained_ti/models.py:17-18` | `from ...engines.protocols import ConstraintMDParser` + `from ...engines.models import ConstraintMetadata` (TYPE_CHECKING) | ✅ 经 engines |
| `enhanced_sampling/constrained_ti/workflow.py:681` | `from ...utils.formats.cp2k.colvar import ColvarMDInfo`（函数体内 lazy） | ⚠️ 绕过 |
| `enhanced_sampling/constrained_ti/correction.py:63` | `from ...utils.formats.bader.acf import load_bader_atoms`（函数体内 lazy） | ⚠️ 绕过 |

`io.py` + `models.py`（TIPointDefinition / TIReport）已经完全 engine-agnostic：
`discover_ti_points(parser="auto", dir_filter=None)` 通过 `ConstraintMDParser`
Protocol + registry sniffing 工作；`TIPointDefinition.metadata` 是
`ConstraintMetadata`，业务从此得到 `xi`、`target_au`、`time_start_fs` 等。

### 6.2 业务实际需要的数据

`io.py + workflow.analyze_ti` 主链路（**已 engines 化**）：
- `TIPointDefinition` 每点 → `metadata: ConstraintMetadata` + `λ(t):
  LambdaSeries` → 四步诊断（ACF / FP block / running avg / Geweke） → `TIReport`
  + 梯形积分自由能。

`workflow.standalone_diagnostics`（`workflow.py:681` lazy，**未经 engines**）：
- 输入是单点的 restart + log 路径（不一定在 ti_target 目录里），用 `ColvarMDInfo`
  组合 metadata + lambda series + 派生 `time_start_fs` / `target_au`，再走 ACF /
  FP / Geweke 诊断（不积分自由能）。
- 业务消费 `md_info.lagrange.collective_shake / md_info.restart.timestep_fs /
  md_info.restart.colvars[colvar_id].target_au / md_info.restart.time_start_fs`。

`correction.py._get_electrode_area`（**未经 engines**）：
- 输入：第一帧 Bader 目录的 `(POSCAR, ACF, POTCAR)`。
- 业务消费：`atoms.cell.array`（用 `np.cross` 计算电极面积）；电荷数据本身在该函数
  里**不**消费。`load_bader_atoms` 在这里被用作"读 POSCAR + 顺手挂电荷"，但
  `_get_electrode_area` 只用 cell。
- 整个 `correction.py.compute_constant_potential_correction` 还间接消费 σ 时序
  （来自 `electrochemical.charge.Bader.trajectory_surface_charge`）+ φ
  （`calibration.mapper.predict(σ)`），但这两个数据源已经在各自业务模块里完成，
  correction 自己只把它们组合起来。

### 6.3 engines 当前覆盖现状

- `analyze_ti` 主链路 ✅ 没有破口。
- `standalone_diagnostics` 跟 slowgrowth 共用同一个 `ColvarMDInfo` 派生视图问题
  （见 5.4，**不重复列**）。
- `_get_electrode_area` 跟 charge 子包共用同一个 Bader 输入问题（见 3.4，**不重
  复列**）。额外的小观察：此处实际只需要 cell，不需要电荷数组——目前为复用
  `load_bader_atoms` 顺带把电荷也读进来，是 acceptable 但不是最小依赖。

### 6.4 Phase 3 候选

- **EXTEND-MODEL**（与 slowgrowth 同源，见 5.4）：`standalone_diagnostics` 用
  `ColvarMDInfo.from_paths(restart, log)` 是 slowgrowth `from_paths` 的同形态。
  Phase 3 处理 5.4 时**一定**会同时影响 6.x；这两个业务节点共享一个 gap。
- **ENGINE-SCOPE**（与 charge 同源，见 3.4）：`_get_electrode_area` 走的是
  `load_bader_atoms`，归属随 charge 节决策。
- **MINIMAL-DEPS 观察**（Phase 5 时可顺手清）：`_get_electrode_area` 只需要 cell，
  不需要电荷。当 charge 数据接入 engines 后，可以让 correction.py 用更轻量的
  "Atoms + cell" 路径替代 `load_bader_atoms`。这是业务侧的清理任务，不属于 engines
  接口设计本身。

---

## 入口层绕过 engines 的现状（cli + agent + scripts）

> 单独成节。这些不是业务模块的数据消费者，但它们直读 utils.formats 的事实，能反过
> 来暴露"Phase 3 必须提供给上层的 facade 形态"。Phase 6 才会迁移它们到 engines。
> Phase 2 只做记录。

### 7.1 当前数据来源

| 文件:行 | import | 用途 |
|---|---|---|
| `cli/_params.py:315` | `from ..utils.formats.cp2k.cell import (CellParseError, parse_abc_from_md_inp, parse_abc_from_restart)` | CLI 交互式 cell 参数采集（`CellAbcParam.collect()`） |
| `cli/_scripts.py:267` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart")` | CLI 410-440 系列脚本菜单展示 SG CV 信息 |
| `cli/_enhanced_sampling.py:44` | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo")` | CLI 30x SG 菜单展示轨迹元数据 |
| `agent/_tasks_legacy.py:677` | exception FQN `"md_analysis.utils.formats.cp2k.colvar.ColvarParseError"` | agent 任务契约的异常映射 |
| `agent/_tasks_ti.py:228` | exception FQN `"md_analysis.engines.protocols.ParserInferenceError"` | agent TI 任务契约的异常映射 ✅ |
| `scripts/TIGen.py:16` | `from ..engines.models import ConstraintMetadata` | TI 工作目录生成的 metadata 类型注解 ✅ |
| `scripts/TIGen.py:19` | `from ..utils.formats.cp2k.colvar import parse_colvar_restart` | TI 工作目录生成时读 SG restart 拿 colvar 信息 |

### 7.2 入口层暴露的 Phase 3 facade 需求

- **cell 解析**（`cli/_params.py` + 5.x water）：入口层需要"按用户提供的 restart/
  md.inp 路径任选一个，返回 cell"——这是 water 在业务层也需要的 facade（见 1.4
  `NEW-FACADE` 候选）。如果 Phase 3 决定加 cell facade，cli 也是消费者。
- **CV 序列预览**（`cli/_scripts.py` + `cli/_enhanced_sampling.py`）：CLI 在
  TIGen / SG 菜单里要展示"轨迹有多少步 / CV 范围 / 时间步"，**直接复用
  ColvarMDInfo / parse_colvar_restart 的派生量**。这与 5.4 / 6.4 的
  EXTEND-MODEL 候选是同一份需求；只要 Phase 3 给出 engine-neutral 派生视图，
  CLI 在 Phase 6 迁移时就能直接换 facade。
- **异常 FQN 字符串**（agent）：`agent/_tasks_legacy.py:677` 仍指向
  `utils.formats.cp2k.colvar.ColvarParseError`，是 utils 层的具体异常类。Phase 3
  如果引入 engine-neutral exception（如 `engines.protocols.ConstraintParseError`），
  agent 的契约会跟着升级；现在不动。
- **TI 工作目录生成**（`scripts/TIGen.py:19`）：脚本 generator 解析 SG restart 来
  snap 到最近帧，使用 `parse_colvar_restart` 直读 utils.formats。**这是 generator
  写文件流程，不是分析流程**，与业务层 charge/water/potential 的 engines 化优先级
  不同；Phase 6 迁移时可以选择继续直读或一并经 engines。

### 7.3 Phase 6 迁移注意

- cli / agent / scripts 的入口层迁移**不应**在 Phase 3-5 提前做（避免与业务层迁移
  互相阻塞）；但 Phase 3 设计 facade 时**必须**把这 7 个调用点列为"未来要消费的下
  游"，确保 facade 签名能覆盖它们的使用场景（不只是业务层）。

---

## 横向汇总

### Phase 3 候选总表

| 标签 | 涉及业务 | 简述 | Phase 2 已记录的关键风险 |
|---|---|---|---|
| `NEW-FACADE` | water (1.4) + cli (7.2) | engine-neutral cell 读取（restart/md.inp/未来 POSCAR） | cell 表达形态（`(a,b,c)` vs 3×3 matrix）；非正交 cell 扩展性 |
| `EXTEND-MODEL` | potential (2.4) | 是否引入"轻量帧"覆盖 Fermi-only / distributed 路径 | 与 `parse_md_out_fermi` dict 接口的迁移成本耦合 |
| `EXTEND-MODEL` | slowgrowth (5.4) + constrained_ti (6.4) + cli (7.2) | engine-neutral "约束 MD 运行数据组合视图"（含派生 target/time series） | VASP 的 SHAKE/RATTLE growth 单位语义可能不与 CP2K 无损对齐；派生量所有权 (dataclass property vs facade 算好后回填) |
| `ENGINE-SCOPE` | charge (3.4) + constrained_ti correction (6.4) | Bader 数据归 engines.vasp / engines.charge / 业务直读三选一 | 输入是 VASP 物理产物，业务语义是 charge frame/trajectory，IndexMapper 等工程产物的归属 |

### Phase 3 设计阶段建议的拍板顺序

1. 先决定 `EXTEND-MODEL`（slowgrowth + constrained_ti）的"约束 MD 运行数据组合视
   图"形态——这是 SG / TI 主流程依赖。
2. 同期决定 `NEW-FACADE`（cell）的去留——这是 water 主流程依赖，且 cli 入口层会消
   费。
3. 再决定 `ENGINE-SCOPE`（charge）——这是 Bader / correction 流程依赖，但牵扯
   `engines.vasp` placeholder 的实现节奏。
4. 最后决定 `EXTEND-MODEL`（potential 轻量帧）——风险是与 `engines/CLAUDE.md`
   "dict 接口有意保留"的历史决策对冲，Phase 4 时再做。

### Phase 5/6 候选迁移点（doc-only，不指导实现）

按"调用点 → 期望走 engines facade"为单位拆成两张表，业务层（Phase 5 主迁移目标）
与入口层（Phase 6 迁移对象）口径分开，避免与第 7 节的入口层叙述重复或混淆。
`utils.formats.common.cube` 这种纯算法工具的业务直读不计入迁移点（按 Phase 1 边
界本就该住在 utils 层，见 2.4 STAY-AS-IS）。

**Phase 5 — 业务层候选迁移点**

| 业务节点 | 调用点数 | 迁移后期望 |
|---|---|---|
| `water.WaterAnalysis._common._parse_abc_from_md_inp` | 1 | 走 engines cell facade（NEW-FACADE） |
| `electrochemical.potential.CenterPotential` 内 `parse_md_out_fermi` | 1 | 走 engines fermi facade（与 EXTEND-MODEL "轻量帧" 决议耦合） |
| `electrochemical.charge.Bader.{BaderData, AtomCharges, SurfaceCharge}` 调 `load_bader_atoms` | 3 处 import × 多次调用（共 7 次） | 视 ENGINE-SCOPE 拍板决定 |
| `enhanced_sampling.slowgrowth.SlowGrowth.{from_paths, from_directory}` 内 `ColvarMDInfo` | 2 | 走 engines composite facade（EXTEND-MODEL） |
| `enhanced_sampling.constrained_ti.workflow.standalone_diagnostics` 内 `ColvarMDInfo` | 1 | 同上 |
| `enhanced_sampling.constrained_ti.correction._get_electrode_area` 调 `load_bader_atoms` | 1 | 视 ENGINE-SCOPE 拍板决定；若 charge 留业务层，独立提"从 POSCAR 读 cell"的 facade |

业务层迁移 commit 边界初步设想：按 `EXTEND-MODEL` / `NEW-FACADE` / `ENGINE-SCOPE`
三组各一笔 commit。

**Phase 6 — 入口层候选迁移点（与第 7 节同源，不重复统计）**

| 入口层节点 | 调用点数 | 形式 |
|---|---|---|
| `cli/_params.py` `CellAbcParam.collect` | 1 | `from ..utils.formats.cp2k.cell import ...` |
| `cli/_scripts.py` `_print_sg_cv_info` | 1 | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart")` |
| `cli/_enhanced_sampling.py` `_print_sg_info` | 1 | `lazy_import("md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo")` |
| `agent/_tasks_legacy.py` ExceptionMapping FQN | 1 | 异常 FQN 字符串 `md_analysis.utils.formats.cp2k.colvar.ColvarParseError` |
| `scripts/TIGen.py` `parse_colvar_restart` | 1 | `from ..utils.formats.cp2k.colvar import ...`（generator 写文件流程，非分析） |

入口层迁移**不在 Phase 5 范围**，避免与业务层迁移互相阻塞；Phase 3 设计 facade 时
需要把这 5 个调用点列为下游消费者一并考虑（见 7.2）。

### Phase 2 不回答的问题（明确留给 Phase 3）

- 任何新 dataclass / facade 的类名、字段名、方法名。
- 任何 engine-neutral dataclass 的所有权（住在 `engines.models` vs 更下层）。
- VASP 侧的 metadata / lambda 单位是否能无损映射到当前 CP2K-driven 的 a.u. / fs
  公式。
- Bader 数据归 `engines.vasp` / `engines.charge` / 业务直读的最终选择。
- Fermi `list[dict]` legacy 接口的迁移时机。
- 入口层（cli / agent / scripts）走 engines 的具体 PR 边界（属于 Phase 6 工程
  问题，不是 Phase 3 设计问题）。
