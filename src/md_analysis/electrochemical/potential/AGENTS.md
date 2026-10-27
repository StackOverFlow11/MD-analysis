# potential — 开发备忘

## 定位

Hartree 势分析工作流：中心势、费米能、电极电位 U vs SHE、φ(z) 剖面、厚度灵敏度扫描。

## 约定

### 绘图模块
- matplotlib 绑定集中在 `_plot.py`（`plot_series_with_cumavg`/`plot_thickness_sensitivity`/`plot_phi_z_profile`），与分析模块完全解耦
- `CenterPotential.py`、`PhiZProfile.py` 本身不直接 import matplotlib，通过 `_plot.*` 委托出图

### cSHE 公式
```
U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP
```
常量定义在 `utils/constants.py`：`DP_A_H3O_W_EV=15.35`, `MU_HPLUS_G0_EV=15.81`, `DELTA_E_ZP_EV=0.35`

### 输入模式（input_mode）
- `input_mode="continuous"`（默认，模式 A）：单一目录下 cube 文件 + md.out
  - cube 发现：`discover_cube_files(pattern)` glob 匹配
  - Fermi 能：从 `md.out` 正则提取（`Fermi energy:` 行在 `STEP NUMBER` 行之前）
  - 原子坐标：从 xyz 轨迹 stream-parse
- `input_mode="distributed"`（模式 B）：`potential_t{time}_i{step}/` 子目录
  - cube 发现：遍历子目录，每个含一个 cube 文件（默认 `sp_potential-v_hartree-1_0.cube`）
  - Fermi 能：从每个子目录的 `sp.out` 提取（仅取最后一条 `Fermi energy:` 行）
  - 原子坐标：从 cube 文件自身读取（`read_cube_atoms()`），天然包含 cell 信息
  - step/time 从目录名正则提取：`_t(\d+)_i(\d+)`（共享自 `utils/io/_frame_discovery.py`）

### 帧数据抽象（canonical 路径在 `engines/`）

`PotentialFrame` 现在是 engine-neutral dataclass，定义在 `md_analysis.engines.models`。两个发现函数也在 Phase 7b2 下沉到 `engines/cp2k.py`：

| canonical（`engines/cp2k.py`） | legacy wrapper（`electrochemical/potential/_frame_source.py`） | 模式 |
|---|---|---|
| `read_continuous_potential_frames(...)` | `discover_continuous_frames(...)` | A：单目录 cube + md.out |
| `read_distributed_potential_frames(...)` | `discover_distributed_frames(...)` | B：`potential_t*_i*` SP 子目录 |

`_frame_source.py` 自 Phase 7b2 起退化为 **thin forwarding wrapper**（~100 行，无业务逻辑），仅保留 `discover_*_frames` 旧名让 `CenterPotential.py` / `PhiZProfile.py` / `__init__.py` 内部 import 继续工作。新代码请直接从 `md_analysis.engines` import canonical 名。

`PotentialFrame` 字段：`step, time_fs, cube_path, header, values, fermi_raw, atoms`（Phase 4 以来未变）。

Fermi 能解析：连续模式 → `engines.cp2k.read_fermi_series(md_out_path)` 返回 `list[FermiRecord]`；分布式模式 → `utils.formats.cp2k.stdout.parse_sp_out_fermi(sp_out_path)`（单点不附 step 元数据，step 由目录名提供）。

两种模式都返回 `list[PotentialFrame]`，下游分析逻辑（`CenterPotential.py` / `PhiZProfile.py`）无差异。

### 分析模式
- `center_mode="interface"`：需要 xyz 轨迹（连续模式）或 cube 原子坐标（分布式模式）做层检测 → 用界面中点作为 slab 中心
- `center_mode="cell"`：用几何中心（`cell_c / 2`）

### φ(z) Profile 居中
- 每帧的 φ(z) 数组通过 `np.roll` 平移，使 slab 中点对齐到 `cell_c / 2`
- **所有帧使用同一个 shift 值**（从第一帧计算）— 见 bug c397c37

### 厚度灵敏度
- 扫描范围默认 3.5 → 15.0 Å，步长 0.5 Å
- 双轴图：左=mean U vs SHE，右=spatial std φ(z) in slab

## 陷阱与历史 Bug

- **Bug c397c37**：之前每帧各自计算 roll shift → 帧间 φ(z) 不对齐。修复后统一用第一帧的 shift
- **Bug a67a7fe**：slab 居中目标从 slab 分数中心改为 `cell_c / 2`
- Fermi energy 从 `md.out` 正则提取时，`Fermi energy:` 行出现在 `STEP NUMBER` 行**之前**
- Cube 文件的 z-grid 跨帧可能不同 → `phi_z_planeavg_analysis` 将所有帧插值到第一帧的 z 网格
- `discover_cube_files` 的 glob pattern 必须匹配 CP2K 的输出命名（如 `*-HARTREE-*.cube`）
- 分布式模式下，未完成计算的子目录（无 cube 文件）会被自动跳过（debug 级别日志）
- 分布式模式下 sp.out 无 `STEP NUMBER` 行（单点计算），`_parse_sp_out_fermi()` 只提取 Fermi 值，step 从目录名获取
- `read_cube_atoms()` 已从 PhiZProfile 私有函数提升为 `formats.common.cube` 公开函数
