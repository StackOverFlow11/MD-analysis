# potential — 开发备忘

## 定位

Hartree 势分析工作流：中心势、费米能、电极电位 U vs SHE、φ(z) 剖面、厚度灵敏度扫描。

## 约定

### cSHE 公式
```
U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP
```
常量定义在 `utils/config.py`：`DP_A_H3O_W_EV=15.35`, `MU_HPLUS_G0_EV=15.81`, `DELTA_E_ZP_EV=0.35`

### 输入模式（input_mode）
- `input_mode="continuous"`（默认，模式 A）：单一目录下 cube 文件 + md.out
  - cube 发现：`discover_cube_files(pattern)` glob 匹配
  - Fermi 能：从 `md.out` 正则提取（`Fermi energy:` 行在 `STEP NUMBER` 行之前）
  - 原子坐标：从 xyz 轨迹 stream-parse
- `input_mode="distributed"`（模式 B）：`potential_t{time}_i{step}/` 子目录
  - cube 发现：遍历子目录，每个含一个 cube 文件（默认 `sp_potential-v_hartree-1_0.cube`）
  - Fermi 能：从每个子目录的 `sp.out` 提取（仅取最后一条 `Fermi energy:` 行）
  - 原子坐标：从 cube 文件自身读取（`read_cube_atoms()`），天然包含 cell 信息
  - step/time 从目录名正则提取：`_t(\d+)_i(\d+)`

### 帧数据抽象（_frame_source.py）
- `PotentialFrame` frozen dataclass：统一两种模式的帧数据
  - 字段：`step, time_fs, cube_path, header, values, fermi_raw, atoms`
- `discover_continuous_frames()` — 模式 A 帧发现
- `discover_distributed_frames()` — 模式 B 帧发现
- `_parse_sp_out_fermi()` — 从单点 sp.out 提取 Fermi 能
- 两种模式返回统一的 `list[PotentialFrame]`，下游分析逻辑无差异

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
- `read_cube_atoms()` 已从 PhiZProfile 私有函数提升为 CubeParser 公开函数
