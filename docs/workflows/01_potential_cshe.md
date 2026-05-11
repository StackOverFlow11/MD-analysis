# 工作流 01 — 电极电势 U vs SHE（cSHE）

[← 回 workflows 索引](../workflows.md) | [← 回主索引](../README.md)

**目标**：从 CP2K AIMD 拿到绝对电极电势 U vs SHE（计算氢电极），含 Fermi 能量、中心 Hartree 势、ΔΨ 修正、厚度敏感性。

---

## 1. 物理量与公式

cSHE（computational SHE）公式：

```
U vs SHE = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP
```

| 量 | 含义 | 单位 | 来源 |
|---|---|---|---|
| `E_Fermi` | 每帧 Fermi 能量 | eV | CP2K stdout（`Fermi energy: ...` 行）|
| `φ_center` | slab 中心区域 Hartree 势的空间平均 | eV | 每帧 cube 文件 plane-averaged 后取中段 |
| `ΔΨ_a(H₃O⁺/w)` | 标准氢电极相互作用项 | eV | 物理常数（`utils/constants.py`） |
| `μ(H⁺,g⁰)` | 气相氢离子化学势 | eV | 同上 |
| `ΔE_ZP` | 零点能修正 | eV | 同上 |

---

## 2. 输入数据结构

支持两种输入模式：**continuous**（连续 AIMD，所有 cube 在同一目录）和 **distributed**（每帧一个 SP 子目录）。

### 2.1 continuous（连续 AIMD）

```
md_run/                          ← 你传给 CLI 的 trajectory dir
├── md.inp                        # CP2K 输入（cell ABC 在这里）
├── md.out                        # CP2K stdout（Fermi 能量在这里）
├── md-pos-1.xyz                  # 轨迹文件（提供原子坐标 + cell）
├── md-1_<step>-vh-1_0.cube       # 每帧 Hartree 势 cube（按 stride）
├── md-1_<step>-vh-1_0.cube
├── ...
└── *.restart                     # 可选；cell ABC 的备用来源
```

cube 文件命名约定：CP2K 默认 `<PROJECT>-<run_id>_<step>-vh-1_0.cube`。本工具按文件名里的 `_<step>` 数值排序。

### 2.2 distributed（每帧一个 SP 子目录）

```
sp_runs/                         ← 你传给 CLI 的 root
├── potential_t1000_i1000/
│   ├── sp.inp
│   ├── sp.out                    # 单点 stdout（Fermi 在这）
│   ├── init.xyz
│   └── sp-v_hartree-1_0.cube
├── potential_t1500_i1500/
│   ├── sp.inp / sp.out / init.xyz / sp-v_hartree-1_0.cube
└── ...
```

目录命名：`potential_t<step>_i<frame>/`，按 `_t(\d+)` 数值排序。

---

## 3. CP2K inp 必须打印的内容

### 3.1 V_HARTREE_CUBE（Hartree 势 cube）

```
&FORCE_EVAL
  &DFT
    &PRINT
      &V_HARTREE_CUBE ON       ! continuous 用 SILENT 也行
         STRIDE 8 8 1            ! 节省磁盘；z 方向 1 必须
         APPEND T
         &EACH
           MD <stride>           ! continuous 模式：MD 步频率
           ! GEO_OPT 0 / JUST_ENERGY 1 — 看你 RUN_TYPE
         &END EACH
         ADD_LAST NUMERIC
      &END V_HARTREE_CUBE
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

> ⚠️ STRIDE 的 z（第三个分量）**必须是 1**。本工具沿 z 做 plane average，z 抽稀会破坏 φ(z) 的连续性。
> ⚠️ `APPEND T` 不强求；不开就每帧一个独立 cube 文件，照样能识别。

### 3.2 Fermi 能量

CP2K **默认在 stdout 打印** `Fermi energy: <value> a.u.`（每个 SCF 收敛后），不需要特别开关。但要确保：

- `&GLOBAL PRINT_LEVEL` 不低于 `LOW`（默认就够）
- 没有 `&FORBIDDEN_INTERACTIONS` 之类压制 SCF 输出

如果你的 `md.out` 里 `grep "Fermi energy:"` 找不到行，说明 SCF 没收敛或者 print level 过低，分析会报 `RuntimeError: No frames with Fermi energy data found`。

### 3.3 cell 信息

CLI 自动从 `md.inp` 的 `&CELL ABC [angstrom] X Y Z` 行解析。或者从 `*.restart` 的 `&CELL` 块。**不需要**额外 print 设置。

---

## 4. 内部数据流示意

### 4.1 continuous 模式

```
md_run/
├── md.inp ──────► CellParser ─────► cell_abc (Å)
├── md.out ──────► _parse_md_out_fermi (regex) ─────► [(step, time_fs, E_F_au)]
├── md-pos-1.xyz ► ASE per-frame iter ──────────────► atoms (incl. cell)
└── *.cube ──────► formats/cube ──┐
                                 │
   ┌─────────────────────────────┘
   │
   ▼
PotentialFrame dataclass
   { step, time_fs, atoms, cube_path, fermi_raw, ... }
   │
   ▼
slab_average_potential_ev(cube, atoms, slab_thickness, axis)
   │   找 metal slab → 取 slab 中段做空间平均
   ▼
φ_center(t) (eV) 时序

   φ_center(t) ─────┐
                    ├─► U(t) vs SHE = -E_F + φ_center + 常数项
   E_F(t) ──────────┘
                    │
                    ▼
   electrode_potential.csv + .png

   φ(z) overlay：每帧的 plane-averaged φ(z) 画一起 → phi_z/phi_z_overlay.png
   thickness sweep：扫不同 slab thickness 看 U 收敛 → thickness_sensitivity/...
```

### 4.2 distributed 模式

```
sp_runs/
└── potential_t<step>_i<frame>/
    ├── sp.out ────────────► fermi via FERMI_RE per directory
    ├── init.xyz ───────────► atoms + cell
    └── sp-v_hartree-1_0.cube ─► formats/cube

    ↑ 每个目录组装一个 PotentialFrame，按 _t(\d+) 排序
    ↓
（同 continuous 后续：slab average → cSHE → CSV/PNG）
```

---

## 5. 执行（CLI 菜单 216 = 全套）

```
$ md-analysis
 Input: 216

 ---------- Full Potential Analysis  (includes 211-215) ----------

 Input mode  1) continuous  2) distributed [1]: 1
 Trajectory dir [.]: ./md_run/
 Cell parameters [auto from md.inp]: <Enter>
   Detected ABC = (10.2239, 10.2239, 26.4220) A
 Slab thickness for averaging (A) [3.0]: 3.0
 Modify advanced parameters? (y/n) [n]: n

 Loading frames: 100%|██████████| 156/156 [00:14<00:00, 11.1it/s]
 Computing center slab potential... done.
 Computing Fermi energy... 156 records
 Combining → U vs SHE... done.
 Computing phi(z) overlay... done.
 Sweeping thickness sensitivity (1.0 → 15.0 A)... done.

 Output:
   center:   output/electrochemical/potential/center/center_slab_potential.csv
   fermi:    output/electrochemical/potential/fermi/fermi_energy.csv
   electrode:output/electrochemical/potential/electrode/electrode_potential.csv
   phi(z):   output/electrochemical/potential/phi_z/phi_z_overlay.png
   thicknessoutput/electrochemical/potential/thickness_sensitivity/...
```

### 单步骤跑法

| 你只要... | 用菜单 |
|---|---|
| Hartree 势中心时序 | 211 |
| Fermi 能量时序 | 212 |
| U vs SHE（自动调 211+212 合并） | 213 |
| φ(z) 多帧叠加图 | 214 |
| 厚度敏感性扫描 | 215 |
| 全套 | 216 |

---

## 6. 输出位置 + 文件含义

```
output/electrochemical/potential/
├── center/
│   ├── center_slab_potential.csv  # t_fs, phi_center_ev, phi_z_std_ev
│   └── center_slab_potential.png
├── fermi/
│   ├── fermi_energy.csv           # step, time_fs, fermi_ev
│   └── fermi_energy.png
├── electrode/
│   ├── electrode_potential.csv    # t_fs, U_vs_SHE_V, ...
│   └── electrode_potential.png
├── phi_z/
│   └── phi_z_overlay.png          # 多帧 φ(z) 叠加
└── thickness_sensitivity/
    ├── sensitivity.csv            # thickness_A, mean_U_V, std_phi_z_eV
    └── sensitivity.png            # 双轴：左 U 右 std
```

电势参考默认 SHE。要切 RHE/PZC，看 [settings.md](../settings.md) 的菜单 931。

---

## 7. 常见坑

- **`No frames with Fermi energy data found`**：`md.out` 里 SCF 没打印 Fermi energy 行。检查 SCF 收敛、PRINT_LEVEL。
- **cube 文件找不到**：continuous 模式按 `*-vh-1_0.cube` glob，文件名要符合 CP2K 默认。distributed 模式找 `sp-v_hartree-1_0.cube`。
- **cell 解析失败**：md.inp 里 `&CELL ABC` 缺，会回退手输。NPT 数据建议手输平均 cell。
- **STRIDE z != 1**：cube 沿 z 抽稀，φ(z) 会出现 step 跳变。重跑 cube。
- **U vs SHE 数值离谱**：默认 cSHE 常数（μ(H⁺), ΔE_ZP）来自 utils/constants.py，是文献常用值。如果你用别的水模型 / 不同 functional，可能需要手工换常数。
- **distributed 模式 sp.out 缺 Fermi**：跟 continuous 一样问题，但每个 SP 单独检查。

---

[← 回 workflows 索引](../workflows.md) | [→ 02_surface_charge.md](02_surface_charge.md)
