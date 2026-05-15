# 工作流 04 — Bader 批量生成 VASP 工作目录

[← 回 workflows 索引](../workflows.md) | [← 回主索引](../README.md)

**目标**：从 CP2K AIMD 轨迹按帧抽出来，生成一批 VASP single-point 工作目录跑 Bader 后处理（产生 [02_surface_charge.md](02_surface_charge.md) 需要的输入）。

---

## 1. 输入数据结构

```
ti_target_<x>/   或   sg_run/   或任意 CP2K AIMD 目录
├── md.inp                    # 用于 cell 解析（auto 模式）
├── md-pos-1.xyz              # ★ 帧来源
├── md-1_<step>.restart       # 备用 cell 来源
└── ...

inp_template/                  # 你提供的 VASP 模板（菜单 911/913 配置）
├── run.sh                    # PBS 提交脚本
└── (POTCAR 通过 VASP_PP_PATH 自动拼接)
```

工具内置 INCAR + KPOINTS 模板（在 `src/md_analysis/scripts/template/`）：

```
INCAR    -- 包含 LCHARG=.TRUE., LAECHG=.TRUE., ENCUT, ISPIN 等
KPOINTS  -- Gamma-only （SP 单点）
```

如果你想用自己的 INCAR / KPOINTS，把它们放在你指定的 `inp_template/` 下就行（CLI 会问"模板目录"）。

---

## 2. 输出数据结构

```
bader_runs/                          ← <output_dir>
├── bader_t1000_i100/                # _t<step>_i<frame_index>
│   ├── POSCAR                        # 从 xyz 第 frame=100 帧（step=1000）
│   ├── INCAR                         # 复制自模板
│   ├── KPOINTS                       # 复制自模板
│   ├── POTCAR                        # 按元素拼接（VASP_PP_PATH）
│   └── run.sh                        # 复制自配置（菜单 911）
├── bader_t1500_i150/
└── ...
```

目录命名 `bader_t<step>_i<frame_index>/`：
- `<step>` = MD 绝对步号（从 xyz comment 行的 `i=N` 提取）
- `<frame_index>` = 帧在 xyz 文件中的位置（0-indexed）

---

## 3. VASP INCAR 关键设置（产生 CHGCAR 给 bader）

下面是 [02_surface_charge.md](02_surface_charge.md) 需要的最小设置（项目模板已含）：

```
ISTART = 0
ICHARG = 2
LCHARG = .TRUE.        # ★ 必需！否则没 CHGCAR，bader 跑不了
LAECHG = .TRUE.        # 推荐；产生 AECCAR0 / AECCAR2 用于 chgsum.pl
ISMEAR = 0
SIGMA  = 0.05
ENCUT  = 450
EDIFF  = 1E-06
PREC   = Normal
LREAL  = Auto
ISPIN  = 2
```

完整模板见 `src/md_analysis/scripts/template/INCAR`。

---

## 4. 内部数据流

```
md-pos-1.xyz                  CP2K xyz 多帧文件
   │
   ▼ ASE.iread() per-frame iterator
   │
   ┌─ frame_selector ─┐
   │  Step range / Time range / Explicit indices
   │  → 选出要处理的 (step, frame_index) 列表
   └──────────────────┘
   │
   ▼ for each (step, frame_index):
   │
   ┌─ atoms = ase.io.read(xyz, index=frame_index) ─┐
   │   atoms.set_cell(cell_abc)                     │ cell 来自 md.inp 解析（formats.cp2k_cell）
   │   atoms.set_pbc(True)                          │
   └────────────────────────────────────────────────┘
   │
   ▼
ase.io.write(POSCAR, atoms, format="vasp")
   │
   ▼ 拼接其他文件
shutil.copy(template/INCAR, dir/INCAR)
shutil.copy(template/KPOINTS, dir/KPOINTS)
shutil.copy(config.vasp_script_path, dir/run.sh)
   │
   ▼ POTCAR 拼接
for element in unique_elements(atoms):
    cat $VASP_PP_PATH/<functional>/<element>/POTCAR >> dir/POTCAR
```

---

## 5. CLI 步骤（菜单 412）

### 5.1 单帧（菜单 411）

```
$ md-analysis
 Input: 411

 ---------- Generate Bader Work Directory (single frame) ----------

 XYZ file: ./md-pos-1.xyz
 Frame index: 100
 Cell parameters [auto]: <Enter>
   Detected ABC = (10.2239, 10.2239, 26.4220) A
 Output dir: ./bader_out/
 Generate POTCAR? (y/n) [y]: y
 VASP submit script [/home/.../run_vasp.sh]: <Enter>

 Generated: ./bader_out/bader_t1000_i100/
   POSCAR, INCAR, KPOINTS, POTCAR, run.sh
```

### 5.2 批量（菜单 412）

```
 Input: 412

 ---------- Batch Generate Bader Work Directories ----------

 XYZ file: ./md-pos-1.xyz
 Cell parameters [auto]: <Enter>
 Frame selection mode:
   1) Step range (start, end, stride)
   2) Time range (fs, fs, n_points)
   3) Explicit indices
 Choice [1]: 1
 Frame start: 0
 Frame end: 1000
 Frame stride: 10
 Output dir: ./bader_runs/
 Generate POTCAR? (y/n) [y]: y

 Generated 100 directories:
   bader_runs/bader_t0_i0/
   bader_runs/bader_t10_i1/
   ...
   bader_runs/bader_t990_i99/
 Each contains: POSCAR, INCAR, KPOINTS, POTCAR, run.sh
```

---

## 6. 提交集群

```bash
for d in bader_runs/bader_*; do
  (cd $d && qsub run.sh)
done
```

跑完后用 `chgsum.pl` + `bader` 工具产生 `ACF.dat`，然后跑 [02_surface_charge.md](02_surface_charge.md) 的 σ 分析。

---

## 7. 常见坑

- **POTCAR 拼接失败**：环境变量 `VASP_PP_PATH` 没设，或者路径下找不到对应元素。先 `echo $VASP_PP_PATH && ls $VASP_PP_PATH` 确认。
- **cell auto 解析失败**：md.inp 里 `&CELL ABC` 缺。手输 ABC 即可。
- **生成的 POSCAR 缺原子**：xyz 文件第一行 `<n_atoms>` 不准（被截断），或者 ASE 读取失败。`ase.io.read(xyz, index=N)` 单独验证一下。
- **运行脚本路径没设**：第一次用先去菜单 **911** 设 VASP 提交脚本路径。
- **POTCAR 元素映射混乱**：xyz 里元素符号大小写要标准（`Cu` 不是 `cu`）。
- **NPT 数据**：cell 每步变化，本工具按首帧 cell 计算每帧 POSCAR。NVT 没问题，NPT 建议跑前先 average cell 手输。

---

[← 03_sg_to_ti.md](03_sg_to_ti.md) | [→ 05_sp_for_dp.md](05_sp_for_dp.md)
