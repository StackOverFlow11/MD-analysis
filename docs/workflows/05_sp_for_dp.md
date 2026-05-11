# 工作流 05 — SP 单点为 DeePMD 训练集抽帧

[← 回 workflows 索引](../workflows.md) | [← 回主索引](../README.md)

**目标**：从 CP2K AIMD 轨迹按帧抽出来，生成 CP2K SP（single-point）工作目录跑高精度能量+力，作为 DeePMD-kit 训练数据。

---

## 1. 跟工作流 04 的关系

| | 04（Bader） | 05（DP SP） |
|---|---|---|
| 抽帧来源 | 同（CP2K AIMD xyz） | 同 |
| 输出 | VASP SP 工作目录 | CP2K SP 工作目录 |
| 模板 | INCAR + KPOINTS（项目内置） | sp.inp（用户提供） |
| 提交脚本 | 菜单 911 设的 VASP | 菜单 912 设的 CP2K |
| 后处理用途 | bader → σ 计算 | 跑完后提取 ETOT + FORCES → DP 训练集 |
| 命名前缀 | `bader_t*_i*/` | `sp_t*_i*/` |

---

## 2. 输入数据结构

```
ti_target_<x>/   或   sg_run/   或任意 CP2K AIMD 目录
├── md.inp                   # cell 解析（auto）
├── md-pos-1.xyz             # ★ 帧来源
└── ...

dp_sp_template/              # 你提供的 CP2K SP inp 模板（菜单 914 配置）
└── dp_sp_template.inp        # 单点 inp，含 &FORCE_EVAL & DFT 设置 + force/energy print
```

⚠️ **不能用 AIMD 的 md.inp 直接当模板**——AIMD 是 `RUN_TYPE MD`，DP 训练数据要 `RUN_TYPE ENERGY_FORCE`（单点能量+力）。

---

## 3. CP2K dp_sp_template.inp 关键设置

```
&GLOBAL
  PROJECT sp                           ! 工具会改成 sp_t<step>_i<frame>
  RUN_TYPE ENERGY_FORCE                ! ★ 必须 single-point
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME ...
    POTENTIAL_FILE_NAME ...
    &MGRID ... &END MGRID
    &QS ... &END QS
    &SCF
      EPS_SCF 1.0E-7                   ! 比 AIMD 严！DP 训练集要高精度
      MAX_SCF 200
      ...
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE &END
    &END XC
    &PRINT
      &FORCES                          ! ★ 必须输出力
        FILENAME force
        ADD_LAST NUMERIC
      &END FORCES
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC <ABC will be replaced>       ! 工具会替换
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME init.xyz         ! ★ 必须，工具会确保这两行存在
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND H ... &END KIND
    &KIND O ... &END KIND
    ! ... 你体系所有元素的 BASIS_SET / POTENTIAL
  &END SUBSYS
&END FORCE_EVAL
```

工具会自动改的项：

| 关键字 | 模板（原） | 输出（改后） |
|---|---|---|
| `PROJECT` | 任意 | `sp_t<step>_i<frame>` |
| `&CELL ABC` | 占位或留空 | 实际 cell（从 md.inp 或 restart 解析） |
| `&TOPOLOGY COORD_FILE_NAME` | 不一定有 | `init.xyz`（确保存在） |
| `&TOPOLOGY COORD_FILE_FORMAT` | 不一定有 | `XYZ`（确保存在） |

---

## 4. 输出数据结构

```
sp_runs/
├── sp_t0_i0/
│   ├── sp.inp               # 模板改造后
│   ├── init.xyz             # 该帧的原子坐标
│   └── run.sh               # CP2K 提交脚本
├── sp_t10_i1/
└── ...
```

跑完后每个目录会多一些 CP2K 输出：

```
sp_t0_i0/
├── sp.inp / init.xyz / run.sh
├── sp.out                   # ★ 总能量 ETOT 在这里
└── sp-force-1_0.xyz         # ★ 每个原子的力（XYZ 格式）
```

---

## 5. 内部数据流

```
md-pos-1.xyz
   │
   ▼ ASE.iread() per-frame iterator + frame_selector
   │
   ▼ for each (step, frame_index):
   │
   atoms = ase.io.read(xyz, index=frame_index)
   atoms.set_cell(cell_abc)
   atoms.set_pbc(True)
   │
   ▼ shutil.copy(template, dir/sp.inp)
   ▼ ase.io.write(dir/init.xyz, atoms, format="xyz")
   │
   ▼ _inp_utils.modify_inp(dir/sp.inp, replacements)
   │     PROJECT → sp_t<step>_i<frame>
   │     ABC     → (a, b, c) Å
   │     COORD_FILE_NAME / COORD_FILE_FORMAT → ensured
   │
   ▼ shutil.copy(config.cp2k_script_path, dir/run.sh)
```

跑完之后（用户手动 / 自动提交），DeePMD-kit 端读取 `sp.out` + `sp-force-1_0.xyz` 作为训练 frame。本工具只负责生成工作目录，**不直接读 sp.out / 不导出 DP-format**——那是 DeePMD-kit 自己的工作。

---

## 6. CLI 步骤（菜单 442）

### 6.1 先把模板路径设好（一次性，菜单 914）

```
$ md-analysis
 Input: 914

 ---------- Set DP SP Inp Template Path ----------

 Current: (not set)
 New path: /home/shaofl/templates/dp_sp_template.inp
   Saved.
```

### 6.2 批量生成（菜单 442）

```
 Input: 442

 ---------- Batch Generate SP Work Directories for DP ----------

 XYZ file: ./md-pos-1.xyz
 Cell parameters [auto]: <Enter>
   Detected ABC = (10.2239, 10.2239, 26.4220) A
 Frame selection mode:
   1) Step range  2) Time range  3) Explicit indices
 Choice [1]: 1
 Frame start: 0
 Frame end: 5000
 Frame stride: 50
 Output dir [./sp_runs/]:
 SP inp template [/home/shaofl/templates/dp_sp_template.inp]: <Enter>
 CP2K submit script [/home/.../run_cp2k.sh]: <Enter>

 Generated 100 directories:
   sp_runs/sp_t0_i0/
   sp_runs/sp_t50_i1/
   ...
   sp_runs/sp_t4950_i99/
 Each contains: sp.inp, init.xyz, run.sh

 Note: existing directories with same _t<step>_i<frame> were skipped (collision check).
```

### 6.3 单帧（菜单 441）

跟 442 一样但只生成一个目录，用于 sanity check 模板是否正确。

---

## 7. 提交集群

```bash
for d in sp_runs/sp_*; do
  (cd $d && qsub run.sh)
done
```

跑完后 DP 训练流程（不在本项目范围）：

```
sp_runs/sp_*/sp.out + sp-force-1_0.xyz
   │
   ▼ DeePMD-kit 自己的脚本（dp prepare / dp train）
   │
DP model
```

---

## 8. 常见坑

- **`SP inp template not configured`**：先去菜单 **914** 设。
- **重名碰撞**：如果 `sp_runs/sp_t<step>_i<frame>/` 已存在，会**跳过**且 warn（不覆盖）。重新生成需要先 `rm -rf sp_runs/`。
- **`RUN_TYPE` 没改成 ENERGY_FORCE**：工具不强制改，你自己模板里要写对。否则跑成 MD 浪费集群时间。
- **`&FORCES` 没开**：DP 训练必须有力，模板里 `&PRINT &FORCES` 段必须存在。
- **`COORD_FILE_NAME init.xyz` 缺失**：工具会自动加上这两行（如果模板里没有）。但元素 `&KIND` 段还是要你模板里写齐。
- **大量 SP 同时提交**：注意集群队列限制 / 磁盘配额。建议每批 100-500 个。
- **DP 训练数据帧间相关性高**：抽帧步长 stride 太小（比如每步抽）训练效率差。建议 stride ≥ 10×τ_corr（τ 来自 SG 或 TI 诊断）。

---

[← 04_bader_batch.md](04_bader_batch.md) | [← 回 workflows 索引](../workflows.md)
