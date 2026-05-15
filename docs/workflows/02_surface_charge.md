# 工作流 02 — 表面电荷 σ + φ 标定

[← 回 workflows 索引](../workflows.md) | [← 回主索引](../README.md)

**目标**：从 VASP Bader 后处理拿到表面电荷密度 σ(t)（双侧），加上 σ→φ 标定外推电极电势。

---

## 1. 物理量

| 量 | 含义 | 单位 |
|---|---|---|
| `σ_aligned`, `σ_opposed` | 法向 +/- 两侧界面的电荷密度（每帧每侧一个数） | μC/cm² |
| `q_atom` | 单原子净电荷（POTCAR ZVAL - Bader 体积积分电荷） | e |
| `area` | 电极表面积（沿表面法向的两 cell 向量叉积模） | Å² |
| `φ` | 外推的电极电势（用预先标好的 σ→φ 映射） | V |

**两种 σ 计算方法**：

| method | 公式 | 适用 |
|---|---|---|
| `counterion` | σ = Σ q(非金属/非水原子) / area | 反离子 + 吸附物模型（电中性单元里只有反离子提供净电荷） |
| `layer` | σ = Σ q(界面金属层原子) / area | 等价但视角不同（金属层吸 e⁻ 平衡反离子） |

物理上等价（电荷守恒），数值上略有差异（取决于电荷归属边界）。

---

## 2. 输入数据结构

每帧一个 VASP single-point 工作目录，集中在一个 root 下：

```
bader_runs/                          ← root
├── bader_t1000_i1000/
│   ├── POSCAR                        # 必需；ASE 读结构
│   ├── ACF.dat                       # 必需；bader 工具产生
│   ├── POTCAR                        # 必需；提供 ZVAL（核电荷）
│   ├── CHGCAR                        # 跑 bader 用，本工具不直接读
│   ├── OUTCAR                        # 不直接读
│   └── (其他 VASP 输出)
├── bader_t1500_i1500/
└── ...
```

目录命名：`bader_t<step>_i<frame>/`，按 `_t(\d+)` 数值排序。

每个目录至少要有 `POSCAR + ACF.dat + POTCAR`，缺任何一个该帧会被跳过 + WARN。

---

## 3. VASP INCAR / 后处理工具要求

### 3.1 VASP INCAR（产生 CHGCAR）

```
LCHARG = .TRUE.        # 必需！否则没有 CHGCAR，bader 跑不了
LAECHG = .TRUE.        # 推荐；输出 AECCAR0 / AECCAR2 用于 chgsum.pl 合成 valence+core
ISMEAR = 0             # 高斯 smearing
SIGMA  = 0.05
PREC   = Accurate
```

参考：`src/md_analysis/scripts/template/INCAR` 是项目内置模板，菜单 411 / 412 直接复制使用。

### 3.2 Bader 工具产生 ACF.dat

VASP 跑完后，用 [Henkelman group bader](http://theory.cm.utexas.edu/henkelman/code/bader/)：

```bash
# 合成 valence + core 电荷密度
chgsum.pl AECCAR0 AECCAR2

# bader 体积积分
bader CHGCAR -ref CHGCAR_sum
# 输出：ACF.dat, BCF.dat, AVF.dat
```

`ACF.dat` 第 5 列（CHARGE）是每个原子的 Bader 体积积分电荷。本工具读这一列 + POTCAR 的 ZVAL 算净电荷：

```
q_atom = ZVAL - ACF.CHARGE
```

### 3.3 POTCAR 元素后缀

POTCAR 第一行 `TITEL = ... <element>_pv` / `_sv` / `_h` 之类后缀会被自动剥离。你不用手改。

---

## 4. 内部数据流示意

### 4.1 σ(t) 时序

```
bader_runs/bader_t<step>_i<frame>/
├── POSCAR ───► ASE.read ──────────► atoms (含 cell + 元素)
├── ACF.dat ──► formats/bader ───────► (atom_idx, bader_charge)
├── POTCAR ───► formats/bader ───────► {element: ZVAL}
                                      │
                                      ▼
                          load_bader_atoms() → atoms 加 .arrays["bader_charge"]
                                      │
                                      ▼
              compute_frame_surface_charge(atoms, method="counterion"/"layer")
                                      │
              ┌───────────────────────┴──────────────────────────┐
              ▼ counterion                                       ▼ layer
         过滤非金属非水原子                            structure.layer 找界面金属层
         按 z > / < midpoint 分两侧                   每层电荷之和
              │                                               │
              └─────────► σ_aligned, σ_opposed (μC/cm²) ◄─────┘

 多帧 → trajectory_surface_charge() → (n_frames, 2) 数组
                                      │
                                      ▼
                       surface_charge.csv + .png
                       列：t_fs, sigma_aligned_uC_cm2, sigma_opposed_uC_cm2
```

### 4.2 σ → φ 标定 + 外推

```
calibration_data.csv                  σ→φ 标定数据
   ├── 列：sigma_uC_cm2, phi_V        （多个不同 σ 跑出来的 (σ, φ) 对）
   │
   ▼
菜单 231/232 → CalibrationFit
   选拟合类型：linear / poly / spline / diff-cap
   │
   ▼
calibration.json
   { "fit_type": "linear", "params": [a, b], "metadata": ... }

 ----------------------------------------------------------------
 σ(t) (从 σ 时序) + calibration.json
   │
   ▼
mapper.predict(sigma_array) → phi_array (V)

   surface_charge.csv 多一列 phi_V
   surface_charge.png 加右轴 φ(t)
   （仅菜单 224 单侧 + φ 模式自动启用此步）
```

---

## 5. 执行步骤

### 5.1 跑 σ 时序（菜单 221 / 222 / 223）

```
$ md-analysis
 Input: 221    # 或 222 / 223

 ---------- Surface Charge (Counterion) ----------

 Bader root dir [.]: ./bader_runs/
 Surface normal axis (a/b/c) [c]: c
 Modify advanced parameters? (y/n) [n]: n

 Discovering bader_t*_i*/ frames... 100 frames found.
 Loading frames: 100%|██████████| 100/100 [00:08<00:00, 12.0it/s]
 Computing surface charge (counterion)...
 Output:
   ./output/electrochemical/charge/counterion/surface_charge.csv
   ./output/electrochemical/charge/counterion/surface_charge.png
```

### 5.2 标定 σ→φ（菜单 231）

需要 CSV 文件，至少 2 个数据点：

```csv
sigma_uC_cm2,phi_V
-3.5,-0.1
-2.0, 0.2
+0.5, 0.6
+2.5, 1.0
```

```
 Input: 231

 ---------- Calibrate from CSV File ----------

 CSV path: ./calibration_data.csv
 sigma column [sigma_uC_cm2]:
 phi column [phi_V]:
 Fit type:
   1) linear   2) polynomial   3) spline   4) differential capacitance
 Choice [1]: 1

 Fit:  φ = 0.182 × σ + 0.451
 R² = 0.998
 Output: ./output/electrochemical/calibration/fit/calibration.json
         ./output/electrochemical/calibration/fit/fit_curve.png
```

### 5.3 σ 时序 + φ 外推（菜单 224）

```
 Input: 224

 ---------- Single-Side Charge + Potential ----------

 Bader root dir: ./bader_runs/
 Side (aligned/opposed) [aligned]: aligned
 Calibration JSON [./output/.../fit/calibration.json]: <Enter>

 surface_charge.csv 加列 phi_V
 surface_charge.png 加右轴 φ(t)
 → ./output/electrochemical/charge/counterion_aligned/
```

---

## 6. 输出位置

```
output/electrochemical/
├── charge/
│   ├── counterion/             # 221
│   ├── layer/                  # 222
│   ├── full/                   # 223（counterion + layer 对比）
│   ├── counterion_aligned/     # 224 单侧 + φ 外推
│   ├── counterion_opposed/
│   ├── tracked/                # 225 追踪指定原子
│   └── counterion_tracking/    # 226 自动检测反离子
└── calibration/
    ├── fit/                    # 231/232
    │   ├── calibration.json
    │   └── fit_curve.png
    └── predict/                # 233
        └── predicted.csv
```

CSV 列约定（共有列）：

```
t_fs                  # 帧时间（fs，如能解析帧 step + dt）
sigma_aligned_uC_cm2  # +法向侧
sigma_opposed_uC_cm2  # -法向侧
sigma_aligned_cumavg  # 累积平均（让收敛肉眼可见）
sigma_opposed_cumavg
phi_aligned_V         # 仅 224 / 已加载 calibration 时
phi_opposed_V
```

---

## 7. 常见坑

- **POTCAR 找不到元素映射**：第一行 `TITEL = PAW_PBE Cu_pv 06Sep2000` — 工具 strip 后缀，认 `Cu`。如果是 `Cu_GW` / `Cu_h` 等罕见后缀也认；如果完全是非标 POTCAR 可能漏。
- **`bader_t*_i*` 命名不对**：缺 `_t<digits>` 段会报错。如果你按别的命名（如 `frame_001`）需要重命名。
- **`ACF.dat` 第 5 列不是 CHARGE**：bader 工具新版默认列序，老版可能不同。手工 `head -2 ACF.dat` 确认。
- **`LCHARG = .FALSE.`**：忘了开就没 CHGCAR，bader 跑不了。重跑 VASP。
- **calibration 数据点太少**：< 3 个点用 linear 还行，poly/spline 容易过拟合。建议 ≥ 5 个 (σ, φ) 数据点。
- **`counterion` vs `layer` 数值差异大**：通常说明界面层定义有问题（金属层聚类容差 921 设的不合适，太松会把次表面当界面）。
- **NPT 数据 cell 变化**：本工具按首帧 cell 计算 area，NPT 中后段 σ 会有偏差。建议 NVT 数据。

---

[← 01_potential_cshe.md](01_potential_cshe.md) | [→ 03_sg_to_ti.md](03_sg_to_ti.md)
