# data_example — 测试样例数据

本目录提供 md_analysis 各模块的最小可复现测试样例。所有数据来自真实研究项目（CuAg₁ SAA / CO2RR），已**瘦身 + 抽稀**以控制体积（总 ~100 MB）。

> 这些数据**只用于测试与文档示例**；不要把分析结果当作物理结论。

---

## 目录索引

| 子目录 | 用途 | 测试覆盖 | 体积 |
|---|---|---|---|
| [`water/`](#water) | 水分析（密度 / 取向 / 吸附层 / 三联图） | CLI 101-105 | 360 KB |
| [`potential/dense/`](#potential-dense) | 连续 MD 模式电势分析（cSHE）— continuous input mode | CLI 211-216（continuous 模式） | 22 MB |
| [`potential/distributed/`](#potential-distributed) | 分布式 SP 电势分析（cSHE）— distributed input mode | CLI 211-216（distributed 模式） | 23 MB |
| [`bader/`](#bader) | Bader 表面电荷分析 | CLI 221-226 | 4 MB |
| [`calibration/`](#calibration) | σ→φ 标定（fit + predict） | CLI 231-233 | 16 KB |
| [`sg/`](#sg) | Slow-Growth 自由能 | CLI 301-302 | 32 MB |
| [`ti/`](#ti) | Constrained TI 收敛诊断 | CLI 311-312 | 9.6 MB |
| [`ti_with_correction/`](#ti_with_correction) | 恒电势修正（TI + bader + calibration） | CLI 313 | 33 MB |

**总计**：~100 MB。

---

## water/

**支持菜单**：CLI 101 / 102 / 103 / 104 / 105。

```
water/
├── md.inp              # CP2K 输入；提供 cell ABC = (10.2239, 10.2239, 26.4220) Å
└── md-pos-1.xyz        # 20 帧水/CuAg 界面轨迹（270 原子，stride=100 抽稀）
```

**来源**：`chem-hpc:/share/home/chem-wangyg/shaofl/projects/03_CuAg_SAA/small/explicit/potential_charge/traj/1k/`

原始 traj 是 1986 帧（34 MB），抽稀为 20 帧（345 KB），等距覆盖 0-9500 fs。

**对应工作流**：[`docs/workflows/`](../docs/workflows/)（水分析独立菜单，无专门 workflow）。

---

## potential/dense/

**支持菜单**：CLI 211 / 212 / 213 / 214 / 215 / 216（**continuous 模式**：所有 cube 文件在同一目录，按帧自然排列）。

```
potential/dense/
├── md.inp                                        # CP2K 输入；cell ABC 在此
├── md.out                                        # CP2K stdout；Fermi 能量在此
├── md-pos-1.xyz                                  # 抽稀轨迹（10 帧 stride=4，对齐 cube 步）
└── md-POTENTIAL-v_hartree-1_<step>.cube          # 10 个均匀分布 cube（step 0/200/400/.../1800）
```

**已抽稀**：原始有 40 个 cube（step 0 → 1950，stride=50），保留 10 个 stride=200。md.inp + md.out 完整保留。

**来源**：从 git commit `8619fe2`（2026-02-13）恢复。

**对应工作流**：[`docs/workflows/01_potential_cshe.md`](../docs/workflows/01_potential_cshe.md)（continuous 段）。

---

## potential/distributed/

**支持菜单**：CLI 211 / 212 / 213 / 214 / 215 / 216（**distributed 模式**：每帧一个 SP 子目录）。

```
potential/distributed/
├── potential_t0_i0/
│   ├── sp.inp                          # CP2K SP 输入
│   ├── sp.out                          # SP stdout（含 Fermi energy 行）
│   ├── init.xyz                        # 单帧结构
│   └── sp_potential-v_hartree-1_0.cube # Hartree 势 cube
├── potential_t500_i500/
├── potential_t1000_i1000/
├── ...                                  # 21 个均匀分布点（t0, t500, ..., t10000）
└── potential_t10000_i10000/
```

**来源**：研究项目 `publish/potential_surface_charge/potential/potential_raw_data/...`（已存在的本地数据）。

**已抽稀**：从原始 201 个 SP 抽 21 个（约 1/10 比例，覆盖完整时间范围）。

**对应工作流**：[`docs/workflows/01_potential_cshe.md`](../docs/workflows/01_potential_cshe.md)。

---

## bader/

**支持菜单**：CLI 221 / 222 / 223 / 224 / 225 / 226，以及 411-412 单帧 / 批量 work-dir 生成的测试基线。

```
bader/
├── single_frame/                       # 单帧样例（CLI 411 或单点分析）
│   ├── POSCAR                          # 结构
│   ├── ACF.dat                         # Bader 体积积分电荷
│   └── POTCAR                          # 提供 ZVAL（核电荷）
└── trajectory/                         # 多帧样例（CLI 22x σ(t) 时序）
    ├── bader_t950_i950/                # 命名约定：bader_t<step>_i<frame>/
    │   └── {POSCAR, ACF.dat, POTCAR}
    ├── bader_t1000_i1000/
    └── bader_t1050_i1050/
```

**已瘦身**：每个目录只保留 `POSCAR + ACF.dat + POTCAR` 三件套，其他 VASP 输出（CHGCAR / WAVECAR / OUTCAR / vasprun.xml / DOSCAR 等 21 个文件）已删除。原始每目录 ~191 MB，现 ~1 MB。

**对应工作流**：[`docs/workflows/02_surface_charge.md`](../docs/workflows/02_surface_charge.md)。

---

## calibration/

**支持菜单**：CLI 231 / 232 / 233。

```
calibration/
├── calibration_data_counterion.csv     # (φ, σ) 数据对，counterion 法
├── calibration_data_layer.csv          # (φ, σ) 数据对，layer 法
└── calibration_counterion.json         # 已 fit 好的 differential-capacitance 标定（counterion）
```

**CSV 格式**：

```csv
potential_V,sigma_uC_cm2
-1.0559,0.0000
-1.4166,-13.5417
...
```

**来源**：`publish/potential_surface_charge/capacitance/differential_capacitance/{counterion,layer}/`。

**对应工作流**：[`docs/workflows/02_surface_charge.md §5.2`](../docs/workflows/02_surface_charge.md)。

---

## sg/

**支持菜单**：CLI 301 / 302，以及 422 批量 TI 生成的输入。

```
sg/
├── angle/              # angle CV
├── constraint/         # 单约束
├── distance/           # 距离 CV
├── distance_combinedCV/ # 组合 CV
├── more_constrain/     # 多约束
└── temp/               # 温度变化
```

每个体系含：

```
sg/<system>/
├── sg.inp                                          # CP2K 输入（&CONSTRAINT &COLLECTIVE）
├── slowgrowth-1.restart                            # 最新 restart（元数据）
├── slowgrowth-1_<step>.restart                     # 中间 restart（部分）
├── slowgrowth-pos-1.xyz                            # 轨迹（生成 TI init.xyz 用）
└── slowgrowth-constraint_force.dat-1.LagrangeMultLog  # ★ λ(t) 时序
```

**对应工作流**：[`docs/workflows/03_sg_to_ti.md`](../docs/workflows/03_sg_to_ti.md)。

---

## ti/

**支持菜单**：CLI 311 / 312，以及 313 的 Phase 1 输入。

```
ti/
└── double_cv/
    ├── 1k/                             # 300 K
    │   ├── ti_target_0.031369/
    │   │   ├── cMD-1.restart                                # 元数据
    │   │   └── cMD-constraint_force.dat-1.LagrangeMultLog   # λ(t)
    │   ├── ti_target_-0.103935/
    │   └── ... 共 9 个约束点
    └── 2k/                             # 200 K
        └── ... 8 个约束点
```

**对应工作流**：[`docs/workflows/03_sg_to_ti.md §7-§9`](../docs/workflows/03_sg_to_ti.md)。

⚠️ **目录命名 `ti_target_<value>` 只是惯例**——重构后（commit `4541369`）TI 流程是内容驱动的（任何含 `*.restart + *.LagrangeMultLog` 的子目录都算 TI 点），ξ 值从 restart 的 `TARGET` 字段读，**不**从目录名解析。

---

## ti_with_correction/

**支持菜单**：CLI 313（恒电势修正）。

```
ti_with_correction/                     # 复刻 ti/double_cv/1k 的 9 个约束点 + 各点加 bader/ 子目录
├── ti_target_0.031369/
│   ├── cMD-1.restart                   # TI 元数据
│   ├── cMD-constraint_force.dat-1.LagrangeMultLog
│   └── bader/                          # ★ 313 修正必需
│       ├── bader_t950_i950/{POSCAR, ACF.dat, POTCAR}
│       ├── bader_t1000_i1000/
│       └── bader_t1050_i1050/
├── ti_target_-0.103935/
└── ... 共 9 个 TI 点
```

每个 TI 点下的 `bader/` 含 3 帧 Bader 数据（复用自 `bader/trajectory/`）。313 需要这些数据计算每点的 σ 系综平均。

**对应工作流**：[`docs/workflows/03_sg_to_ti.md §10`](../docs/workflows/03_sg_to_ti.md)。

需要配套：[`calibration/calibration_counterion.json`](#calibration)（σ→φ 标定）。

---

## 怎么跑通端到端测试

最快的 sanity check：从根目录启动 CLI，然后试每一个菜单：

```bash
$ md-analysis

# 水分析（菜单 105）
 Input: 105
 XYZ trajectory file: data_example/water/md-pos-1.xyz
 ...

# 表面电荷（菜单 221）
 Input: 221
 Bader root dir: data_example/bader/trajectory
 ...

# σ→φ 标定（菜单 231）
 Input: 231
 CSV path: data_example/calibration/calibration_data_counterion.csv
 ...

# Slow-Growth quick plot（菜单 301）
 Input: 301
 Restart file: data_example/sg/distance/slowgrowth-1.restart
 ...

# Constrained TI（菜单 312）
 Input: 312
 TI root directory: data_example/ti/double_cv/1k
 ...

# 恒电势修正（菜单 313）
 Input: 313
 TI root directory: data_example/ti_with_correction
 Calibration JSON: data_example/calibration/calibration_counterion.json
 ...
```

输出会落在 `./output/<module>/<sub>/`。

---

## 数据真实性免责

这些样例是**研究项目的中间数据**（CuAg₁ SAA / CO2 → COOH 体系），具有真实物理意义但**精度不足以发表**（部分 σ_A > 5 meV，部分 TI 点 N_eff < 50）。它们适合：

- ✅ 单元测试 / 集成测试 / 文档示例
- ✅ 学习 md_analysis 使用流程
- ❌ 不适合：作为物理参考、复现发表结果

发表级数据请参考论文附录或联系 fenglinshao02@gmail.com。
