# 输入约定 + 排错速查

[← 回索引](README.md)

---

## 输入数据约定

### 1. 帧目录命名规则

需要按帧聚合的命令（21x distributed / 22x charge / 31x correction）依赖固定的目录命名：

| 用途 | 命名 | 排序 |
|---|---|---|
| Bader（VASP 后处理） | `bader_t<step>_i<frame>/` | 按 `_t(\d+)` **数值排序** |
| 分布式 SP 电势 | `potential_t<step>_i<frame>/` | 同上 |
| TI 约束点 | 任意（`ti_target_*` / `pt_001` / 都行） | 内容驱动（看是否含 `*.restart` + `*.LagrangeMultLog`） |
| SG | 单一 traj，不分目录 | 不适用 |

⚠️ 字典序排序会把 `t1000` 排到 `t100` 后面 / `t99` 前面 — 但本工具用正则 `_t(\d+)` 提整数排，**不会**踩这个坑；前提是文件名里有 `_t<digits>_i<digits>` 段。

### 2. Interface 标签

沿表面法向的两个金属-水界面：

- `normal_aligned` — 法向 +轴（如 +c）那一侧的界面
- `normal_opposed` — 法向 -轴（如 -c）那一侧的界面

水分析的 `start_interface` 参数、表面电荷的双侧 σ 列名都用这套。

> 旧版本曾用 `low_c` / `high_c`，已统一替换。

### 3. 单位约定

| 物理量 | 单位 | 在哪用 |
|---|---|---|
| 距离 | Å | 全局，包括输入 cell_abc |
| 能量 | eV | 全局对外；内部 Hartree → eV 自动转换 |
| 时间 | fs | 全局 |
| 分数坐标 | [0, 1) | structure.layer 内部 |
| TI / SG 内部 | a.u.（CP2K 默认） | restart 中所有 ξ / λ / TARGET / TARGET_GROWTH |
| 表面电荷 | μC/cm² | 22x 输出 |
| 电势 | V | 21x / 23x 输出 |

⚠️ **TI 工作目录生成（422）的 TARGET 值用 a.u.**，不是用户单位。`TARGET [angstrom] 1.5` 这种行会被改写成 `TARGET <au_value>`，**单位标注被去掉**。这是为了规避复杂 CV（比如配位数）量纲难表达的问题。

### 4. cell_abc 解析优先级

CLI 问 cell 时默认 "auto"，按以下顺序找：

1. 命令显式传值（最高优先级）
2. 同目录下 `*.restart` 文件的 `&CELL ABC` 块
3. 同目录或父目录下 `md.inp` 文件的 `ABC [angstrom] ...` 行
4. 兜底回退：手输

如果你的轨迹是 NPT 跑出来（cell 变化），auto 模式只取**首帧**的 cell — 不一定准。NPT 数据建议手输平均 cell。

### 5. CP2K 文件最低要求

| 命令组 | 必需文件 |
|---|---|
| 21x continuous | `*.cube` 文件（每帧）+ `md.out` + `md-pos-1.xyz` |
| 21x distributed | `potential_t*_i*/sp.out + sp-v_hartree-1_0.cube` |
| 22x | `bader_t*_i*/{POSCAR, ACF.dat, POTCAR}` |
| 30x SG | `slowgrowth-1*.restart` + `slowgrowth-constraint_force.dat-1.LagrangeMultLog` |
| 31x TI | `ti_target_*/{*.restart, *.LagrangeMultLog}`（命名可改，看内容） |
| 313 修正 | TI 输入 + `ti_target_*/bader/bader_t*_i*/{POSCAR, ACF.dat, POTCAR}` + `calibration.json` |

POTCAR 的元素后缀（`_pv` / `_sv`）会被自动剥离，不用手改。

---

## 常见错误与排查

### "ColvarParseError: No &MD block found"

**症状**：30x / 31x / 422 解析 restart 失败。

**原因**：你给的不是 `*-1_<step>.restart` 文件，而是 `_RESTART.wfn` 之类的 checkpoint。

**修法**：找正确的 `slowgrowth-1_<step>.restart` 或 `cMD-1_<step>.restart`。

### "LagrangeMultLog file is too short / empty"

**原因**：CP2K 作业还没产生足够步，或者作业崩了 LagrangeMultLog 没写。

**修法**：检查 CP2K stdout，看是否 SCF 不收敛 / cell 失败 / 内存超限。

### LagrangeMultLog 中出现 `***`

**原因**：CP2K 数值溢出（约束力过大），通常发生在 SCF 没收敛或起始结构太离谱。

**自动处理**：`***` → `nan`，下游不崩。

**修法**：看 `ti_diagnostics_<i>.png` 的 ACF 面板有没有 NaN gap。多的话重跑这一段（更小步长 / 更稳定起点）。

### "Inconsistent dt across TI points"

**原因**：你混了不同 dt 的 TI 点（比如有的 0.5 fs 有的 1.0 fs）。

**修法**：所有 TI 点必须同 dt。重跑 dt 错的那几个点。

### TI 诊断 `passed=False`

按 `failure_reasons` 字段定位：

| reason | 含义 | 修法 |
|---|---|---|
| `acf_neff_low` | N_eff < 10 地板（自相关时间过长 vs 总帧数，τ_corr 估计失效） | 延长 traj，或加 equilibration |
| `block_plateau_not_reached` | F&P block 平均没找到平台（采样不够） | 同上 |
| `running_avg_drift` | running average 漂移过大 | 加 equilibration（弃前期非平衡段） |
| `geweke_z_high` | 平稳性 z > 阈值（前后段不一致） | 同上 |

可以先试 **`auto_equilibration=True`**：CLI 312 / 313 提问时回答 `y`，自动二分弃前半。

### "No constraint-point directories found"

**原因**：32 / 313 在 root 下找不到任何能识别的目录。

**新行为**（2026-05-10 重构后）：
- 默认 `auto` 模式 + 内容判断 → 子目录里有 `*.restart` + `*.LagrangeMultLog` 才算
- 如果你目录命名是 `pt_001/` / `lambda_3.5/` 这种自定义，**不影响识别**
- 如果是空目录或没这些文件 → 报错

**修法**：检查目录确实包含 CP2K 输出文件；或者目录命名特别（用 `dir_filter` 选）。

### "POTCAR not found / ValueError on element"

**原因**：环境变量 `VASP_PP_PATH` 没设，或者元素映射缺失。

**修法**：

```bash
export VASP_PP_PATH=/path/to/vasp/potcar/database
```

POTCAR 数据库要按 ASE 期望的目录结构（`<VASP_PP_PATH>/<functional>/<element>/POTCAR`）。

### "VASP submit script not configured"

**原因**：菜单 411 / 412 用 911 设的脚本路径，第一次没设。

**修法**：菜单 911 设一次。

### "Cell parameters could not be determined"

**原因**：auto 模式找不到 `*.restart` 也没 `md.inp`。

**修法**：手输 ABC（CLI 会 fallback 提问）。或者把 restart / md.inp 放到轨迹同目录。

### Mac 下 `md-analysis` 命令找不到

**原因**：pip install 没把 console script 装进 PATH。

**修法**：

```bash
# 看 entry point 装到哪：
pip show -f md-analysis | grep bin/
# 或者直接调模块：
python -m md_analysis.cli
```

### 跑了一半中断 / 掉电

**自动续算**：CLI 不做断点续算（每个菜单是一锤子买卖）。

**部分输出已写入**：CSV 和 PNG 会在分析完成时一次性写出，所以中断的话半成品文件不会留下来（不会被误读为成品）。

**集群作业断点**：CP2K 自身的 restart 机制走 `EXT_RESTART`；本工具不参与。集群作业重跑后再过来用 CLI 分析即可。

---

## 怎么报 issue

GitHub 仓库 `StackOverFlow11/MD-analysis`。issue 时附上：

1. 用的菜单代码
2. 完整错误信息（含 traceback）
3. 输入文件结构（`tree -L 2 your_dir/`）
4. CLI 版本（`md-analysis` 启动 banner 上的 Version）

---

[← 回索引](README.md)
