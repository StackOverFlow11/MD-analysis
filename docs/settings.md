# 配置与默认值

[← 回索引](README.md)

CLI 的默认值（提交脚本路径 / inp 模板 / 分析参数等）持久化在 `~/.config/md_analysis/config.json`。改它有两种方式：菜单 9xx，或者直接编辑 JSON 文件。

---

## Settings 子菜单

`Input: 9` 进入 Settings 顶层。下分四组：

```
 ---------- Settings ----------

 90) Show / Reset
 91) Script Paths
 92) Analysis Defaults
 93) Potential Output

   0) Back / Exit
```

也可以直接输三位代码（如 `911`）跳到具体设置项。

---

## 90 — Show / Reset

| 代码 | 用途 |
|---|---|
| **900** | 显示当前所有配置值（含来源：用户配置 vs 硬编码默认） |
| **909** | 重置所有默认值（清空 `config.json` 里的覆盖） |

```
 Input: 900

 ---------- Show Current Configuration ----------

 Configuration file: /home/shaofl/.config/md_analysis/config.json
 (status: exists)

 [Script Paths]
   VASP submit script:       /home/shaofl/.local/bin/run_vasp.sh    (user)
   CP2K submit script:       /home/shaofl/.local/bin/run_cp2k.sh    (user)
   SP inp template:          /home/shaofl/templates/sp.inp          (user)
   DP SP inp template:       (not set)                              (default)

 [Analysis Defaults]
   Layer clustering tol (Å):  0.5     (default)
   Z bin width (Å):           0.1     (default)
   Theta bin width (deg):     5.0     (default)
   Water O-H cutoff (Å):      1.30    (default)

 [Potential Output]
   Reference scale:           SHE     (default)
   pH:                        0.0     (default)
   Temperature (K):           298.15  (default)
   φ_PZC (V):                 (not set)
```

---

## 91 — Script Paths（提交脚本 / 模板路径）

| 代码 | 用途 |
|---|---|
| **911** | VASP 提交脚本路径（411/412 自动用） |
| **912** | CP2K 提交脚本路径（421/422/431/432/441/442 自动用） |
| **913** | CP2K SP inp 模板（431/432 用） |
| **914** | DeePMD SP inp 模板（441/442 用） |

```
 Input: 911

 ---------- Set VASP Submission Script Path ----------

 Current: (not set)
 New path (or empty to clear): /home/shaofl/.local/bin/run_vasp.sh
   Saved to ~/.config/md_analysis/config.json
```

第一次跑批量生成（412 / 422 / 432 / 442）建议先把对应路径设好，不然每次都要手输。

---

## 92 — Analysis Defaults（数值默认值）

| 代码 | 配置键 | 含义 | 默认值 | 命令影响 |
|---|---|---|---|---|
| **921** | `layer_tol_a` | 金属层周期聚类容差 | 0.5 Å | 1xx 全部、224 |
| **922** | `z_bin_width_a` | Z 轴密度 / 取向 bin 宽 | 0.1 Å | 1xx 全部 |
| **923** | `theta_bin_deg` | θ 角分布 bin 宽 | 5.0° | 104, 105 |
| **924** | `water_oh_cutoff_a` | O-H 键判定截止 | 1.30 Å | 102-105 |

```
 Input: 921

 ---------- Layer Clustering Tolerance (A) ----------

 Current: 0.5  (default)
 New value (or empty to keep current): 0.6
   Saved.
```

> 改了之后下次跑分析自动用新值；不重启 CLI 也生效（每次 `execute` 都会重新读 config）。

---

## 93 — Potential Output（电势输出参考）

| 代码 | 用途 |
|---|---|
| **931** | 设置电势输出参考（211-216 的 U 输出会按这个换算） |

```
 Input: 931

 ---------- Set Potential Output Reference ----------

 Current reference: SHE

 1) SHE  (vs. computational SHE, default)
 2) RHE  (= SHE - 0.0592 × pH)
 3) PZC  (= SHE - φ_PZC)

 Choice [1]: 2
 pH [0.0]: 7.0
 Temperature (K) [298.15]: 298.15
   Saved: reference=RHE, pH=7.0, T=298.15 K
```

PZC 模式还会问 `φ_PZC` 的值（V vs SHE）；预先用 232 / 231 + 233 估出来填进去。

---

## 直接编辑 config.json

`~/.config/md_analysis/config.json` 是 plain JSON，可以直接编辑：

```json
{
  "vasp_script_path": "/home/shaofl/.local/bin/run_vasp.sh",
  "cp2k_script_path": "/home/shaofl/.local/bin/run_cp2k.sh",
  "sp_inp_template_path": "/home/shaofl/templates/sp.inp",
  "layer_tol_a": 0.6,
  "z_bin_width_a": 0.1,
  "theta_bin_deg": 5.0,
  "water_oh_cutoff_a": 1.30,
  "potential_reference": "RHE",
  "potential_ph": 7.0,
  "potential_temperature_k": 298.15
}
```

CLI 启动时如果 JSON 解析失败 → 报 warning + 走默认值（不崩）。

---

## 配置 vs 物理常量

不要混淆：

| `~/.config/md_analysis/config.json` | `src/md_analysis/utils/constants.py` |
|---|---|
| 用户偏好 + 默认值 | 物理常数（`HA_TO_EV`, `AU_TIME_TO_FS`, `BOHR_TO_ANG`, ...） |
| 可改 | **不要改**（CODATA 值，写死） |
| 影响菜单默认提示 | 影响所有数值计算 |

---

[← 回索引](README.md)
