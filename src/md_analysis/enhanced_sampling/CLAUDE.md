# enhanced_sampling — 开发备忘

## 定位

增强采样工作流包。当前包含 slowgrowth 子包。不从 `md_analysis.__init__` re-export — 需直接 `from md_analysis.enhanced_sampling.slowgrowth import ...`。

TIGen 工作目录生成在 `scripts/TIGen.py`（不在此包中）。

## 子目录

| 目录 | 用途 |
|---|---|
| `slowgrowth/` | SG 数据结构、积分、绘图 → `slowgrowth/CLAUDE.md` |
| `constrained_ti/` | 约束 TI 收敛诊断（ACF + F&P block avg + Geweke） → `constrained_ti/CLAUDE.md` |

## 慢增长 (SG) → 热力学积分 (TI) 工作流

### 整体流程

1. **SG 粗扫**：沿反应坐标做 slow-growth，得到粗略自由能曲线，判断反应是否值得深入
2. **TI 精采样**：在反应坐标上选 7–13 个点（TARGET 值），每个点做固定约束 MD（`TARGET_GROWTH = 0`），采集 Lagrange 乘子的系综平均 ⟨λ⟩
3. **数值积分**：对 ⟨λ⟩ vs ξ 积分得到精确 ΔA
4. **迭代细化**：检查 ⟨λ⟩ 是否趋近 0 来定位极值/鞍点/平台区间，在需要的地方补点（SG 看到的初末态未必合理，需要 TI 验证）

### TIGen 模块设计要点

#### TARGET 指定方式（两种）

- **时间模式** `(t_initial_fs, t_final_fs, n_points)`：在时间区间内 linspace → 通过 `compute_target_series` 映射到 CV(a.u.) → snap 到最近的轨迹帧
- **数值模式**：直接给一组 a.u. 值 → snap 到最近的轨迹帧
- **不支持自定义单位**：避免复杂量纲转换（如配位数 CV 量纲难以处理），统一使用 CP2K 默认 a.u.

#### 帧匹配与 TARGET snap

- SG 轨迹帧有固定输出间隔（如每 5 MD steps 一帧）
- 用 `compute_target_series(restart, n_steps)` 重建每个 step 的 CV(a.u.)
- 轨迹帧的 step 信息从 `atoms.info["i"]` 获取
- 若期望 TARGET 落在两帧之间 → snap 到最近帧，实际 TARGET 取该帧对应的 CV 值

#### inp 文件修改规则

- SG inp 文件由用户指定路径（**文件名不一定以 `.inp` 结尾**，用户可能 `mv sg.inp sginp1` 后续算）
- 修改项：
  - `PROJECT` → `cMD`
  - `TARGET [unit] val` → `TARGET <snapped_au>`（**去掉单位标注**，用 a.u.）
  - `TARGET_GROWTH [unit] val` → `TARGET_GROWTH 0`（**去掉单位标注**）
  - `STEPS` → 用户指定值（默认 10000）
- 确保 `&TOPOLOGY` 中存在 `COORD_FILE_NAME init.xyz` 和 `COORD_FILE_FORMAT XYZ`（续算时用户可能已删除这两行）
- 去掉单位的原因：restart 文件中全部使用 CP2K 默认 a.u.，且复杂 CV 的量纲可能难以表达

#### 输出结构

- 每个 TARGET 点生成一个目录：`ti_target_<value>/`
- 目录内容：`init.xyz`（从 SG 轨迹抽取的帧）+ `cMD.inp`（修改后的输入文件）

#### 单帧 vs 批量

- `generate_ti_workdir()`：单个 TARGET 点（用于迭代补点场景）
- `batch_generate_ti_workdirs()`：批量生成（支持时间模式和数值模式）

### TI 收敛诊断模块（constrained_ti/）

已实现，详见 `constrained_ti/CLAUDE.md`。核心流程：
1. ACF → τ_corr, N_eff, SEM_auto
2. Flyvbjerg-Petersen block averaging → SEM_block（pow2 + δSEM 平台检测）
3. Running average drift check
4. Geweke stationarity test
5. sem_final: F&P plateau → ACF fallback
