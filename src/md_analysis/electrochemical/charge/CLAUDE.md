# charge — 开发备忘

## 定位

Bader 电荷分析：计算金属-水界面的表面电荷密度（μC/cm²），指定原子电荷追踪，counterion 自动检测追踪。

## 约定

### 两种表面电荷计算方法
- **counterion**：σ = −Σ(counterion net charge) / area — 电荷中性原理，取所有非水非金属物种的电荷**取反**
- **layer**：σ = Σ(interface layer net charge) / area — 直接求和界面金属层的净电荷

### atoms.info 键
`compute_frame_surface_charge()` **原地修改** `atoms.info`，设置 4 个键：
- `surface_charge_density_e_A2` → `[aligned, opposed]`（e/Å²）
- `surface_charge_density_uC_cm2` → `[aligned, opposed]`（μC/cm²）
- `n_charged_atoms_per_surface` → `[aligned, opposed]`
- `charge_per_surface_e` → `[aligned, opposed]`

### Frame 目录
- 模式：`bader_t*_i*`（`DEFAULT_DIR_PATTERN`）
- 排序：按 `_t(\d+)` 正则提取的数值升序
- 每个目录需含 `POSCAR`、`ACF.dat`、`POTCAR`

### 单位换算
- `E_PER_A2_TO_UC_PER_CM2 = 1602.176634`（精确 CODATA 2018 值）

## 陷阱与历史 Bug

- **Bug 1f76032（关键）**：counterion 方法之前漏了取反 — σ 应该是 counterion 电荷的**负值**（电荷中性原理）
- `n_surface_layers` 参数仅影响 layer 方法，默认=1（仅最外层）
- Counterion 方法中"非水非金属"的判定依赖 `detect_water_molecule_indices()` — 如果水拓扑检测失败，counterion 物种会被误判
- `compute_frame_surface_charge` 要求 atoms 已有 `bader_net_charge` 数组（由 `load_bader_atoms` 设置）
- Charged atom 到界面的分配使用 MIC（minimum image convention）周期性距离，不是简单距离

## 子目录

| 目录 | 用途 |
|---|---|
| `Bader/` | Bader 电荷分析子包 → `Bader/CLAUDE.md` |
