# scripts — 开发备忘

## 定位

自动化工作目录生成：VASP Bader 单点（BaderGen）和 CP2K 约束 MD（TIGen）。不从 `md_analysis.__init__` re-export。

## 约定

### BaderGen
- 生成 VASP 目录：POSCAR + INCAR + KPOINTS + POTCAR（可选）+ script.sh（可选）
- POSCAR 注释行包含 IndexMapper 编码（双射映射元数据）
- 批量目录命名：`bader_t{time}_i{step}`，从 XYZ 注释行 `atoms.info` 提取
- POTCAR 通过 `subprocess` 调用 `vaspkit 103`（需要 vaspkit 在 PATH 中）
- 模板文件通过 `importlib.resources` 访问 `template/` 目录

### TIGen
- 修改 SG 的 inp 文件生成约束 MD 输入：
  - `PROJECT` → `cMD`
  - `TARGET` → snapped CV 值（去掉 `[unit]` 标注，用 bare a.u.）
  - `TARGET_GROWTH` → `0`（去掉 `[unit]` 标注）
  - `STEPS` → 用户值（默认 10000）
  - 确保 `&TOPOLOGY` 有 `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ`
- Frame snapping：目标 CV → 找轨迹中 CV 最近的帧 → 用该帧的实际 CV 作为 TARGET
- 批量两种模式：numeric（直接给 a.u. 值）/ time（linspace 时间范围 → 映射到 CV）
- 不支持自定义单位（避免配位数等复杂 CV 的量纲转换问题）

## 陷阱与历史 Bug

- **inp 文件名不固定**：用户可能 `mv sg.inp sginp1`，所以 TIGen 接受任意路径
- **续算场景**：用户可能已删除 `&TOPOLOGY` 中的 `COORD_FILE_NAME`/`COORD_FILE_FORMAT` — TIGen 会自动补回
- TARGET regex 必须用负向前瞻 `(?!_GROWTH)` 避免误匹配 `TARGET_GROWTH`
- `_modify_inp_for_ti` 中 STEPS 替换仅限 `&MD` 块内（避免误改 `MAX_SCF` 等其他 STEPS）

## 子目录

| 目录 | 用途 |
|---|---|
| `utils/` | IndexMapper（CP2K↔VASP 索引映射）→ `utils/CLAUDE.md` |
| `template/` | VASP 模板文件（INCAR, KPOINTS），通过 importlib.resources 访问 |
