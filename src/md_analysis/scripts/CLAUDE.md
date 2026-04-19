# scripts — 开发备忘

## 定位

自动化工作目录生成：VASP Bader 单点（BaderGen）、CP2K 约束 MD（TIGen）、CP2K 单点电势（PotentialGen）、CP2K 单点 DP 训练（SpGen）。不从 `md_analysis.__init__` re-export。

## 约定

### BaderGen
- 生成 VASP 目录：POSCAR + INCAR + KPOINTS + POTCAR（可选）+ script.sh（可选）
- POSCAR 注释行包含 IndexMapper 编码（双射映射元数据）
- 批量目录命名：`bader_t{time}_i{step}`，从 XYZ 注释行 `atoms.info` 提取
- POTCAR 通过 `subprocess` 调用 `vaspkit 103`（需要 vaspkit 在 PATH 中）
- 模板文件通过 `importlib.resources` 访问 `template/` 目录
- **Agent 任务** `bader_gen_batch`（CLI 412）：新 wrapper
  `generate_bader_batch_with_report` 在 `batch_generate_bader_workdirs` 之上
  加 xyz/script 预检 + cell_abc 长度预检 + 结构化 report（`BaderGenBatchReport`
  含 `workdirs` / `n_frames` / `frame_indices` / `steps` / `times_fs` /
  `generate_potcar`，`to_dict()` 直接 JSON 可序列化）。`element_order` 接受
  list/tuple 两种。**只准备目录，不提交 VASP 作业、不解析 Bader 结果**；
  作业提交归 pbs-auto、结果解析归 `charge_*` 类 API。
- `BaderGenError` 在 agent 契约里映射为 `error_type="validation"`：典型触发
  路径是 `generate_potcar=True` 但 `vaspkit` 不在 `PATH` 上（环境前置条件
  未满足），agent 的修复策略是"修正输入"，与 `validation` 一致；`analysis`
  留给分析本身的数值/逻辑失败。

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

### PotentialGen
- 从 MD 轨迹抽帧生成 CP2K 单点电势计算目录：init.xyz + sp.inp + script.sh（可选）
- sp.inp 由用户提供模板（体系相关，不内嵌到包中），通过 `KEY_SP_INP_TEMPLATE_PATH` 持久化或运行时指定
- CELL ABC 自动替换：从 restart/md.inp 读取 cell 后替换模板中的 `&CELL ABC` 行
- `&TOPOLOGY` 确保 `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ`
- 批量目录命名：`potential_t{time}_i{step}`，与分析模块 `_frame_source.py` 的 `_SP_DIR_RE` 匹配
- 批量模式下 inp 只解析一次（cell 替换 + topology 检查），复用 modified text 写入每个目录
- 共享 inp 修改逻辑由 `_inp_utils.py` 提供（`replace_cell_abc` / `ensure_topology_init_xyz` / `modify_inp_for_sp`）

### SpGen
- **用途区分**：用于 DeePMD 训练数据收集（前端），和 PotentialGen 的 Hartree 电势分析完全独立
- 从 MD 轨迹抽帧生成 CP2K 单点计算目录，输出同样含 init.xyz + sp.inp + script.sh
- **独立配置键** `KEY_DP_SP_INP_TEMPLATE_PATH`：用户可同时持久化"电势分析模板"（含 `V_HARTREE CUBE`）和"DP 训练模板"（含 `PRINT FORCES`），不用来回切换
- **独立目录前缀** `sp_t{time}_i{step}`：与 `potential_t*_i*` 分开，避免电势分析 `discover_distributed_frames()`（默认 `dir_pattern="potential_t*_i*"`）误读
- `_frame_discovery.py` 的 `FRAME_DIR_STEP_TIME_RE` 只匹配 `_t\d+_i\d+` 段，对前缀无约束
- 复用 `_inp_utils.py` 的 cell/topology 修改逻辑，和 PotentialGen 共享
- CLI 菜单：441（单帧）、442（批量），独立 MenuGroup "44 DeePMD SP Preparation"
- Settings 菜单：914 `SetDpSpInpTemplateCmd`
- **Agent 任务**：`sp_gen_batch` 注册到 agent 模块（`_handlers.py`），通过 `batch_generate_sp_workdirs` 直通。PotentialGen 仍是 CLI-only；BaderGen / TIGen 已各自暴露为带完整 contract 的 agent 任务（`bader_gen_batch` / `ti_gen_batch`），通过各自的 `*_with_report` wrapper 直通
- 后端链路：SP 算完 → `dpdata` 或 `cp2kdata` 插件转训练集 → DeePMD-kit 训练

### 共享 helper：`_inp_utils.py`
- 私有模块（前导下划线，不 re-export）
- 纯字符串处理，无 I/O，无副作用
- 函数：`replace_cell_abc(inp_text, a, b, c)`、`ensure_topology_init_xyz(inp_text)`、`modify_inp_for_sp(inp_text, cell_abc)`
- 正则：`_CELL_BLOCK_RE`、`_ABC_LINE_RE`、`_TOPOLOGY_BLOCK_RE`、`_COORD_FILE_NAME_RE`、`_COORD_FILE_FORMAT_RE`
- 由 PotentialGen 和 SpGen 共享使用，避免双份维护

### 共享 helper：`_frame_selector.py`
- 私有模块，统一轨迹帧切片逻辑
- 由 **BaderGen / PotentialGen / SpGen** 三个 batch 函数共享使用；TIGen **不适用**（TIGen 的 time mode 是"linspace → 映射到 CV → nearest 帧"的目标点语义，与区间切片不同）
- 核心抽象：`FrameSelection` frozen dataclass + `iter_selected_frames(xyz_path, selection)` iterator
- `mode` 字段是显式字符串判别器（`"index"` / `"time"`），不是隐式 None 触发 — 目的是让 agent 的 JSON Schema 能清晰暴露 enum 选择
- **Index mode**（默认）：`frame_start/end/step` 零基索引切片，与重构前行为一致
- **Time mode**：要求同时提供 `time_start_fs` / `time_end_fs` / `time_step_fs`；贪心匹配：每个目标 `t_k = t_start + k*t_step`，yield 首个 `time >= t_k` 的帧，然后 `next_target = time + t_step`；区间边界 `[t_start, t_end]` 包含；缺 `atoms.info['time']` 元数据 → 直接 `FrameSelectionError`（早失败）
- **单帧定位**：`resolve_single_frame(xyz_path, *, mode, frame=0, time_fs=None, time_tol_fs=1e-6)` — Single Cmd 专用。Time mode 下做最近邻搜索（假设 time 单调），返回 `(idx, atoms, warnings)` 三元组；`|actual - requested| > tol` 时 warnings 非空，由 CLI/Agent 分别处理
- 错误类型：`FrameSelectionError(MDAnalysisError)` — `dispatch()` 层捕获为 `ERROR_VALIDATION`
- 向后兼容：所有 batch 函数的 `mode` 默认 `"index"`，旧调用无需修改
- CLI 默认 `mode="time"`（通过 `ChoiceParam` 的 `default="time"` 触发）；Python API 默认 `mode="index"`

## 子目录

| 目录 | 用途 |
|---|---|
| `utils/` | IndexMapper（CP2K↔VASP 索引映射）→ `utils/CLAUDE.md` |
| `template/` | VASP 模板文件（INCAR, KPOINTS），通过 importlib.resources 访问 |
