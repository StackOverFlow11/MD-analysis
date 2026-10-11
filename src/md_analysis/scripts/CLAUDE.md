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
- `BaderGenError` 典型触发路径是 `generate_potcar=True` 但 `vaspkit`
  不在 `PATH` 上（环境前置条件未满足）——属"修正输入"类失败，区别于
  分析本身的数值/逻辑失败。

### TIGen
- 修改 SG 的 inp 文件生成约束 MD 输入（parse → 树编辑 → generate，见下"共享 helper"）：
  - `PROJECT_NAME`（规范键 `project_name`）→ `cMD`
  - `TARGET` → snapped CV 值（内部单位裸值；模板里的显式 `[unit]` 标注由 `inp.py` 预处理 W2 换算剥除）
  - `TARGET_GROWTH` → `0`（全部 `&COLLECTIVE` 块，同样内部单位裸值）
  - `STEPS` → 用户值（默认 10000，树路径 `motion.md.steps`）
  - 确保 `&TOPOLOGY` 有 `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ`（无则新建 section）
- **行为变更（refactor/cp2k-inp-parser）**：`colvar_id=None` 语义 = primary（第一个 `&COLLECTIVE` 块），与读侧 `colvars.primary` 对齐；旧 regex 实现写所有块（latent bug，真实语料零回归面）
- Frame snapping：目标 CV → 找轨迹中 CV 最近的帧 → 用该帧的实际 CV 作为 TARGET
- 批量两种模式：numeric（直接给 a.u. 值）/ time（linspace 时间范围 → 映射到 CV）
- 不支持自定义单位（避免配位数等复杂 CV 的量纲转换问题）
- **`generate_ti_batch_with_report`（workflow 后端，Phase 6.6）**：
  - `colvar_id: int | None = None`（None=primary）已 threading 进
    `_plan_ti_targets` + 两次 `_load_trajectory_cv` + per-target
    `generate_ti_workdir`。**纯 threading**：snapping / `ti_target_<cv>`
    dirname / inp 修改规则均不变
  - `overwrite: bool = False`（默认 = pre-write collision guard，撞已存在
    `ti_target_*` 抛 `ValueError`）。`overwrite=True` **仅跳过该 guard**，
    后续仍走 `generate_ti_workdir`（`mkdir(exist_ok=True)` 逐文件覆盖
    `cMD.inp`/`init.xyz`/`script.sh`）—— **绝不 `rmtree`/清目录**，撞目录
    内其他文件保留
  - `verbose: bool = False`：tqdm 进度（与 legacy `batch_generate_ti_workdirs`
    同款）。CLI 422 传 `True`
  - 旧 `batch_generate_ti_workdirs`（legacy，无 collision guard、无 preflight）
    迁移后 CLI 不再调用，但仍是 **public Python API**（`scripts/__init__.py`
    `__all__` + 直接单测）。**弃用计划见 `context4agent/requirements/`
    `short_term.md`**；真正删除须为独立一笔 commit + 用户明确批准 + commit
    message 标 `API BREAK` + 同步 scripts public docs/tests。本阶段不删。

## 陷阱与历史 Bug

- **inp 文件名不固定**：用户可能 `mv sg.inp sginp1`，所以 TIGen 接受任意路径
- **续算场景**：用户可能已删除 `&TOPOLOGY` 中的 `COORD_FILE_NAME`/`COORD_FILE_FORMAT` — TIGen 会自动补回
- TARGET/TARGET_GROWTH 在 cit 树中是两个独立规范键（`target` / `target_growth`），不存在 regex 误匹配问题；写树一律用规范小写键
- `STEPS` 走树路径 `motion.md.steps`，天然只改 `&MD` 块，不会误伤 `MAX_SCF` 等其他同名 keyword
- **树键必须用规范名**（如 `project_name` 而非别名 `PROJECT`）：别名键 generate 不报错但产生重复 keyword，re-parse 抛 `NameRepetitionError`（silent 陷阱）
- **生成物 normalize（用户已接受）**：注释/空行丢失、keyword 大写、跨名同级 section 文本序反转（LIFO，非语义序）；语义序（无参重复 section 的 list 序）由 `inp._stabilize_tree` + 生成后自校验硬性保证；W1 范围展开不可逆（`LIST 7..10` → 全量整数）

### PotentialGen
- 从 MD 轨迹抽帧生成 CP2K 单点电势计算目录：init.xyz + sp.inp + script.sh（可选）
- sp.inp 模板三级解析：显式参数 > `KEY_SP_INP_TEMPLATE_PATH` 持久化配置 > 包内置默认电势模板（`template/sp_potential.inp`，合法 CP2K 输入、无占位符，`&KIND` 块按帧元素用 `ensure_kind_blocks` 生成，默认 `DZVP-MOLOPT-SR-GTH` / `GTH-PBE`）；用户模板按树校验（必须有 `&V_HARTREE_CUBE`、`&SMEAR`，`&KIND` 覆盖全部元素）
- CELL ABC 替换：从 restart/md.inp 读取 cell 后经 `inp.set_cell_abc` 写树（`force_eval.subsys.cell.abc`）；A/B/C 向量形态的模板 warning + no-op（语义等价迁移）
- `&TOPOLOGY` 确保 `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ`
- 批量目录命名：`potential_t{time}_i{step}`，与分析模块 `_frame_source.py` 的 `_SP_DIR_RE` 匹配
- 批量模式下 parse 一次 → 树编辑一次（cell + topology，与帧无关）→ generate 一次，复用生成文本写入每个目录
- 共享 inp 树编辑逻辑由 `_inp_tree_edit.py` 提供（`ensure_topology_init_xyz` / `ensure_kind_blocks` 等），底层 parse/generate 走 `utils.formats.cp2k.inp`

### SpGen
- **用途区分**：用于 DeePMD 训练数据收集（前端），和 PotentialGen 的 Hartree 电势分析完全独立
- 从 MD 轨迹抽帧生成 CP2K 单点计算目录，输出同样含 init.xyz + sp.inp + script.sh
- **独立配置键** `KEY_DP_SP_INP_TEMPLATE_PATH`：用户可同时持久化"电势分析模板"（含 `V_HARTREE CUBE`）和"DP 训练模板"（含 `PRINT FORCES`），不用来回切换
- **独立目录前缀** `sp_t{time}_i{step}`：与 `potential_t*_i*` 分开，避免电势分析 `discover_distributed_frames()`（默认 `dir_pattern="potential_t*_i*"`）误读
- `_frame_discovery.py` 的 `FRAME_DIR_STEP_TIME_RE` 只匹配 `_t\d+_i\d+` 段，对前缀无约束
- 复用 `_inp_tree_edit.py` 的 topology 编辑与 `inp.set_cell_abc`（cell 更新），和 PotentialGen 共享；批量 parse-once / generate-once 文本复用
- CLI 菜单：441（单帧）、442（批量），独立 MenuGroup "44 DeePMD SP Preparation"
- Settings 菜单：914 `SetDpSpInpTemplateCmd`
- **Workflow facade 后端**：SpGen / BaderGen / TIGen 各有 `*_with_report` wrapper，由 `workflows.scripts.run_sp_batch` / `run_bader_batch` / `run_ti_batch` 调用（facade 与 wrapper 参数 1:1，返回 `WorkflowResult`，artifacts 为工作目录路径、extra 携带 report）。PotentialGen 仅 CLI 入口，无 workflow facade。
- 后端链路：SP 算完 → `dpdata` 或 `cp2kdata` 插件转训练集 → DeePMD-kit 训练

### 共享 helper：`_inp_tree_edit.py`
- 私有模块（前导下划线，不 re-export），工作流语义层：知道单点/约束 MD 模板该有哪些 section、生成器拥有哪些 keyword、`&KIND` 默认基组/赝势
- 两层划分：格式原语（`parse_inp` / `generate_inp` / `get_section` / `set_keyword` / `iter_collectives` / `set_cell_abc`）在 `utils/formats/cp2k/inp.py`；本模块只放工作流语义 helper
- 函数（全部原地编辑 cit 树 dict）：`ensure_topology_init_xyz(tree)`、`set_project(tree, name)`、`set_md_steps(tree, steps)`、`set_collective_target(tree, target_au, colvar_id=None)`、`zero_target_growth(tree)`、`ensure_kind_blocks(tree, elements, *, basis_set, potential)`
- **树编辑纪律（双轨）**：TIGen 每目标 TARGET 不同——parse 一次成 base_tree，每目标 `copy.deepcopy` 后编辑再 generate；SpGen / PotentialGen 编辑与帧无关——parse 一次 → 编辑一次 → generate 一次 → 文本复用
- 由 TIGen / PotentialGen / SpGen 共享使用，避免双份维护
- `set_collective_target` 的 `colvar_id=None` = primary（第一个 COLLECTIVE 块）；`zero_target_growth` 清零全部块；二者写内部单位裸值

### 共享 helper：`_frame_selector.py`
- 私有模块，统一轨迹帧切片逻辑
- 由 **BaderGen / PotentialGen / SpGen** 三个 batch 函数共享使用；TIGen **不适用**（TIGen 的 time mode 是"linspace → 映射到 CV → nearest 帧"的目标点语义，与区间切片不同）
- 核心抽象：`FrameSelection` frozen dataclass + `iter_selected_frames(xyz_path, selection)` iterator
- `mode` 字段是显式字符串判别器（`"index"` / `"time"`），不是隐式 None 触发 — 便于在 schema / 结构化输出中清晰暴露 enum 选择
- **Index mode**（默认）：`frame_start/end/step` 零基索引切片，与重构前行为一致
- **Time mode**：要求同时提供 `time_start_fs` / `time_end_fs` / `time_step_fs`；贪心匹配：每个目标 `t_k = t_start + k*t_step`，yield 首个 `time >= t_k` 的帧，然后 `next_target = time + t_step`；区间边界 `[t_start, t_end]` 包含；缺 `atoms.info['time']` 元数据 → 直接 `FrameSelectionError`（早失败）
- **单帧定位**：`resolve_single_frame(xyz_path, *, mode, frame=0, time_fs=None, time_tol_fs=1e-6)` — Single Cmd 专用。Time mode 下做最近邻搜索（假设 time 单调），返回 `(idx, atoms, warnings)` 三元组；`|actual - requested| > tol` 时 warnings 非空，由 CLI 处理
- 错误类型：`FrameSelectionError(MDAnalysisError)` — 属输入校验类失败
- 向后兼容：所有 batch 函数的 `mode` 默认 `"index"`，旧调用无需修改
- CLI 默认 `mode="time"`（通过 `ChoiceParam` 的 `default="time"` 触发）；Python API 默认 `mode="index"`

## 子目录

| 目录 | 用途 |
|---|---|
| `utils/` | IndexMapper（CP2K↔VASP 索引映射）→ `utils/CLAUDE.md` |
| `template/` | VASP 模板文件（INCAR, KPOINTS），通过 importlib.resources 访问 |
