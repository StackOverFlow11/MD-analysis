# `md_analysis.utils` 内部实现准则（当前实现口径）

> 适用范围：`src/md_analysis/utils/`，Phase 2-4 后划分为三个子层：
> - `formats/`：单文件解析器（`cube.py` / `bader.py` / `cp2k_cell.py` / `cp2k_colvar.py` / `cp2k_stdout.py` / `cp2k_xyz.py` / `vasp_{report,outcar,locpot}.py` 占位）
> - `structure/`：几何 / 化学语义（`layer.py` / `water.py` / `cluster.py`）
> - `io/`：路径与 IO 调度（`_frame_discovery.py` / `_io_helpers.py` / `cell_resolver.py`）
> - 顶层：物理常量 `constants.py`。
>
> 目标：在明确物理口径与输入输出契约的前提下，提供可复用、可测试、可维护的底层实现。

> ⚠️ **Phase 9 部分同步状态**：本文档全文做了 sed 路径替换，但下方按旧目录结构组织的章节内容尚未按 formats/structure/io 重写。开发者请以 `src/md_analysis/utils/CLAUDE.md` 为准。

## 1. 职责分层与边界

### `constants.py`

- 只承载默认参数与常量：
  - 元素集合
  - 分箱步长
  - 几何阈值
  - 单位换算相关常量（`HA_TO_EV`、`BOHR_TO_ANG`、cSHE 常量）
  - 轴映射与面积向量索引（`AXIS_MAP`、`AREA_VECTOR_INDICES`）
  - 界面标签常量（`INTERFACE_NORMAL_ALIGNED`、`INTERFACE_NORMAL_OPPOSED`）
  - 电荷方法常量（`CHARGE_METHOD_COUNTERION`、`CHARGE_METHOD_LAYER`）
- 禁止写任何业务计算逻辑。

### `__init__.py` — 子模块直接导入约定

- `utils/__init__.py` 仅包含包定位说明和 `__all__ = []`，不 re-export 任何符号。
- 原则：任何调用方（包内/测试/外部用户）必须按 **子模块直接路径** 导入，例如：
  - `from md_analysis.utils.constants import HA_TO_EV`
  - `from md_analysis.utils.structure.layer import detect_interface_layers`
- 理由：集中 re-export hub 要为每次子模块变更同步维护，且相对子模块直接路径无额外信息。子模块直接路径使调用站点自文档化，并消除 `__init__.py` 作为单点故障源。
- 子模块内部仍应声明 `__all__`（见各 `interface_exposure.md`）以表达公开意图。

### `_io_helpers.py`

- 私有模块（`_` 前缀），不属于公开 API，不出现在 `__all__` 中。
- 提供全模块共享的 I/O 与数值 helper：
  - `_cumulative_average(values)` — 1-D 数组逐元素累积平均
  - `_write_csv(path, rows, fieldnames)` — 将字典行列表写入标准 CSV 文件（自动创建父目录）
  - `_write_csv_from_arrays(path, columns)` — 将命名 numpy 数组按列写入标准 CSV 文件。`columns` 为 `dict[str, np.ndarray]`，所有数组须等长
- 两个 CSV 写入函数输出格式完全一致（标准 CSV：首行为逗号分隔列名，无 `#` 前缀），使用方按数据形态选择：dict rows 用 `_write_csv`，numpy arrays 用 `_write_csv_from_arrays`
- 不依赖 ASE 或其他上层模块。

### `_frame_discovery.py`

- 私有模块（`_` 前缀），不属于公开 API。
- 统一了 Bader (`bader_t{time}_i{step}`) 与 SP Potential (`potential_t{time}_i{step}`) 帧目录的发现与排序逻辑。
- 公开符号：
  - `FRAME_DIR_STEP_TIME_RE` — 编译正则 `_t(\d+)_i(\d+)`
  - `extract_step_time_from_dirname(dirname)` → `(step, time_fs)` 或 `None`
  - `discover_frame_dirs(root, glob_pattern, *, required=True)` → 按 `(time_fs, step)` 升序排列的 `list[Path]`
- 两种 frame 目录约定下，`time_fs` 和 `step` 对真实轨迹均单调，元组排序保证稳定。
- 被 `electrochemical/charge/Bader/_frame_utils.py`（Bader 目录）和 `electrochemical/potential/_frame_source.py`（SP 目录）共同调用，取代原本三处各自的正则与排序实现。

### `structure/cluster.py`

- 负责 1D 周期性聚类、最大间隙检测、间隙中点计算。
- 被 `structure/layer.py` 和 `md_analysis.potential.CenterPotential` 调用。
- 不依赖 ASE 或其他上层模块。

### `cube.py`

- 负责 Gaussian cube 文件的读取与解析。
- 负责 plane-averaged φ(z) 计算和 slab-averaged potential 计算。
- 负责 cube 文件发现与切片（`discover_cube_files`）。
- 输入：cube 文件路径或已解析的 header + values。
- 输出：`CubeHeader` 数据结构、φ(z) 数组、`(phi_center_ev, info)` 等。
- CP2K cube 文件约定：z 为 fastest-running index，reshape 为 `(nx, ny, nz)`。
- 单位：cube 内部为 Bohr/Hartree，输出转为 Angstrom/eV。
- `discover_cube_files(cube_pattern, *, workdir, frame_start, frame_end, frame_step)` — 在 `workdir` 下按 glob 匹配 cube 文件，词法排序后按 `[frame_start:frame_end:frame_step]` 切片返回。无匹配时抛出 `FileNotFoundError`。被 `potential/CenterPotential.py` 调用，取代原来内联的 glob+检查逻辑。

### `structure/layer.py`

- 负责金属层识别、界面层标记、法向符号判定。
- 负责层级数据结构与摘要输出。
- 不负责水分子拓扑识别、角度 PDF 统计。
- 使用 `constants.py` 中的 `AXIS_MAP` 常量（取代原有的模块局部 `_AXIS_MAP` 字典）。
- 使用 `constants.py` 中的 `INTERFACE_NORMAL_ALIGNED`/`INTERFACE_NORMAL_OPPOSED` 常量作为界面标签。
- `circular_mean_fractional()` 委托给 `cluster._circular_mean(values, period=1.0)`。

### `formats/cp2k/inp.py`

- CP2K `.inp` 文本 parse/generate 封装层，基于 cp2k-input-tools（cit，git-pinned 依赖）+ pint；**全包唯一直接 import cit/pint 的模块**，两者均惰性单例（`lru_cache`，首次调用才付 import 开销）。
- `parse_inp(text_or_path, *, xmlspec=None) -> dict`：预处理管线 **CRLF 归一 → G1（`@` 预处理指令门禁，拒绝而非改写）→ W3（LIST 行逗号归一）→ W1（整数范围 `N..M` 展开，含 `\` 续行）→ W2（TARGET/TARGET_GROWTH 显式单位换算为内部单位并剥 token）→ 引号感知注释清洗**，再交给 cit simplified parser。
- `generate_inp(tree, *, xmlspec=None) -> str`：`_stabilize_tree` 预反转 list-of-dicts（无参重复 section）抵消 generator 的 LIFO 反转；生成后 re-parse 与调用方原树全等自校验（键序不敏感、list 序敏感），不一致抛 `InpGenerateError`。
- 树约定：键 = 规范名小写（别名 parse 时归一，如 `PROJECT` → `project_name`；读写一律用规范名）；`_` 键 = section 参数；`*` 键 = 默认 keyword 行；keyword 值落在 XML 默认单位（约束类 TARGET/TARGET_GROWTH 为 CP2K 内部单位，与 `.restart` 写出值同口径）。
- 格式原语 helper：`get_section` / `get_keyword` / `set_keyword` / `iter_collectives` / `set_cell_abc`（A/B/C 向量形态模板 warning + no-op）；工作流语义 helper 不在此层（见 `scripts/_inp_tree_edit.py`）。
- 异常：`InpParseError` / `InpGenerateError`（均继承 `MDAnalysisError`）。

### `formats/cp2k/cell.py`

- 统一的 CP2K cell 参数解析模块，负责从 `.restart` 和 `md.inp` 文件解析正交 cell 参数；两个入口都先经 `formats.cp2k.inp.parse_inp` 解析成 cit 树再读树。
- `parse_abc_from_restart`：读 `force_eval.subsys.cell.{a,b,c}` 向量（angstrom 语义）的对角元素，经共享 helper `_abc_ang_from_cell_tree` 验证正交性（离对角元素 < 1e-6 Å）。
- `parse_abc_from_md_inp`：读 `force_eval.subsys.cell.abc` 三标量；cit 自动把任意显式单位 token 换算到默认单位（angstrom），`ABC [bohr]` / 裸值形态同样正确（旧 regex 只匹配 `ABC [angstrom]`，此为行为超集）。ABC 三标量原理上无法检测非正交，正交校验仅 restart/向量路径有。
- 两个函数返回类型一致：`(float, float, float)`，可互换使用；"parse 成功但树缺块"场景错误文案与旧实现逐字一致，parse 本身失败则抛 `inp.InpParseError`（带 cit 行号）。
- `_abc_ang_from_cell_tree(tree, src_path)` 由 restart 路径与 colvar 路径共享（同一棵树），保证正交政策处处一致。
- 统一异常类型 `CellParseError`。
- 不依赖 ASE。

### `formats/cp2k/colvar.py`

- 负责 CP2K COLVAR（集合变量约束）模拟的 restart 文件和 LagrangeMultLog 文件解析；restart 经 `formats.cp2k.inp.parse_inp` 解析成 cit 树后按树路径直读（`global.project_name`、`motion.md.*`、`motion.constraint.collective.*` 内部单位直读、`motion.constraint.lagrange_multipliers.filename`、`motion.constraint.fixed_atoms.list`）。
- cell 参数复用 `formats.cp2k.cell._abc_ang_from_cell_tree()`（同一棵树），避免重复实现，正交政策与 cell 路径一致。
- 多个 `&COLLECTIVE` 块在树中为 list 形态，天然支持多约束场景，结果包装为 `ColvarInfo` 容器；校验次序钉死：PROJECT_NAME → MD → CONSTRAINT → COLLECTIVE → cell（lagrange/fixed_atoms 无块不抛）。
- 不解析 COLVAR 定义（CV 类型由 CP2K 决定，种类繁多）。
- `&FIXED_ATOMS` 的 `LIST` 行：逗号归一（W3）、`N..M` 范围展开（W1）、`\` 续行由 inp.py 预处理完成，本模块仅显式 `sorted` 保证升序 tuple 契约。
- LagrangeMultLog 格式自动检测：单约束（每行一个值）vs 多约束（多行矩阵块）；纯文本逐行解析，与 cit 无关。
- 不完整末尾处理：末尾只有 Shake 无配对 Rattle 时静默丢弃。
- 不依赖 ASE。

### `structure/water.py`

- 负责水分子标记、氧索引提取、质量密度与取向统计。
- 负责 c 轴窗口角度 PDF 统计。
- 不负责界面层判别与金属层聚类。

## 2. 坐标与几何口径（必须一致）

### 2.1 z/c 方向统计统一口径

- z 分箱基于分数坐标第三轴与 `|c|`：
  - `z = f_c * |c|`
- 取向参考方向使用晶胞 `c_unit = c / |c|`，不使用笛卡尔固定 z 轴替代。

### 2.2 周期性与区间规则

- 分数坐标统一使用 wrap 后值（`[0, 1)`）。
- c 轴窗口区间采用左闭右开：
  - `start < end`：`[start, end)`
  - `start > end`：跨边界窗口 `[start, 1) U [0, end)`
  - `start == end (mod 1)`：全区间

### 2.3 分箱边界规则

- z 分箱数：`nbins = ceil(Lz / dz)`
- 最后一段边界直接落在 `Lz`，避免累计误差
- bin 归属使用 `searchsorted(..., side="right") - 1` 再 clip，确保边界可控

## 3. `layer.py` 实现准则

### 3.1 输入与默认值

- `metal_symbols=None` 时，必须回退到 `DEFAULT_METAL_SYMBOLS`
- `normal` 仅支持 `"a" / "b" / "c"` 轴标签（自定义向量不支持，ValueError）
- 界面层默认每侧 1 层（共 2 层）；下游模块（`charge/BaderAnalysis.py`）通过 `n_surface_layers` 参数支持每侧多层

### 3.2 分层与界面判定

- 先按法向投影做 1D 聚类
- 界面层仅保留直接面向非金属环境的一层（normal_aligned/normal_opposed 各一层）
- 对 `"a" / "b" / "c"` 轴采用分数坐标 MIC 判法向；自定义向量法向不支持（ValueError）
- `Layer.center_frac`：层中心的分数坐标（[0,1)），由 `circular_mean_fractional` 计算
- `Layer.interface_label`：`"normal_aligned"` | `"normal_opposed"` | `None`
- `metal_layers_sorted` 排列：`[normal_aligned, slab_interior…, normal_opposed]`（穿过 slab 遍历）
- `SurfaceDetectionResult` 提供 `.interface_normal_aligned()` 和 `.interface_normal_opposed()` 访问器

### 3.3 数据结构约束

- `Layer` 与 `SurfaceDetectionResult` 保持 dataclass 不可变语义
- 对外返回应优先使用 tuple，避免可变容器泄漏

## 4. `water.py` 实现准则

### 4.1 水分子标记

- O-H 连通判据基于 MIC 距离与 `oh_cutoff_A`
- 分配策略为"按距离排序的贪心匹配"
- 每个 H 在最终分配中最多归属一个 O

### 4.2 质量密度 z 分布

- 一个 O 计作一个 H2O 分子
- 质量密度单位固定为 `g/cm^3`
- 体积来自 `A_xy * Δz`，其中 `A_xy = ||a × b||`

### 4.3 取向加权 z 分布

- 角度向量为 H-O-H 角平分线（O 为起点）
- `cos(theta)` 必须按 `dot(bisector, c_unit) / |bisector|`
- 结果单位固定为 `g/cm^3`（公式：$\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin,cm^3}}$）

### 4.4 c 窗口角度 PDF

- 输入窗口为分数坐标 c 区间
- 输出长度为 `180 / ndeg`
- `ndeg` 允许浮点；要求 `180 / ndeg` 为整数（容差校验）
- PDF 单位固定为 `degree^-1`
- 若窗口内无样本，返回全零数组（不抛错）

## 5. 参数管理准则

- 默认参数只能从 `constants.py` 读取，不得在函数内部重复硬编码
- 新增默认参数时必须：
  1. 写入 `constants.py`
  2. 在调用链中显式接入
  3. 更新 `data_contract.md` 对应条目

## 6. 异常与错误信息准则

### 6.1 异常分类

- 所有领域异常均继承 `MDAnalysisError`（定义在 `md_analysis.exceptions`）
- 几何/界面类异常优先使用 `SurfaceGeometryError`
- 水拓扑/取向类异常优先使用 `WaterTopologyError`
- 参数合法性问题使用 `ValueError`

### 6.2 错误信息要求

- 报错信息必须包含最小定位信息（参数名、索引、符号名或条件）
- 禁止吞错；捕获异常时需保留原问题语义

## 7. 数值稳定性准则

- 所有归一化前必须检查向量范数是否为 0
- `cos(theta)` 在反三角函数前需 clip 到 `[-1, 1]`
- 浮点整除判定（如 `180/ndeg`）使用容差，不使用纯整数 `%`

## 8. 性能与可扩展性准则

- 优先 vectorized 处理批量数据，循环仅用于必要几何步骤
- 大轨迹统计应采用逐帧累积策略，避免一次性持有全部帧
- 输出结构优先使用 numpy 数组，保持与后续统计/绘图链路兼容

## 9. 文档与测试同步准则

发生以下改动时，必须同步测试与文档：

- 新增/修改公开函数
- 默认参数变化（如 `dz`、`ndeg`、阈值）
- 统计口径变化（单位、方向、分箱规则、窗口规则）

同步范围至少包括：

- `test/unit/utils/`
- `test/integration/utils/`
- `context4agent/architecture/modules/data_contract.md`
- `context4agent/architecture/modules/glossary_units.md`（若涉及术语/单位）

## 10. 提交前检查清单

- [ ] 模块职责是否越界（facade/实现/config 是否混杂）
- [ ] c 方向口径是否全链路一致
- [ ] 默认参数是否全部来自 `constants.py`
- [ ] 异常类型与报错信息是否清晰
- [ ] 单元与集成测试是否通过
- [ ] 契约文档是否同步更新
