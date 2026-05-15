# engines — 开发备忘

## 定位

CP2K / VASP 等计算引擎的门面层，把每个引擎的具体文件组合（CP2K 的 `*.restart` + `*.LagrangeMultLog`、`md.out` + cube；VASP 的 REPORT/OUTCAR/LOCPOT 等）转换成 engine-neutral 的 frozen dataclass，供上层业务模块（`enhanced_sampling`、`electrochemical`）消费。

新引擎扩展点 = 实现 `ConstraintMDParser` Protocol + `register_parser("name", ParserCls)` + 在 `engines/<name>.py` 加门面函数；上层业务代码不动。

## 依赖方向

```
engines/  →  utils/{formats, structure, io}, utils/constants, exceptions
```

**不允许**反向 import `electrochemical` / `water` / `enhanced_sampling` / `scripts` / `cli` / `agent`。验证：
```
rg "from\s+\.\.+\.?(electrochemical|water|enhanced_sampling|cli|scripts|agent)" src/md_analysis/engines/
# (must be empty)
```

## 模块结构

| 文件 | 职责 |
|---|---|
| `__init__.py` | 包级 facade：re-export Protocol + registry + 公共 dataclass + CP2K 门面函数；触发默认 CP2K parser 注册 |
| `protocols.py` | `ConstraintMDParser` Protocol、`ParserInferenceError`、registry（`register_parser` / `get_parser` / `infer_parser` / `resolve_parser` / 私有 `_REGISTRY`） |
| `models.py` | engine-neutral frozen dataclass：`ConstraintInfo` / `ColvarInfo` / `ConstraintMetadata` / `LambdaSeries` / `ConstraintRun`（约束 MD 族）；`PotentialFrame` / `CenterPotentialScalarFrame` / `FermiRecord`（potential 族）；`CellSpec`（cell 族） |
| `cp2k.py` | `CP2KParser`（Protocol 实现）+ 12 个模块级 facade：约束 MD 目录级 (`read_constraint_metadata` / `read_lambda_series` / `read_constraint_run`) + 约束 MD 文件级 (`read_constraint_metadata_from_restart` / `read_lambda_series_from_log` / `read_constraint_run_from_files`) + 约束 MD 派生 (`compute_target_series`) + Fermi 文件 (`read_fermi_series`) + Cell 文件 (`read_cell`) + Potential 帧 (`read_continuous_potential_frames` / `read_distributed_potential_frames`) + Potential 标量 (`read_center_potential_scalar_frame`) |
| `vasp.py` | `VASPParser` 占位类，三个 Protocol 方法都 raise NotImplementedError。**不自动注册**到 `_REGISTRY`，避免 `infer_parser` 误派发到 stub |

## 公共 API（`from md_analysis.engines import …`）

| 符号 | 类型 | 说明 |
|---|---|---|
| `ConstraintMDParser` | Protocol | 约束 MD parser 协议（三方法：`is_constraint_directory` / `parse_metadata` / `parse_lambda_series`）|
| `ParserInferenceError` | Exception | `infer_parser` 找不到匹配 parser 时抛出 |
| `CP2KParser` | class | CP2K 实现；可直接构造也可走 registry |
| `ConstraintInfo` | dataclass | 单个 CV constraint 参数（colvar_id / target_au / target_growth_au / intermolecular）|
| `ColvarInfo` | dataclass | CV constraints 集合（`primary` / `__len__` / `__getitem__`）|
| `ConstraintMetadata` | dataclass | 约束 MD 元数据（project_name / step_start / time_start_fs / timestep_fs / total_steps / colvars / lagrange_filename / cell_abc_ang / fixed_atom_indices）|
| `LambdaSeries` | dataclass | λ(t) 时间序列（shake / rattle / n_steps / n_constraints；`collective_shake` / `collective_rattle` 派生）|
| `ConstraintRun` | dataclass | 约束 MD 单点 composite（`metadata` + `lambda_series` + 派生 `n_steps` / `steps` / `times_fs` / `target_series_au()`）|
| `CellSpec` | dataclass | 引擎中立 cell 描述（`cell_matrix_ang` (3,3) 行向量 + `pbc`；`abc_ang` / `is_orthorhombic` 派生）|
| `PotentialFrame` | dataclass | 单帧势数据（step / time_fs / cube_path / header / values / fermi_raw / atoms） |
| `CenterPotentialScalarFrame` | dataclass | 标量级 slab-averaged potential（step / time_fs / center_source / center_z_ang / slab_thickness_ang / phi_center_ev / fermi_level_ev / phi_z_std_ev / n_slices） |
| `FermiRecord` | dataclass | 单条 Fermi 记录（step / time_fs / fermi_raw），含 `from_legacy_dict` 桥接 |
| `register_parser(name, factory)` | fn | 注册新引擎适配器 |
| `get_parser(name)` | fn | 按名查 parser（case-insensitive）|
| `infer_parser(directory)` | fn | 嗅探目录，返回第一个识别它的 parser |
| `resolve_parser(parser \| str)` | fn | 把字符串名或实例都规范化成 parser 实例（`"auto"` 是 caller 责任，本函数会拒绝）|
| `read_constraint_metadata(directory)` | fn | CP2K 约束 MD 目录 → `ConstraintMetadata` |
| `read_lambda_series(directory)` | fn | CP2K 约束 MD 目录 → `LambdaSeries` |
| `read_constraint_run(directory)` | fn | CP2K 约束 MD 目录 → `ConstraintRun` |
| `read_constraint_metadata_from_restart(path)` | fn | CP2K `*.restart` 文件 → `ConstraintMetadata` |
| `read_lambda_series_from_log(path)` | fn | CP2K `*.LagrangeMultLog` 文件 → `LambdaSeries` |
| `read_constraint_run_from_files(restart, log)` | fn | `(restart, log)` 文件对 → `ConstraintRun` |
| `compute_target_series(metadata, n_steps, *, colvar_id=None)` | fn | 重建 ξ(t) 目标序列（CP2K SHAKE/RATTLE 公式）|
| `read_fermi_series(md_out_path)` | fn | CP2K `md.out` → `list[FermiRecord]` |
| `read_cell(path)` | fn | CP2K `*.restart` 或 `md.inp` → `CellSpec`（suffix 嗅探）|
| `read_continuous_potential_frames(...)` | fn | 连续 MD cube 序列 → `list[PotentialFrame]`（mode A）|
| `read_distributed_potential_frames(...)` | fn | 分布式 SP 子目录 → `list[PotentialFrame]`（mode B）|
| `read_center_potential_scalar_frame(frame, *, center_z_ang, slab_thickness_ang, ...)` | fn | `PotentialFrame` + slab geometry → `CenterPotentialScalarFrame`（标量级 reduce；显式拒绝 `center_z_ang=None`）|

## 关键设计决策

### dataclass canonical 名

- 约束元数据：`ConstraintMetadata`
- λ(t) 时间序列：`LambdaSeries`
- 单点 composite：`ConstraintRun`

物理位置全部在 `engines.models`；早期重构期一度保留的 CP2K-名兼容 alias 已在命名清理时移除，业务方直接使用 canonical 名。

### Registry 注册时机

`engines/__init__.py` 在包级 import 时调用 `register_parser("cp2k", CP2KParser)`。**VASP 不在这里注册**——VASP placeholder 还没有实现，`infer_parser` 在 VASP 标记目录上必须抛 `ParserInferenceError`，而不是返回 stub。

### `parse_md_out_fermi` 是 utils 解析器输出形态(dict)

`utils.formats.cp2k.stdout.parse_md_out_fermi` 仍返回 `list[dict]` —— 它是底层 parser，dict 是 parser contract 输出形态（仅供 engines facade 内部使用 + 测试 pin parser 形态）。`engines.cp2k.read_fermi_series` 是 typed facade：内部调 parser → `FermiRecord.from_legacy_dict` 把每行转 `list[FermiRecord]`。

业务层（`electrochemical.potential.CenterPotential.fermi_energy_analysis`）已直接消费 `read_fermi_series` 的 typed 输出（业务首次连业务消费完成）。

### 与 utils/formats 的方向

VASP 占位文件（`utils/formats/vasp_{report,outcar}.py`）需要 `ConstraintMetadata` / `LambdaSeries` / `FermiRecord` 做类型注解。canonical 路径在 `engines.models`，但 `utils/formats` 运行时不能反向 import `engines/`，所以这些注解通过 `TYPE_CHECKING` 块 + 字符串引用实现。

## Phase 7-8 历史脉络

- Phase 5a：建 `engines/` 骨架，把 `enhanced_sampling/_parsers.py` 的 Protocol + CP2KParser 搬过来；`PotentialFrame` 从 `electrochemical/potential/_frame_source.py` 搬到 `engines/models.py`；VASP placeholder 加好。
- 重构期 dataclass rename：CP2K 风格 colvar/lagrange 命名统一改为 engine-neutral canonical 名（`ConstraintMetadata` / `LambdaSeries` / `ConstraintRun`，物理位置 `engines.models`）；过渡期保留的 CP2K-名兼容 alias 已在命名清理时移除。
- Phase 6：删 `enhanced_sampling/_parsers.py` shim，业务代码改从 `engines` import。
- formats 抽取期间：抽 `utils/formats/cp2k/stdout.py` 和 `utils/formats/cp2k/xyz.py`（CP2K stdout / xyz 解析；Phase 1 后路径从平铺 `cp2k_*` 迁到 `cp2k/` 子包）。
- Phase 7b1：加 `read_constraint_metadata` / `read_lambda_series` / `read_fermi_series` 薄 facade；引入 `FermiRecord`。
- Phase 7b2：下沉 `read_continuous_potential_frames` / `read_distributed_potential_frames` 到 `engines.cp2k`；`electrochemical/potential/_frame_source.py` 退化成 thin forwarding wrapper。
- Phase 8：补 `utils/formats/vasp_*.py` 三个占位文件 + 配套测试，确认 VASP 不进 auto-discovery。
