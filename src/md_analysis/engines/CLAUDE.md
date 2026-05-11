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
| `models.py` | engine-neutral frozen dataclass：`ConstraintMetadata`、`LambdaSeries`、`PotentialFrame`、`FermiRecord` |
| `cp2k.py` | `CP2KParser`（Protocol 实现）+ 5 个模块级 facade：`read_constraint_metadata` / `read_lambda_series` / `read_fermi_series` / `read_continuous_potential_frames` / `read_distributed_potential_frames` |
| `vasp.py` | `VASPParser` 占位类，三个 Protocol 方法都 raise NotImplementedError。**不自动注册**到 `_REGISTRY`，避免 `infer_parser` 误派发到 stub |

## 公共 API（`from md_analysis.engines import …`）

| 符号 | 类型 | 说明 |
|---|---|---|
| `ConstraintMDParser` | Protocol | 约束 MD parser 协议（三方法：`is_constraint_directory` / `parse_metadata` / `parse_lambda_series`）|
| `ParserInferenceError` | Exception | `infer_parser` 找不到匹配 parser 时抛出 |
| `CP2KParser` | class | CP2K 实现；可直接构造也可走 registry |
| `ConstraintMetadata` | dataclass | 约束 MD 元数据（target / growth / timestep / fixed atoms / cell）|
| `LambdaSeries` | dataclass | λ(t) 时间序列（shake / rattle / n_steps / n_constraints） |
| `PotentialFrame` | dataclass | 单帧势数据（step / time_fs / cube_path / header / values / fermi_raw / atoms） |
| `FermiRecord` | dataclass | 单条 Fermi 记录（step / time_fs / fermi_raw），含 `from_legacy_dict` 桥接 |
| `register_parser(name, factory)` | fn | 注册新引擎适配器 |
| `get_parser(name)` | fn | 按名查 parser（case-insensitive）|
| `infer_parser(directory)` | fn | 嗅探目录，返回第一个识别它的 parser |
| `resolve_parser(parser \| str)` | fn | 把字符串名或实例都规范化成 parser 实例（`"auto"` 是 caller 责任，本函数会拒绝）|
| `read_constraint_metadata(directory)` | fn | CP2K 约束 MD 目录 → `ConstraintMetadata` |
| `read_lambda_series(directory)` | fn | CP2K 约束 MD 目录 → `LambdaSeries` |
| `read_fermi_series(md_out_path)` | fn | CP2K `md.out` → `list[FermiRecord]` |
| `read_continuous_potential_frames(...)` | fn | 连续 MD cube 序列 → `list[PotentialFrame]`（mode A）|
| `read_distributed_potential_frames(...)` | fn | 分布式 SP 子目录 → `list[PotentialFrame]`（mode B）|

## 关键设计决策

### dataclass 名（Phase 5b rename）

- `ColvarRestart` → `ConstraintMetadata`
- `LagrangeMultLog` → `LambdaSeries`

旧名作为 module-level alias 保留在 `utils/formats/cp2k_colvar.py`，让旧测试不报错；canonical 名只通过 `engines.models` 暴露。

### Registry 注册时机

`engines/__init__.py` 在包级 import 时调用 `register_parser("cp2k", CP2KParser)`。**VASP 不在这里注册**——VASP placeholder 还没有实现，`infer_parser` 在 VASP 标记目录上必须抛 `ParserInferenceError`，而不是返回 stub。

### `parse_md_out_fermi` 返回 dict 不动

`engines.cp2k.read_fermi_series` 通过 `FermiRecord.from_legacy_dict` 把 `parse_md_out_fermi` 的 `list[dict]` 转 typed model；但 `parse_md_out_fermi` 本身仍返回 `list[dict]`。理由：`electrochemical.potential.CenterPotential.py` 还用 dict-style 访问（`r["step"]` 等），强迫迁移会牵涉 caller，Phase 9 范围只到 engines 接口，CenterPotential 迁移留给后续 entrance refactor。

### 与 utils/formats 的方向

VASP 占位文件（`utils/formats/vasp_{report,outcar}.py`）需要 `ConstraintMetadata` / `LambdaSeries` / `FermiRecord` 做类型注解。canonical 路径在 `engines.models`，但 `utils/formats` 运行时不能反向 import `engines/`，所以这些注解通过 `TYPE_CHECKING` 块 + 字符串引用实现。

## Phase 7-8 历史脉络

- Phase 5a：建 `engines/` 骨架，把 `enhanced_sampling/_parsers.py` 的 Protocol + CP2KParser 搬过来；`PotentialFrame` 从 `electrochemical/potential/_frame_source.py` 搬到 `engines/models.py`；VASP placeholder 加好。
- Phase 5b：dataclass rename。
- Phase 6：删 `enhanced_sampling/_parsers.py` shim，业务代码改从 `engines` import。
- Phase 7a：抽 `utils/formats/cp2k_stdout.py` 和 `utils/formats/cp2k_xyz.py`（CP2K stdout / xyz 解析）。
- Phase 7b1：加 `read_constraint_metadata` / `read_lambda_series` / `read_fermi_series` 薄 facade；引入 `FermiRecord`。
- Phase 7b2：下沉 `read_continuous_potential_frames` / `read_distributed_potential_frames` 到 `engines.cp2k`；`electrochemical/potential/_frame_source.py` 退化成 thin forwarding wrapper。
- Phase 8：补 `utils/formats/vasp_*.py` 三个占位文件 + 配套测试，确认 VASP 不进 auto-discovery。
