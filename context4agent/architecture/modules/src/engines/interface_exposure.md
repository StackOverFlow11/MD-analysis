# `md_analysis.engines` 接口暴露约定（当前实现）

> 对应代码：`src/md_analysis/engines/`
>
> 本文档定义 `md_analysis.engines` 的符号级公开接口与暴露边界。
> 详细开发说明见 `src/md_analysis/engines/CLAUDE.md`。

## 1. 接口角色定义

`md_analysis.engines` 是 CP2K / VASP 等计算引擎的门面层：

- 输入：每个引擎的原始文件组合（CP2K 的 `*.restart` + `*.LagrangeMultLog`、`md.out` + cube；VASP 的 REPORT/OUTCAR/LOCPOT 等）
- 输出：engine-neutral 的 frozen dataclass（`ConstraintMetadata` / `LambdaSeries` / `PotentialFrame` / `FermiRecord`）
- 上游业务（`electrochemical`、`enhanced_sampling`、`water`）只面向 dataclass 编程，与具体引擎解耦

## 2. 公开符号（`from md_analysis.engines import …`）

`engines/__init__.py` 的 `__all__`：

### 2.1 Protocol + registry

- `ConstraintMDParser` — 约束 MD parser 协议
- `ParserInferenceError` — `infer_parser` 找不到匹配 parser 时抛出
- `CP2KParser` — CP2K 实现类
- `register_parser(name, factory) -> None`
- `get_parser(name) -> ConstraintMDParser`
- `infer_parser(directory) -> ConstraintMDParser`
- `resolve_parser(parser_or_name) -> ConstraintMDParser`

### 2.2 Engine-neutral dataclass

- `ConstraintMetadata` — 约束元数据（含 `project_name` / `step_start` / `time_start_fs` / `timestep_fs` / `total_steps` / `colvars: ColvarInfo` / `lagrange_filename` / `cell_abc_ang` / `fixed_atom_indices`）
- `LambdaSeries` — λ(t) 时间序列（含 `shake` / `rattle` / `n_steps` / `n_constraints` + `collective_shake` / `collective_rattle` properties）
- `PotentialFrame` — 单帧势数据（`step` / `time_fs` / `cube_path` / `header: CubeHeader` / `values: np.ndarray` / `fermi_raw` / `atoms: ase.Atoms | None`）
- `FermiRecord` — 单条 Fermi 记录（`step` / `time_fs` / `fermi_raw`） + `from_legacy_dict(d) -> FermiRecord` 桥接

### 2.3 CP2K 门面（模块级函数）

- `read_constraint_metadata(directory: str | Path) -> ConstraintMetadata`
- `read_lambda_series(directory: str | Path) -> LambdaSeries`
- `read_fermi_series(md_out_path: str | Path) -> list[FermiRecord]`
- `read_continuous_potential_frames(cube_pattern, *, workdir, md_out_path, xyz_path, center_mode, metal_elements, fermi_unit, frame_start, frame_end, frame_step) -> list[PotentialFrame]`
- `read_distributed_potential_frames(root_dir, *, dir_pattern, cube_filename, sp_out_filename, center_mode, metal_elements, layer_tol_ang, frame_start, frame_end, frame_step, verbose) -> list[PotentialFrame]`

### 2.4 历史遗留 dataclass alias（仍可用，但请使用 canonical 名）

`utils.formats.cp2k_colvar` 模块仍 export：

- `ColvarRestart = ConstraintMetadata`（Phase 5b 之前的旧名）
- `LagrangeMultLog = LambdaSeries`

Runtime 是同一个类（`is` 检查为 True）；新代码请用新名。

## 3. 暴露但 NOT 自动注册的引擎

`engines.vasp.VASPParser` 实现 `ConstraintMDParser` Protocol 形状但每个方法 raise `NotImplementedError`，**不**在 `engines/__init__.py` 调 `register_parser`。理由：`infer_parser` 在 VASP 标记目录上必须抛 `ParserInferenceError`，而不是返回 stub。

显式构造合法：

```python
from md_analysis.engines.vasp import VASPParser
p = VASPParser()        # OK, no I/O
p.parse_metadata(path)  # raises NotImplementedError
```

`utils.formats.vasp_{report,outcar,locpot}` 三个占位模块也是同样的"placeholder + 显式 NotImplementedError"模式。

## 4. 推荐导入方式

```python
# 业务工作流（推荐）
from md_analysis.engines import (
    ConstraintMDParser, ParserInferenceError,
    CP2KParser, infer_parser, resolve_parser,
    PotentialFrame, FermiRecord,
    read_constraint_metadata, read_lambda_series,
    read_fermi_series,
    read_continuous_potential_frames,
    read_distributed_potential_frames,
)

# 扩展新引擎适配器
from md_analysis.engines.protocols import register_parser
register_parser("vasp", VASPParser)
```

## 5. 依赖方向

```
engines/  →  utils/{formats, structure, io}, utils/constants, exceptions
```

**禁止**反向 import `electrochemical` / `water` / `enhanced_sampling` / `scripts` / `cli` / `agent`。开发期间用以下命令守护：

```bash
rg "from\s+\.\.+\.?(electrochemical|water|enhanced_sampling|cli|scripts|agent)" \
   src/md_analysis/engines/
# (must be empty)
```

## 6. 稳定性

- Protocol + registry + dataclass + CP2K facade：**Evolving**（utils/engines 重构仍在收尾，字段名稳定但接口可能再调整）
- VASP placeholder：**Unstable**（Phase 8 占位，实际实现时签名可能变动）
