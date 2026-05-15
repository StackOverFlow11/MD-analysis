# `md_analysis.engines` 内部实现准则（当前实现口径）

> 适用范围：`src/md_analysis/engines/`
> 详细开发说明见 `src/md_analysis/engines/CLAUDE.md`。

## 1. 职责边界

- **是**：把每个引擎的具体文件组合转换成 engine-neutral dataclass；提供 `ConstraintMDParser` Protocol + registry；通过模块级 facade 函数（`read_*`）简化常用调用。
- **不是**：业务分析逻辑（积分、收敛诊断、绘图、CSV 写出）。这些归 `electrochemical` / `enhanced_sampling` / `water`。

## 2. 模块结构

| 文件 | 职责 |
|---|---|
| `__init__.py` | 包级 facade：re-export 公共 API + 触发默认 CP2K parser 注册 |
| `protocols.py` | `ConstraintMDParser` Protocol、`ParserInferenceError`、registry 函数、私有 `_REGISTRY` |
| `models.py` | engine-neutral dataclass：`ConstraintMetadata` / `LambdaSeries` / `PotentialFrame` / `FermiRecord`（前两者从 `utils.formats.cp2k_colvar` re-export） |
| `cp2k.py` | `CP2KParser` 类 + 5 个模块级 facade 函数 |
| `vasp.py` | `VASPParser` 占位类，三个 Protocol 方法都 raise NotImplementedError；**不**自动注册 |

## 3. 依赖方向

```
engines/  →  utils/{formats, structure, io}, utils/constants, exceptions
```

**禁止**反向 import `electrochemical` / `water` / `enhanced_sampling` / `scripts` / `cli` / `agent`。每个 phase 后扫描：

```bash
rg "from\s+\.\.+\.?(electrochemical|water|enhanced_sampling|cli|scripts|agent)" \
   src/md_analysis/engines/
# (must be empty)
```

## 4. 关键设计决策

### 4.1 Protocol-first，新引擎零业务改动

新增引擎 = 实现 `ConstraintMDParser` 三方法 + 调 `register_parser("name", ParserCls)`。`enhanced_sampling.constrained_ti.io` 和 `electrochemical.potential` 等上层模块**不动**。

### 4.2 默认注册时机

`engines/__init__.py` 在包级 import 时调用：

```python
register_parser("cp2k", CP2KParser)
```

这样 `infer_parser` / `get_parser("cp2k")` 默认可用。VASP 不在这里注册——placeholder 尚未实现，让 `infer_parser` 在 VASP 标记目录上抛 `ParserInferenceError`，比返回 stub 更诚实。

### 4.3 dataclass 单一真相源

- `ConstraintMetadata` 和 `LambdaSeries` 的 frozen dataclass 定义在 `utils.formats.cp2k_colvar`，`engines.models` 仅 re-export。理由：dataclass 本身是 parser 的返回类型，留在 parser 模块更内聚。
- `PotentialFrame` 和 `FermiRecord` 定义在 `engines.models`（不在 `utils.formats`），因为它们是多个文件解析结果的组合产物（cube + xyz + stdout），engine 层更合适。

### 4.4 VASP placeholder 的依赖反向规避

`utils.formats.vasp_{report,outcar}.py` 的类型注解需要引用 `engines.models` 的 dataclass，但 `utils.formats` 不能在运行时反向 import `engines.models`（会破坏层次方向）。解决方案：

```python
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ...engines.models import ConstraintMetadata, LambdaSeries

def parse_vasp_report_metadata(report_path) -> "ConstraintMetadata":  # string literal
    raise NotImplementedError(...)
```

`TYPE_CHECKING` 块在运行时为 `False`，import 不会执行；字符串注解让静态分析（mypy）仍能解析类型。

### 4.5 `read_fermi_series` 与遗留 dict 形状的桥接

`utils.formats.cp2k_stdout.parse_md_out_fermi` 历史上返回 `list[dict]`，被 `electrochemical.potential.CenterPotential` 用 dict-style 访问（`r["step"]` 等）。Phase 7b1 引入 `FermiRecord` typed model 后，没有强制迁移 caller，理由：

- `parse_md_out_fermi` 仍返回 `list[dict]`（pin 在 `test_parse_md_out_fermi_still_returns_dict`）
- `engines.cp2k.read_fermi_series` 通过 `FermiRecord.from_legacy_dict` 转 typed model
- 新代码应该用 `read_fermi_series`；CenterPotential 的迁移留给后续 entrance refactor

### 4.6 Potential-frame discovery 下沉（Phase 7b2）

`read_continuous_potential_frames` / `read_distributed_potential_frames` 的真实实现在 `engines.cp2k.py`；`electrochemical.potential._frame_source.py` 退化成 thin forwarding wrapper（~100 行，无业务逻辑），保留 `discover_continuous_frames` / `discover_distributed_frames` 旧名以兼容现有 `CenterPotential.py` / `PhiZProfile.py` / `__init__.py` 内部 import。

行为不变量（受 parity test 守护）：
- `frame_start / frame_end / frame_step` 切片语义
- 缺失 `md.out` / `sp.out` / `xyz` 时的 fallthrough（`fermi_raw=None`、跳过 interface 检测等）
- metal-element 自动检测的 fallback 链
- mode B 的 `FileNotFoundError`：root 不存在 / 模式无匹配

## 5. 测试入口

- `test/unit/engines/test_facade.py` — 包级 facade 烟雾测试 + Phase 7b1/7b2 行为测试 + Phase 8 VASP placeholder 测试
- `test/unit/engines/test_parsers.py` — CP2KParser Protocol 兼容 + registry + sniffer（含 fake-engine 注册的清理）
