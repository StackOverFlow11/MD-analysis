# 主入口重构计划：`main.py` facade + `workflows/` 分层

> 当前取向：本轮不维护 legacy agent 层兼容性。`main.py` 应设计为用户脚本、notebook、CLI 都能依赖的**程序化主入口**，但不承载具体业务实现。

---

## 问题判断

当前 `src/md_analysis/main.py` 同时承担了三类职责：

1. 对外暴露 `run_water_analysis`、`run_potential_analysis`、`run_charge_analysis` 等编程入口。
2. 保存 agent-facing report dataclass 和 `*_with_report` wrapper。
3. 实现部分组合 workflow，例如 `run_all`。

这导致两个问题：

- **冗长**：入口文件超过千行，里面混有输出路径、artifact 校验、agent report、具体 workflow 调度。
- **不完整**：CLI 已覆盖 water / potential / charge / calibration / enhanced_sampling / scripts，但 `main.py` 主要只覆盖 water / potential / charge；`run_all` 实际也只跑 water + potential。

本轮重构的目标不是继续往 `main.py` 填函数，而是把它变成一个薄入口。

---

## 与底层重构计划的关系

本计划应排在 `current_reconstructions.md` 之后执行。理由是：入口层依赖底层模块路径和 engine 边界；如果先改 `main.py` / CLI，再做 utils + engines 重构，会重复修改同一批 import 和 `lazy_import` 字符串。

| 冲撞点 | `current_reconstructions.md` | 本计划 | 建议处理 |
|---|---|---|---|
| `electrochemical/potential/_frame_source.py` | Phase 7 把散落 CP2K 解析收敛到 `engines/cp2k.py` | Phase 2 搬迁 `run_potential_analysis` 到 `workflows/potential.py` | 可顺序执行；文件不同但语义相关，先固定底层 |
| CLI `lazy_import` 字符串 | Phase 2-4 更新旧 `utils.*Parser` 路径 | Phase 6 让 CLI 调用 `workflows` / `main` facade | 不并行改；先完成 utils 路径，再集中改 workflow 入口 |
| `*_with_report` wrapper | 不主动处理 | Phase 3 删除或迁移 agent-facing wrapper | 入口计划独占；在 agent 决策后执行 |

推荐顺序：

1. 先执行 `current_reconstructions.md`：固定 utils / engines / CP2K 解析边界。
2. 再执行本计划：重构 `main.py`、`workflows/` 和 CLI 调用。

---

## 目标设计

### 目录结构

```text
src/md_analysis/
├── main.py
└── workflows/
    ├── __init__.py
    ├── models.py
    ├── water.py
    ├── potential.py
    ├── charge.py
    ├── calibration.py
    ├── enhanced_sampling.py
    ├── scripts.py
    └── composite.py
```

### 分层职责

| 文件 | 职责 |
|---|---|
| `main.py` | 公开 facade，只 import / re-export 稳定入口函数 |
| `workflows/models.py` | 统一返回类型，如 `WorkflowResult`、`ArtifactMap` |
| `workflows/water.py` | 水结构相关完整 workflow |
| `workflows/potential.py` | 电势、Fermi、电极电势、phi(z)、厚度敏感性 workflow |
| `workflows/charge.py` | Bader surface charge、tracked charge、counterion charge workflow |
| `workflows/calibration.py` | charge-potential calibration fit / predict workflow |
| `workflows/enhanced_sampling.py` | slow-growth、constrained TI、constant-potential correction workflow |
| `workflows/scripts.py` | Bader / TI / potential / DeePMD SP 工作目录生成 workflow |
| `workflows/composite.py` | 跨模块组合 workflow，例如界面分析总流程 |

### `main.py` 目标形态

`main.py` 最终应接近下面这种薄 facade：

```python
"""Public programmatic entry points for md_analysis workflows."""

from __future__ import annotations

from .workflows.models import WorkflowResult
from .workflows.water import run_water_three_panel
from .workflows.potential import run_potential_full
from .workflows.charge import (
    run_surface_charge,
    run_tracked_charge,
    run_counterion_charge,
)
from .workflows.calibration import (
    run_calibration_fit,
    run_calibration_predict,
)
from .workflows.enhanced_sampling import (
    run_slowgrowth_plot,
    run_ti_full_analysis,
    run_ti_constant_potential_correction,
)
from .workflows.scripts import (
    run_bader_batch,
    run_ti_batch,
    run_potential_batch,
    run_sp_batch,
)
from .workflows.composite import run_interface_analysis

__all__ = [
    "WorkflowResult",
    "run_water_three_panel",
    "run_potential_full",
    "run_surface_charge",
    "run_tracked_charge",
    "run_counterion_charge",
    "run_calibration_fit",
    "run_calibration_predict",
    "run_slowgrowth_plot",
    "run_ti_full_analysis",
    "run_ti_constant_potential_correction",
    "run_bader_batch",
    "run_ti_batch",
    "run_potential_batch",
    "run_sp_batch",
    "run_interface_analysis",
]
```

---

## API 命名原则

1. **函数名必须诚实**
   - 当前 `run_all` 实际只跑 water + potential，不应继续叫 `run_all`。
   - 推荐改为 `run_interface_analysis` 或 `run_water_potential_analysis`。

2. **不要保留 agent 语义**
   - 删除或迁走依赖旧 `main.py` report dataclass 的 `*_with_report`。
   - 脚本准备类 wrapper 如果已经是 contract-backed，不应按名字盲删；应迁移到 `workflows/scripts.py` 并改成普通 `run_*` 入口。
   - 如果结构化返回值有价值，所有 `run_*` 统一返回 `WorkflowResult`，而不是只给 agent wrapper 使用。

3. **CLI 与 notebook 共用入口**
   - CLI 只负责交互输入和参数收集。
   - CLI 命令执行时调用 `md_analysis.main.run_*` 或 `md_analysis.workflows.*.run_*`。
   - 不在 CLI 模块里复制业务流程逻辑。

4. **workflow 管输出契约**
   - 每个 workflow 模块负责自己的输出目录和 artifact keys。
   - `main.py` 不应该知道具体 CSV/PNG 文件名。

5. **底层模块保持科学计算职责**
   - `water`、`electrochemical`、`enhanced_sampling`、`scripts` 仍负责真实分析和生成逻辑。
   - `workflows` 只是把参数、输出目录、artifact、metadata 组织成稳定入口。

---

## 统一返回类型

建议先使用一个通用结果类型；简单 workflow 只返回 `WorkflowResult`，复杂 workflow 可以在 `extra` 中携带强类型 report：

```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class WorkflowResult:
    name: str
    output_dir: Path
    artifacts: dict[str, Path] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    extra: Any | None = None

    def to_dict(self) -> dict[str, object]:
        extra = (
            self.extra.to_dict()
            if hasattr(self.extra, "to_dict") else self.extra
        )
        return {
            "name": self.name,
            "output_dir": str(self.output_dir),
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "metadata": dict(self.metadata),
            "extra": extra,
        }
```

`artifacts` 用来放真实文件路径，例如 CSV、PNG、TXT。  
`metadata` 用来放轻量信息，例如帧数、method、是否跑了 Fermi、是否跑了 thickness sweep、输入模式等。

复杂 workflow，例如 TI full analysis，允许保留自己的强类型 dataclass，因为它的诊断字段会被上层用来区分失败类别。推荐规则：

- 简单 workflow：只使用 `WorkflowResult.metadata`。
- 复杂 workflow：返回 `WorkflowResult`，并把强类型 report 放入 `WorkflowResult.extra`。
- 强类型 report 自己提供 `to_dict()`，方便后续 agent / API 层重新接入。
- 不把所有复杂字段塞进无结构的 `dict[str, Any]`，避免丢失类型提示和 schema 线索。

---

## Agent 处理策略

本轮入口重构不维护旧 legacy agent 兼容性，但是否保留少数成熟 contract-backed task 需要在执行前拍板。

| 选项 | 含义 | 影响 |
|---|---|---|
| A | 删除整个 `agent/` 包 | 最干净；14 个 task、dispatch、TaskContract 全删 |
| B | 保留 contract-backed 的脚本/复合任务，删除依赖旧 `*_with_report` 的 legacy task | 推荐；保留已成熟资产，避免旧 report wrapper 拖住入口重构 |
| C | 全保留 agent，但 handler 全部改调 workflows | 工作量最大；本轮不建议 |

建议采用 **B**：保留已经独立成熟、且不依赖旧 `main.py` report dataclass 的 task；依赖 `*_with_report` 的 legacy task 直接删或停用。这样既不为旧 agent 兼容性牺牲架构，也不浪费已经成型的 contract-backed 工作。

需要特别区分：

- 分析类 `*_with_report`：作为旧 agent 边界删除或迁移。
- 脚本准备类 wrapper：迁到 `workflows/scripts.py`，改名为普通 `run_*`，返回 `WorkflowResult`。

---

## 建议公开入口

### Water

- `run_water_three_panel(...) -> WorkflowResult`

### Potential

- `run_center_potential(...) -> WorkflowResult`
- `run_fermi_energy(...) -> WorkflowResult`
- `run_electrode_potential(...) -> WorkflowResult`
- `run_phi_z_profile(...) -> WorkflowResult`
- `run_thickness_sensitivity(...) -> WorkflowResult`
- `run_potential_full(...) -> WorkflowResult`

### Charge

- `run_surface_charge(...) -> WorkflowResult`
- `run_tracked_charge(...) -> WorkflowResult`
- `run_counterion_charge(...) -> WorkflowResult`
- `run_charge_full(...) -> WorkflowResult`

### Calibration

- `run_calibration_fit(...) -> WorkflowResult`
- `run_calibration_predict(...) -> WorkflowResult`

### Enhanced Sampling

- `run_slowgrowth_quick_plot(...) -> WorkflowResult`
- `run_slowgrowth_publication_plot(...) -> WorkflowResult`
- `run_ti_single_diagnostics(...) -> WorkflowResult`
- `run_ti_full_analysis(...) -> WorkflowResult`
- `run_ti_constant_potential_correction(...) -> WorkflowResult`

### Scripts / Preparation

- `run_bader_single(...) -> WorkflowResult`
- `run_bader_batch(...) -> WorkflowResult`
- `run_ti_single(...) -> WorkflowResult`
- `run_ti_batch(...) -> WorkflowResult`
- `run_potential_single(...) -> WorkflowResult`
- `run_potential_batch(...) -> WorkflowResult`
- `run_sp_single(...) -> WorkflowResult`
- `run_sp_batch(...) -> WorkflowResult`

### Composite

- `run_interface_analysis(...) -> WorkflowResult`
- 可选：`run_all(...) -> WorkflowResult`，但只有在它真的覆盖完整分析面时才保留这个名字。

---

## Phase 0：确认当前入口使用面

**目的**：先知道哪些地方依赖 `main.py`，避免盲删。

- [ ] 扫描引用：`rg "md_analysis.main|from \\.\\.main|lazy_import\\(\"md_analysis.main\"" src test docs context4agent`
- [ ] 列出 CLI 当前调用 `main.py` 的菜单项
- [ ] 列出 test 当前直接 import `main.py` 的测试
- [ ] 标记 agent-only 引用，本轮不维护

**Acceptance**：知道哪些引用要迁移，哪些可以删除。

---

## Phase 1：新增 `workflows/models.py`

**目的**：先建立统一返回类型，不动现有 workflow 行为。

- [ ] 新建 `src/md_analysis/workflows/__init__.py`
- [ ] 新建 `src/md_analysis/workflows/models.py`
- [ ] 添加 `WorkflowResult`
- [ ] 添加 artifact 校验 helper，例如 `require_artifacts_exist(result)`
- [ ] 添加最小单元测试

**Acceptance**：新类型可 import；不影响旧入口。

---

## Phase 2：迁移非 agent 的 leaf workflows

**目的**：把当前 `main.py` 里的普通入口拆到对应模块。

- [ ] `run_water_analysis` 迁到 `workflows/water.py`，改名 `run_water_three_panel`
- [ ] `run_potential_analysis` 迁到 `workflows/potential.py`，改名 `run_potential_full`
- [ ] `run_charge_analysis` 迁到 `workflows/charge.py`，改名 `run_surface_charge`
- [ ] `run_tracked_charge_analysis` 迁到 `workflows/charge.py`，改名 `run_tracked_charge`
- [ ] `run_counterion_charge_analysis` 迁到 `workflows/charge.py`，改名 `run_counterion_charge`
- [ ] 返回值逐步改成 `WorkflowResult`
- [ ] `main.py` 暂时 re-export 新函数

**Acceptance**：water / potential / charge 的主要编程入口从 `workflows` 提供；`main.py` 变薄。

---

## Phase 3：删除 agent-facing wrappers

**目的**：不再维护依赖旧 `main.py` report dataclass 的 agent-facing wrapper。

- [ ] 删除 `ChargeSurfaceReport`
- [ ] 删除 `TrackedChargeReport`
- [ ] 删除 `CounterionChargeReport`
- [ ] 删除 `WaterThreePanelReport`
- [ ] 删除 `PotentialFullReport`
- [ ] 删除 `RunAllReport`
- [ ] 删除或迁移分析类 `*_with_report`
- [ ] 脚本准备类 `*_with_report` 不按名字盲删；迁到 `workflows/scripts.py` 并改成普通 `run_*`
- [ ] 删除或停用依赖旧 report wrapper 的 agent contract tests

**Acceptance**：`main.py` 和 `workflows` 中不再出现 agent-facing 命名。

---

## Phase 4：补齐 CLI 已有菜单对应的 workflow 入口

**目的**：让程序化入口覆盖 CLI 的主要能力。

- [ ] `workflows/calibration.py`
  - `run_calibration_fit`
  - `run_calibration_predict`
- [ ] `workflows/enhanced_sampling.py`
  - slow-growth quick / publication plot
  - TI single diagnostics
  - TI full analysis
  - constant-potential correction
- [ ] `workflows/scripts.py`
  - Bader single / batch
  - TI single / batch
  - potential single / batch
  - SP single / batch
- [ ] 每个入口返回 `WorkflowResult`

**Acceptance**：CLI 的主要菜单项都有对应程序化入口。

---

## Phase 5：重命名 composite workflow

**目的**：修正 `run_all` 名称不准确的问题。

- [ ] `run_all` 改为 `run_interface_analysis` 或 `run_water_potential_analysis`
- [ ] 放入 `workflows/composite.py`
- [ ] 明确它只组合 water + potential，除非显式增加 charge / calibration / TI
- [ ] 暂不把 docs 里的五条主线 workflow 全部做成 composite；先完成 leaf workflow 和 `run_interface_analysis`
- [ ] 第二轮可评估是否增加公开 composite：
  - cSHE potential
  - surface charge + calibration
  - SG to TI
  - Bader batch
  - SP for DP
- [ ] 如果需要短期迁移，可在 `main.py` 保留 `run_all = run_interface_analysis`，并标注 deprecation
- [ ] 文档同步改名

**Acceptance**：不存在“名字像全量分析，实际只跑一部分”的入口。

---

## Phase 6：让 CLI 调用 workflows/main facade

**目的**：CLI 只负责交互和参数收集，业务执行统一走 workflow 入口。

- [ ] 检查 `src/md_analysis/cli/_water.py`
- [ ] 检查 `src/md_analysis/cli/_potential.py`
- [ ] 检查 `src/md_analysis/cli/_charge.py`
- [ ] 检查 `src/md_analysis/cli/_calibration.py`
- [ ] 检查 `src/md_analysis/cli/_enhanced_sampling.py`
- [ ] 检查 `src/md_analysis/cli/_scripts.py`
- [ ] 将菜单命令执行逻辑统一改为调用 `md_analysis.main` 或 `md_analysis.workflows`
- [ ] 避免 CLI 直接知道底层 CSV/PNG 细节

**Acceptance**：CLI 和 notebook/script 用户走同一批 workflow 函数。

---

## Phase 7：瘦身 `main.py`

**目的**：把 `main.py` 固定为稳定公开 facade。

- [ ] 删除具体 workflow 实现
- [ ] 删除具体 report dataclass
- [ ] 只保留 imports、`__all__`、简短模块 docstring
- [ ] 视情况保留短期兼容 alias
- [ ] 补充 `src/md_analysis/CLAUDE.md` 中的入口说明

**Acceptance**：`main.py` 不再超过约 100 行；新入口清楚可读。

---

## Phase 8：测试与文档同步

- [ ] 更新 `test/integration/test_main.py`
- [ ] 新增或更新 `test/unit/workflows/`
- [ ] 更新 `docs/workflows.md`
- [ ] 更新 `docs/README.md`
- [ ] 更新 `context4agent/architecture/modules/src/interface_exposure.md`
- [ ] 更新 `context4agent/requirements/short_term.md`
- [ ] 跑 `pytest test/unit -q`
- [ ] 可行时跑 `pytest test/integration -q`

**Acceptance**：文档、测试、入口函数名一致。

---

## 迁移风险

| 风险 | 处理方式 |
|---|---|
| `main.py` 被 CLI lazy import | 先让 `main.py` re-export 新函数，再逐步改 CLI |
| 旧测试直接 import `run_*_analysis` | Phase 2 后保留短期 alias，Phase 8 统一改测试 |
| `run_all` 名称变更影响文档 | Phase 5 集中改名，并解释新 composite 范围 |
| 返回类型从 `dict[str, Path]` 改为 `WorkflowResult` | 可先兼容 `WorkflowResult.artifacts`；必要时提供 `to_dict()` |
| scripts 类 workflow 有外部副作用 | 明确只生成工作目录，不提交 PBS/SLURM 作业；提交作业必须另设确认入口 |
| CLI 改造与 utils 路径改造撞车 | 先执行 `current_reconstructions.md`；本计划 Phase 6 再集中改 CLI |
| `WorkflowResult` 与 legacy dict 迁移期共存 | 在 `WorkflowResult` 上临时加 `from_legacy_dict()`，入口迁完后删除 |

---

## 建议验收命令

```bash
rg "with_report|Agent-safe|agent-facing" src/md_analysis/main.py src/md_analysis/workflows
rg "run_all" src test docs context4agent
rg "md_analysis.main" src test docs context4agent
pytest test/unit -q
pytest test/integration -q
```

第一条命令应无结果，表示主入口不再带 agent 语义。  
第二条命令用于确认 `run_all` 是否已改名或明确保留。  
第三条命令用于检查哪些位置仍依赖 `main.py` facade。

---

## 建议暂停点

- Phase 2 后：普通 water / potential / charge 入口已经迁入 `workflows`
- Phase 3 后：agent wrapper 已清除，`main.py` 体积会明显下降
- Phase 5 后：composite 命名修正，API 语义更清楚
- Phase 7 后：`main.py` 变成薄 facade，适合做一次整体 review
