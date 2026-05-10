# 实现 + 注册 + 测试 + 验收

[← 回 README](README.md) | [← interface_contract.md](interface_contract.md)

---

## 1. 文件放哪

建议两个新文件：

```
src/md_analysis/utils/RestartParser/VaspParser.py   # 解析 INCAR / ICONST / REPORT
src/md_analysis/enhanced_sampling/_parsers_vasp.py  # VASPParser 类 + 注册
```

**为什么这样分**：

- `VaspParser.py`（在 utils/RestartParser 下）只负责"把 VASP 文件 → ColvarRestart / LagrangeMultLog"，是**纯解析 + 单位换算**，可独立单元测试
- `_parsers_vasp.py`（在 enhanced_sampling 下）薄类 + 注册，跟现有 `_parsers.py` 平级

也可以全塞 `_parsers_vasp.py`，但解析逻辑稍多就拆开比较清爽。

---

## 2. VASPParser 实现骨架

```python
# src/md_analysis/enhanced_sampling/_parsers_vasp.py

from pathlib import Path

from ..utils.RestartParser.ColvarParser import ColvarRestart, LagrangeMultLog
from ..utils.RestartParser.VaspParser import (   # 你写的
    parse_vasp_metadata,
    parse_vasp_report,
    find_vasp_files,
)
from ._parsers import register_parser


class VASPParser:
    name = "vasp"

    def is_constraint_directory(self, directory: Path) -> bool:
        if not directory.is_dir():
            return False
        try:
            find_vasp_files(directory)  # 检查 INCAR + ICONST + REPORT 都在
        except FileNotFoundError:
            return False
        return True

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        files = find_vasp_files(directory)
        return parse_vasp_metadata(
            incar=files.incar,
            iconst=files.iconst,
            report=files.report,  # 用第一行拿当前 ξ
        )

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
        files = find_vasp_files(directory)
        return parse_vasp_report(files.report)


# 模块导入时自动注册
register_parser("vasp", VASPParser)
```

`parse_vasp_metadata` / `parse_vasp_report` 内部要做 [interface_contract.md §3](interface_contract.md#3-单位转换--关键章节) 列出的单位换算。

---

## 3. 让 import 触发注册

`register_parser("vasp", VASPParser)` 必须**真正执行过**才能让 sniffer 识别 VASP 数据。两种做法二选一：

### 做法 A：用户显式 import（不强制注册）

用户调用前显式：

```python
from md_analysis.enhanced_sampling._parsers_vasp import VASPParser  # 触发注册
```

或调 CLI 时不需要做任何事——但要求 CLI 启动时 import 这个模块。

### 做法 B：`enhanced_sampling/__init__.py` 顶层 import（自动注册）

在 `src/md_analysis/enhanced_sampling/__init__.py` 加一行：

```python
from . import _parsers_vasp  # noqa: F401  — auto-register VASP parser
```

代价：每次 `import md_analysis.enhanced_sampling` 就把 VASPParser 加载进 registry，启动稍慢一点点但可忽略。

**推荐做法 B**，对用户最透明。

---

## 4. 单元测试模板

新建 `test/unit/enhanced_sampling/test_parsers_vasp.py`：

```python
"""Unit tests for VASPParser."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.enhanced_sampling._parsers import (
    ConstraintMDParser,
    infer_parser,
)
from md_analysis.enhanced_sampling._parsers_vasp import VASPParser


REPO_ROOT = Path(__file__).resolve().parents[3]
VASP_FIXTURE = REPO_ROOT / "data_example" / "vasp_ti" / "ti_target_0.5"

_skip = pytest.mark.skipif(
    not VASP_FIXTURE.is_dir(),
    reason="VASP TI fixture missing — add data_example/vasp_ti/",
)


# ---------- Protocol conformance ----------

def test_satisfies_protocol():
    parser = VASPParser()
    assert isinstance(parser, ConstraintMDParser)
    assert parser.name == "vasp"


# ---------- Real-data parsing ----------

@_skip
def test_recognises_directory():
    assert VASPParser().is_constraint_directory(VASP_FIXTURE)


@_skip
def test_metadata_target_keeps_native():
    """target_au keeps VASP native value, no Å→Bohr / deg→rad
    conversion. Algorithm layer is unit-agnostic on ξ."""
    meta = VASPParser().parse_metadata(VASP_FIXTURE)
    # If REPORT first row's cv column reads 0.5 (Å in VASP native),
    # target_au must store 0.5 — NOT 0.5/BOHR_TO_ANG.
    assert meta.colvars.primary.target_au == pytest.approx(0.5, abs=1e-3)
    assert meta.timestep_fs > 0
    # TI: TARGET_GROWTH should be zero (constraint fixed).
    assert meta.colvars.primary.target_growth_au == pytest.approx(0.0)


@_skip
def test_lambda_series_shape():
    log = VASPParser().parse_lambda_series(VASP_FIXTURE)
    assert log.n_steps > 0
    assert log.collective_shake.shape == (log.n_steps,)


@_skip
def test_lambda_sign_and_scale():
    """shake = -g_vasp / HA_TO_EV  (Hartree / ξ_native).
    Numerical scale check: |shake| median should fall in
    1e-5 .. 1e-1 Hartree-per-(ξ_native), assuming chemistry-typical
    g_vasp ~ 0.001 .. 5 eV/(ξ_native)."""
    log = VASPParser().parse_lambda_series(VASP_FIXTURE)
    median_abs = np.median(np.abs(log.collective_shake))
    assert 1e-5 < median_abs < 1.0, (
        f"|shake| median = {median_abs:.3e} a.u. — likely missing "
        "÷ HA_TO_EV or wrong REPORT column."
    )


# ---------- Negative cases ----------

def test_rejects_empty_dir(tmp_path):
    assert VASPParser().is_constraint_directory(tmp_path) is False


# ---------- Sniffer integration ----------

@_skip
def test_infer_parser_picks_vasp():
    parser = infer_parser(VASP_FIXTURE)
    assert parser.name == "vasp"
```

跑：

```bash
pytest test/unit/enhanced_sampling/test_parsers_vasp.py -v
```

---

## 5. 端到端测试

### 5.1 TI 全流程

```python
from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    run_ti_full_from_root,
)

# parser="auto" 会嗅探出 vasp（前提：含 VASP 文件）
report = run_ti_full_from_root(
    root_dir="path/to/vasp_ti_root/",
    output_dir="output/",
    epsilon_tol_ev=0.01,
    equilibration=500,
)

print(f"ΔA = {report.delta_A_eV:+.4f} ± {report.sigma_A_eV:.4f} eV")
```

或 CLI 菜单 312：

```
$ md-analysis
 Input: 312
 TI root directory: path/to/vasp_ti_root/
 ...
```

输出目录结构跟 CP2K 完全一样（`output/enhanced_sampling/constrained_ti/...`）。

### 5.2 SG 全流程

```python
from md_analysis.enhanced_sampling.slowgrowth.SlowGrowth import SlowgrowthFull

sg = SlowgrowthFull.from_directory("path/to/vasp_sg_run/")  # parser="auto"
print(f"Final ΔA = {sg.free_energy_au[-1] * 27.2114:.4f} eV")
```

或 CLI 菜单 301 / 302。

### 5.3 关键 sanity check：ΔA 数量级

跑你最熟悉的化学反应（一个 SG 或一组 TI 点），看 ΔA 是不是落在合理范围（~ 0.1 - 5 eV）。如果跑出 ΔA = 27 eV 或 0.001 eV，**单位换算错了**——回 [interface_contract.md §3](interface_contract.md#3-单位转换--关键章节) 复查。

---

## 6. 验收标准

### 必过测试

1. **Protocol conformance**：`isinstance(VASPParser(), ConstraintMDParser)` → True
2. **元数据正确**：`timestep_fs == POTIM`；`target_au` 跟 ICONST 给的目标值经过单位换算后一致
3. **λ 时序形状正确**：单 CV 时 `(N,)`；多 CV 时 `(N, K)`
4. **λ 单位对**：`np.median(np.abs(shake))` 落在 1e-5 ~ 10 a.u. 之间（meV/Å 量级换算后的合理值）

### 必过 sanity check

5. **TI 全流程跑通**：菜单 312 / `run_ti_full_from_root` 不报错，CSV / PNG 都生成
6. **SG 全流程跑通**：菜单 301 / `SlowgrowthFull.from_directory` 不报错
7. **ΔA 数量级合理**：在你预期的化学反应范围（典型 0.1 - 5 eV）
8. **CP2K 路径不破**：现有 `pytest test/unit/enhanced_sampling/` 全绿（708 → 708+，零回归）

### 加分项

9. **数值交叉验证**：CP2K vs VASP 同体系同 functional 跑同一个反应坐标，每点 ⟨λ⟩ 量级一致（不强求数值完全一致——不同 functional / 不同 PP 不可能完全一致；量级对、形状趋势对就行）
10. **Multi-CV 支持**：`n_constraints > 1` 的 ICONST 也能正确解析
11. **ICONST 多种 CV 类型**：distance / angle / coordination 都过

---

## 7. PR 流程

1. fork `StackOverFlow11/MD-analysis`
2. 从 `development_feature_agent` 分支 checkout 一个新分支：`feat/vasp-io-adapter`
3. 实现 + 测试
4. 跑全套 unit suite：`pytest test/unit/ -q` 必须 ≥ 708 全绿
5. PR 标题：`feat(enhanced_sampling): add VASPParser for engine-agnostic TI/SG`
6. PR 描述里说清楚：
    - **VASP 版本**（5.x / 6.x，及具体子版本号）
    - **测过的体系**（什么反应、几个 TI 点 / 多长 SG / 用了哪种 CV 类型）
    - **REPORT 文件具体哪一列**当作 λ（列号 + 列名 + 单位）
    - **符号方向确认**（你怎么验证的 dA/dξ 正负方向跟 CP2K 一致）
    - **λ 换算公式**（应为 `shake = -g_vasp / HA_TO_EV`；`HA_TO_EV` 从 `md_analysis.utils.constants` import）
    - 任何 sanity check 数值（ΔA 跟手算 / CP2K 对比）

---

## 8. 你需要从我（fenglin）那里拿的

理论上你只要 README + interface_contract + 这份 integration_steps 三份文档就够了。但你可能需要：

- **空一台 chem-hpc 的账号 / 测试时间**：跑一两个 VASP TI demo 当 fixture 落到 `data_example/vasp_*/`
- **VASP REPORT 真实样本**：我手头没有；你跑出来后给我两份（TI 一份、SG 一份），我归档到 `data_example/vasp_*/`
- **PBS 提交脚本模板**：现有 CP2K / VASP 模板可参考 `~/.config/md_analysis/config.json` 里的 `cp2k_script_path` / `vasp_script_path`

---

## 9. 交付时间建议

不限。但分阶段：

- **第 1 周**：把 `VaspParser.py` 解析函数写出来，跑通 unit tests（用 mock REPORT 字符串，不需要真 VASP 数据）
- **第 2 周**：跑一个真 VASP TI demo，过 end-to-end + ΔA sanity check
- **第 3 周**：PR + review

中途 stuck 发邮件 fenglinshao02@gmail.com 或开 GitHub issue。

---

[← 回 README](README.md) | [← interface_contract.md](interface_contract.md)
