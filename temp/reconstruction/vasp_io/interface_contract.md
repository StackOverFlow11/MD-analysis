# 接口契约

[← 回 README](README.md)

下面"必填"指 TI/SG 算法实际会用到的字段；其他是 metadata，填合理占位即可。

---

## 1. Protocol：你要实现的三个方法

定义在 `src/md_analysis/enhanced_sampling/_parsers.py`：

```python
from typing import Protocol, runtime_checkable

@runtime_checkable
class ConstraintMDParser(Protocol):
    name: str  # class attribute, lowercase identifier (e.g. "vasp")

    def is_constraint_directory(self, directory: Path) -> bool:
        """Return True iff *directory* contains the files this parser
        recognises. Cheap (file existence check only, no content parse)."""

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        """Parse engine metadata: dt, target value, growth rate, etc.
        Cheap (KB-sized file). Used at discovery time to sort points by ξ."""

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
        """Parse the full constraint-force time series.
        Heavy. Called only when analysis actually needs the data."""
```

**关键约束**：
- 三个方法是 Protocol 接口（duck typing），不需要继承基类
- 参数 `directory` 是 `pathlib.Path`，单个约束点目录
- `name` 必须是小写、纯 ASCII（用作 registry key）
- **方法不能有副作用**（不写盘、不打印、不改环境）；纯读 + 解析

---

## 2. 数据契约：你要填的两个 dataclass

两个 dataclass 都在 `src/md_analysis/utils/RestartParser/ColvarParser.py`，都是 `frozen=True`，**字段定义不要改**——你 import 用就行。

### 2.1 `ColvarRestart`（元数据）

```python
@dataclass(frozen=True)
class ColvarRestart:
    project_name: str
    step_start: int
    time_start_fs: float
    timestep_fs: float
    total_steps: int
    colvars: ColvarInfo
    lagrange_filename: str | None
    cell_abc_ang: tuple[float, float, float]
    fixed_atom_indices: tuple[int, ...] | None
```

#### 必填字段（TI/SG 算法直接用到）

**这些字段都从 REPORT 一个文件就能拿出来**，不需要解析 INCAR / ICONST。

| 字段 | 物理含义 | 数值 | REPORT 来源 |
|---|---|---|---|
| `step_start` | 起始 MD 步号 | 整数 | 第一行 step 列 |
| `time_start_fs` | 起始时间 (fs) | float | 第一行 time(fs) 列 |
| `timestep_fs` | dt (fs) | float | 后行 time − 前行 time |
| `colvars.primary.target_au` | 当前 ξ 值 | **保留 VASP native 数值，不换算** | 第一行 cv 列 |
| `colvars.primary.target_growth_au` | dξ/dt 增长率 | per-(a.u. time)，见 §3.3 | 前两行 cv 差 + 时间差 |

#### 可填占位字段

| 字段 | 兜底值 | 影响 |
|---|---|---|
| `project_name` | `""` 或目录名 | 仅 metadata |
| `total_steps` | 0 | 仅 metadata（用户提示用） |
| `lagrange_filename` | `None` | 不影响 |
| `cell_abc_ang` | `(0.0, 0.0, 0.0)` 或从 POSCAR 读 | TI/SG 不读；恒电势修正模块（313）会读 |
| `fixed_atom_indices` | `None` | TI/SG 不读 |

### 2.2 `ConstraintInfo`（CV 信息）

```python
@dataclass(frozen=True)
class ConstraintInfo:
    colvar_id: int
    target_au: float
    target_growth_au: float
    intermolecular: bool

@dataclass(frozen=True)
class ColvarInfo:
    constraints: tuple[ConstraintInfo, ...]
    # 必须至少有一个；.primary 返回 constraints[0]
```

| 字段 | 物理含义 | 备注 |
|---|---|---|
| `colvar_id` | CV 编号；从 1 开始 | 单 CV 直接填 1 |
| `target_au` | 当前 CV 目标值 ξ | 见"单位转换" |
| `target_growth_au` | dξ/dt 增长率 | TI = 0；SG ≠ 0 |
| `intermolecular` | 是否分子间约束 | TI/SG 不读，填 `False` |

### 2.3 `LagrangeMultLog`（约束力时序）

```python
@dataclass(frozen=True)
class LagrangeMultLog:
    shake: np.ndarray
    rattle: np.ndarray
    n_steps: int
    n_constraints: int
```

**形状规则**：

| n_constraints | shake.shape | rattle.shape |
|---|---|---|
| 1 | `(N_steps,)` | 同左 |
| > 1 | `(N_steps, K)`，第一列是 primary CV | 同左 |

| 字段 | 含义 | 备注 |
|---|---|---|
| `shake` | 主约束力 λ(t) | 见"单位转换" |
| `rattle` | 速度修正阶段乘子 | VASP 没有 RATTLE 概念，填 `np.zeros_like(shake)` |
| `n_steps` | 时序长度 | |
| `n_constraints` | CV 数 | |

---

## 3. 单位换算 — 只动 λ，不动 ξ

### 3.1 共轭量原则

```
ξ        : 反应坐标，任意单位 U（VASP native）
dA/dξ    : 自由能梯度，单位 = 能量 / U
∫ dA/dξ × dξ → 能量      # U 自动消掉
```

算法层 unit-agnostic，**根本不知道 ξ 是 Å 还是 deg 还是 CN**。所以 parser 也不应该 "聪明地" 把 ξ 从 Å 换成 Bohr——这种"修正"不仅没必要，还会引入换算 bug 和跨 CV 类型的 if/else 分支。

**正确做法**：

> ξ 保留 VASP REPORT 给你的 native 数值，不动。
> 只把 λ（约束力 / 自由能梯度）做一次性换算，让 ∫λdξ 的输出落在能量 a.u.（Hartree）量级。

### 3.2 λ 的统一换算公式（唯一一条规则）

设 VASP REPORT 给出的自由能梯度（eV / ξ_native），记作 `g_vasp`：

```
shake[t] = - g_vasp[t] / HA_TO_EV
         单位: -Hartree / ξ_native
```

两个动作：

| 动作 | 为什么 |
|---|---|
| 除以 `HA_TO_EV`（= 27.211386245988） | eV → Hartree，让 ∫λdξ 直接输出 Hartree 数值 |
| 取负号 | TI/SG 内部最后还会取一次反（`dA/dξ = -⟨shake⟩`，CP2K 约定）；提前反一次让两次反号抵消 → 最终 dA/dξ 数值 = +g_vasp/HA_TO_EV，符号正确 |

### 3.3 时间相关字段

| 字段 | 公式 | 备注 |
|---|---|---|
| `timestep_fs` | REPORT 后行 time − 前行 time | fs，直接用 |
| `time_start_fs` | REPORT 第一行 time | fs，直接用 |
| `target_growth_au` | `(dξ_native_per_fs) × AU_TIME_TO_FS` | 字段语义是 "per a.u.-time"；TI 时 = 0；SG 时 = REPORT 前两行 cv 差除以两行 time 差，再 × AU_TIME_TO_FS |

`AU_TIME_TO_FS = 0.02418884326585` 在 `md_analysis.utils.constants`。

> 字段名带 `_au` 是历史命名习惯，不代表 ξ 必须用 atomic unit。`target_au` 字段实际语义是"primary CV 的当前值"，单位由你 parser 选定（推荐 native）。`target_growth_au` 字段语义是"per a.u.-time"——这里的 a.u. 指**时间单位**，不是 ξ 单位，所以这个 _au 是真要乘 AU_TIME_TO_FS 的。

### 3.4 端到端自洽性验证（写完 parser 跑过一次确认）

| 步 | 量 | 单位 | 来源 |
|---|---|---|---|
| 1 | shake | -Hartree / ξ_native | parser 填 |
| 2 | TI 内部：dA/dξ = -⟨shake⟩ | +Hartree / ξ_native | algorithm |
| 3 | ΔA = ∫(dA/dξ) dξ | Hartree | algorithm |
| 4 | CSV / 绘图：ΔA × HA_TO_EV | eV | output layer |

只要你 §3.2 那一行换算对了，1 → 4 全链路单位自洽，ξ 是什么单位完全不参与。

### 3.5 一个 corner case：VASP REPORT 是否已经是 dA/dξ？

不同 VASP 版本 REPORT 列约定可能不同。有的列是 `lambda`（拉格朗日乘子原始值），有的列是 `|z|^(1/2)*(lambda+GkT)`（带 mass-metric tensor 校正后的实际 dA/dξ）。**TI/SG 想要的是后者**。

确认方法：找一个简单 case（纯距离 CV，跑一段约束 MD），看 ⟨g_vasp⟩ 量级。化学反应 dA/dξ 量级一般 0.01 ~ 5 eV/(ξ_native 单位)。如果你拿到的列数值离谱（< 1e-4 或 > 100），说明拿错列了。

> 错了 ΔA 整体反号或量级偏，肉眼能看出来。**PR 描述里写明 VASP 版本 + 你确认过的列号 + 列名**。

---

## 4. CP2K 实现作参考

`src/md_analysis/enhanced_sampling/_parsers.py` 里 `CP2KParser` 不到 100 行：

```python
class CP2KParser:
    name = "cp2k"

    def is_constraint_directory(self, directory: Path) -> bool:
        try:
            self._find_restart(directory)
            self._find_log(directory)
        except FileNotFoundError:
            return False
        return True

    def parse_metadata(self, directory: Path) -> ColvarRestart:
        return parse_colvar_restart(self._find_restart(directory))

    def parse_lambda_series(self, directory: Path) -> LagrangeMultLog:
        return parse_lagrange_mult_log(self._find_log(directory))
```

`parse_colvar_restart` / `parse_lagrange_mult_log` 是 CP2K 文件具体解析（regex 解析 `&MD` / `&CONSTRAINT` / `&COLLECTIVE` 块）。你写 `parse_vasp_metadata` / `parse_vasp_report` 时**不需要**复用——结构跟 CP2K 文件完全不同。

---

## 5. 你**不**需要碰的东西（重要）

| 模块 | 不要碰 |
|---|---|
| 算法层 4 步诊断（ACF / F&P / running avg / Geweke） | 全部 unit-agnostic |
| TI 梯形积分 + SEM 传播 | unit-agnostic |
| SG midpoint 积分 | 跟 timestep_fs / target_growth_au 数值匹配即可 |
| dataclass 字段名 / `_au` 后缀 | 历史约定，不是单位依赖 |
| CSV 列名 / 绘图标签 | 写死了 "eV" / "a.u." 标签——你 parser 让数值落对了，标签自动正确 |
| `discover_ti_points` / `load_ti_series` | 完全 engine-agnostic |
| TaskContract / agent dispatch | 自动支持新 parser，不动 |

---

[← 回 README](README.md) | [→ integration_steps.md](integration_steps.md)
