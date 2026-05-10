# vasp_io — 给 VASP 适配器实现者的交付规格

## 这是什么

md-analysis 现有的 TI（Thermodynamic Integration）和 SG（Slow-Growth）分析模块完全可以复用，**只要你能把 VASP 的约束 MD 输出文件解析成两个标准 dataclass**。

写完这部分，TI 全流程（4 步收敛诊断 + ΔA 积分 + 恒电势修正）和 SG（midpoint 累积积分 + 双轴绘图）都自动支持 VASP 数据，**md-analysis 计算层 / 算法层一行不动**。

---

## 一个核心 framing：你只需要交付一对**共轭量**

自由能积分本质上是：

```
ξ        : 反应坐标，任意单位 U（Å / Bohr / 弧度 / 无量纲 CN ...）
dA/dξ    : 自由能梯度，单位 = 能量/U
∫ dA/dξ × dξ → 能量      # U 自动消掉
```

**算法层完全 unit-agnostic**：

| 模块 | 对单位 / engine 的依赖 |
|---|---|
| ACF 自洽截断（Sokal 1997） | 零，纯 numpy |
| Flyvbjerg-Petersen block average | 零 |
| running-average drift / Geweke | 零 |
| 梯形积分 + SEM 传播 | 零 |
| SG midpoint 积分 | 微弱（硬编码符号 + per-step 假设） |

也就是说：**ξ 是什么单位、CN 用 CP2K 的有理函数还是 VASP 的 fermi-dirac switching、λ 来自哪个引擎——算法都不在乎**。算法只要拿到一组 (ξ_array, λ_array) 在同一组反应坐标上、彼此共轭单位匹配，就能跑出正确的 ΔA。

唯一的"软约束"：**最终 ΔA 数值要落在 eV 量级**，不然下游绘图轴标签 / CSV 列名带的"eV"标签会撒谎。

要做到这点你只需要做**一件事**：把 λ 从 `eV/ξ_native` 换成 `-Hartree/ξ_native`（`λ = -g_vasp / HA_TO_EV`），然后 ξ 保留 VASP native 数值不动。这样 ∫λdξ 自然产出 -Hartree，输出层 × HA_TO_EV 把符号吃回去得到 +eV。完整推导见 [interface_contract.md §3.4](interface_contract.md#34-端到端自洽性验证写完-parser-跑过一次确认)。

> 换句话说：dataclass 字段名带 `_au`（atomic unit）是历史命名习惯，**不是单位依赖**。`target_au` 字段填 VASP native 数值（Å / deg / CN），算法层照样跑对。

---

## 你不需要担心的事

读完上面那段你应该想到：

| 担忧 | 实际情况 |
|---|---|
| "VASP 的 CN switching 跟 CP2K 不同，ξ 不可比" | 这是**物理问题**，不是数据层问题。同一引擎内部 ξ 自洽即可。跨引擎 ΔA 形状不可直接对比，但 invariant（如总 ΔA、TS 位置）可以——这是用户做物理分析时要意识到的事。 |
| "VASP 多 CV 组合（线性叠加）量纲乱" | 同上。VASP 内部把组合 CV 当成"单一抽象 ξ"看待 + 输出对应的"抽象 λ"，工具不需要知道里面的拆分。 |
| "VASP λ 单位是 eV，跟 CP2K Hartree 不一样" | parser 里做一次性换算 `λ = -g_vasp / HA_TO_EV`（除以常数 + 翻符号），算法层不感知。详见 [interface_contract.md §3](interface_contract.md#3-单位换算--只动-λ不动-ξ)。 |
| "VASP ξ 是 Å / deg / 无量纲，要不要换成 Bohr / 弧度才能跟 CP2K 同尺度？" | **不用**。算法层 unit-agnostic，ξ 保留 VASP native 数值，∫λdξ 自动消单位。 |
| "VASP 没法填 ColvarRestart 的某些字段（如 PROJECT_NAME）" | 那些字段是 metadata，TI/SG 真正用到的字段很少，[interface_contract.md](interface_contract.md) 里有最小必填集 + 兜底值表。 |
| "CN 是无量纲，dA/dξ 怎么积出能量" | 无量纲也是单位，∫(eV) × (dimensionless) = eV，自然成立。 |

---

## 阅读顺序

1. **[interface_contract.md](interface_contract.md)** — 完整接口规范
    - `ConstraintMDParser` Protocol 三方法签名
    - `ColvarRestart` 元数据 dataclass 字段表（每个字段的物理含义、最小必填集、兜底值）
    - `LagrangeMultLog` 约束力时序 dataclass
    - λ 的唯一换算公式（共轭量原则）
    - **不需要做的事**：ξ 单位推断 / 算法层 / 字段名 / 4 步诊断 / 积分逻辑全部不要碰

2. **[integration_steps.md](integration_steps.md)** — 实现 + 测试 + PR 流程
    - 文件放哪、注册策略
    - 单元测试模板（含 mock 不依赖真 VASP 数据）
    - 端到端 sanity check（跑通 TI 312 / SG 301）
    - 验收标准

---

## 你需要的背景

**TI / SG 物理直觉**（如果还没接触过）：

- **TI**：在反应坐标 ξ 上选 K 个固定点，每个点跑一段约束 MD（CV 卡死 = ξ_k），统计约束力 ⟨λ_k⟩，最后 ΔA = -∫⟨λ⟩dξ ≈ -Σ w_k⟨λ_k⟩。
- **SG**：CV 线性匀速从初态拉到末态（一段长 traj），每步累积 ∫λdξ，得 ΔA(ξ) 曲线。

**两者只需要一个核心物理量**：约束力 λ(t) 时间序列。其他都是元数据（dt, ξ₀, dξ/dt, t₀）。

**VASP 约束 MD 端**（请你自己再核对手册细节）：

- VASP 用 `ICONST` 文件定义约束（蓝月集合 / Blue Moon ensemble）—— **这是你跑 VASP 的事，本工具 parser 不需要解析 ICONST**
- 约束力 λ + cv(t) + step + time 都写在 `REPORT` 文件，每个 MD step 一行
- **本工具 parser 只读 REPORT 一个文件**就够（dt 从行间时间差推出、target 从第一行 cv 列拿、growth 从前两行 cv 差推）；不依赖 INCAR 的 POTIM 等

> ⚠️ REPORT 文件每一列的具体含义、单位、符号约定，不同 VASP 版本可能有出入。**请以你实际跑的 VASP 版本手册为准**，并在 PR 描述里说明你 parser 用的是哪一列作为 λ。

---

## 上下文链接

- 现有 Protocol 实现：`src/md_analysis/enhanced_sampling/_parsers.py`（CP2K 实现可作参考，不到 100 行）
- TI 入口：`src/md_analysis/enhanced_sampling/constrained_ti/`
- SG 入口：`src/md_analysis/enhanced_sampling/slowgrowth/SlowGrowth.py`
- dataclass 定义：`src/md_analysis/utils/RestartParser/ColvarParser.py`
- 中文用户文档：`docs/`（CLI 用法、5 条主线工作流、菜单 reference、排错）
- TI/SG IO 重构记录：commit `4541369`（"engine-agnostic IO via ConstraintMDParser Protocol"）

有问题随时找 fenglinshao02@gmail.com。
