# md-analysis 用户文档

面向使用者的 CLI 操作手册。如果你拿到 CP2K MD 数据想跑标准分析（水分布、电极电势、表面电荷、自由能等），**从这里开始**。

> 这一份是中文文档；项目根的 `README.md` 是英文的包介绍 / 架构概览。两份不冲突，互补。

---

## 阅读顺序

1. **[quickstart.md](quickstart.md)** — 5 分钟跑通第一个分析（水三联图）。先看这个。
2. **[workflows.md](workflows.md)** — 按"想做什么"组织的 5 条主线工作流：
    - 算电极电势（cSHE）
    - 表面电荷 σ + φ 标定
    - Slow-Growth 粗扫 → TI 精算
    - Bader 批量生成 VASP 工作目录
    - SP 单点为 DeePMD 训练集抽帧
3. **[menu_reference.md](menu_reference.md)** — 每个菜单代码（101 / 211 / 311 ...）的精简 reference：用途 + 必填参数 + 输出位置。
4. **[settings.md](settings.md)** — `9xx` 设置菜单 + `~/.config/md_analysis/config.json`。
5. **[pitfalls.md](pitfalls.md)** — 输入数据约定（帧目录命名 / interface 标签 / 单位）+ 常见错误排查。

---

## 核心概念速查

| 名词 | 含义 |
|---|---|
| **菜单代码** | `101` / `213` / `412` 等三位数字。在 CLI 输入即可直达，无需逐级进菜单。 |
| **interface 标签** | 沿表面法向的两个金属-水界面：`normal_aligned`（+轴指向）/ `normal_opposed`（-轴指向）。 |
| **帧目录命名** | Bader 用 `bader_t<step>_i<frame>/`，分布式 SP 电势用 `potential_t<step>_i<frame>/`。按 `_t(\d+)` 数值排序。 |
| **单位** | 距离 Å，能量 eV（内部 Hartree 自动转换），分数坐标 [0,1)，时间 fs，TI/SG 内部用 a.u. |
| **cSHE 参考** | 计算氢电极。U vs SHE = -E_Fermi + φ_center + ΔΨ - μ(H⁺) - ΔE_ZP。 |
| **TI 约束力** | CP2K LagrangeMultLog 文件的 SHAKE 乘子 λ；TI 中 dA/dξ = -⟨λ⟩。 |

---

## 非交互场景

CLI 是给人交互用的。如果你想在脚本里调用、批量跑多体系、或者从 Jupyter 调用，应该用：

- **`md_analysis.workflows` 的 `run_*()`** — 程序化入口，跟 CLI 菜单一一对应；返回 `WorkflowResult`（`md_analysis.main` 是同一批名字的薄 re-export facade）
- **`md_analysis.agent.dispatch(task, params)`** — agent-friendly 统一调度入口（可序列化、JSON Schema）

详见 [workflows.md 末尾"非交互场景出口"](workflows.md#非交互场景出口)。

---

## 反馈

文档错漏 / CLI 行为跟文档不一致 → GitHub issue（仓库 `StackOverFlow11/MD-analysis`）。
