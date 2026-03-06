# 2026-03-06 Claude Code Prompt — Rename Charge Outputs to `aligned/opposed`

你现在在 `MD-analysis` 项目中工作。请围绕 `src/md_analysis/charge/BaderAnalysis.py` 完成一次**明确的命名迁移**：

- 项目已正式决定：`charge` 模块中两列表面电荷输出的语义，不再使用 `bottom/top`
- 统一改为：`aligned/opposed`
- 这里的 `aligned/opposed` 对应 `LayerParser` 中稳定的界面标签：
  - `normal_aligned`
  - `normal_opposed`

## 背景

当前 `development` 分支已经修复了一个关键 bug：`BaderAnalysis` 不再按 `center_frac` 排序输出，而是按稳定的界面标签输出，因此两列的真实语义已经是：

- index `0` = aligned
- index `1` = opposed

但是项目里仍然残留了大量旧命名：

- `sigma_bottom`
- `sigma_top`
- `n_bottom`
- `n_top`
- `charge_per_surface_e` 中的 bottom/top 注释
- CSV 列名中的 `sigma_bottom_*` / `sigma_top_*`
- CLI 输出中的 bottom/top
- `context4agent` 文档中的 bottom/top
- 可能还有 README、测试断言、图例、注释、docstring

因此需要一次**完整、系统、语义一致**的重命名与文档对齐。

## 目标

请完成以下工作：

1. 将 `charge` 模块对外可见的两列语义统一改名为 `aligned/opposed`
2. 确保代码实现、docstring、CSV 列名、CLI 输出、图例、测试、`context4agent` 文档全部一致
3. 不要只改文档，也不要只改实现；必须做完整对齐
4. 保持 `PBC shift` 不换列这一语义稳定性

## 需要重点检查和修改的范围

请至少检查并按需修改以下位置：

### 代码

- `src/md_analysis/charge/BaderAnalysis.py`
- `src/md_analysis/charge/__init__.py`
- `src/md_analysis/main.py`
- `src/md_analysis/CLI.py`

### 测试

- `test/unit/charge/test_charge_analysis.py`
- `test/integration/charge/test_charge_trajectory.py`
- 任何因 CSV 列名或 CLI 文案变化而需要同步的测试

### 文档与契约

- `context4agent/architecture/modules/src/charge/interface_exposure.md`
- `context4agent/architecture/modules/src/charge/implementation_guidelines.md`
- `context4agent/architecture/modules/data_contract.md`
- `context4agent/architecture/modules/glossary_units.md`（如果术语定义需要同步）
- `README.md`（如果对外示例或说明出现 bottom/top）
- 任何已经公开承诺 `bottom/top` 的说明文字

## 实现要求

### 1. 统一命名语义

以下命名应系统迁移到 `aligned/opposed` 风格：

- `sigma_bottom` -> `sigma_aligned`
- `sigma_top` -> `sigma_opposed`
- `sigma_bottom_cumavg_*` -> `sigma_aligned_cumavg_*`
- `sigma_top_cumavg_*` -> `sigma_opposed_cumavg_*`
- `n_bottom` / `n_top` -> `n_aligned` / `n_opposed`
- `q_bottom` / `q_top` -> `q_aligned` / `q_opposed`

如果某些内部临时变量仍保留 `bottom/top`，但其真实语义已经是 `aligned/opposed`，也请一并改掉，避免后续误导。

### 2. 对外输出必须同步改名

请特别注意这些地方：

- `atoms.info["surface_charge_density_e_A2"]`
  - 值顺序仍是长度为 2 的数组，但文档语义应改为 `[σ_aligned, σ_opposed]`
- `atoms.info["surface_charge_density_uC_cm2"]`
  - 同样改为 `[σ_aligned, σ_opposed]`
- `atoms.info["n_charged_atoms_per_surface"]`
  - 文档语义改为 `[n_aligned, n_opposed]`
- `atoms.info["charge_per_surface_e"]`
  - 文档语义改为 `[Σq_aligned, Σq_opposed]`

### 3. CSV 与图表

请将 `surface_charge_analysis()` 输出中的 CSV 字段从 `bottom/top` 改为 `aligned/opposed`，并同步：

- 写 CSV 的字段名
- 读 CSV 的 CLI 汇总代码
- PNG 图例
- verbose 打印文案
- 所有与这些列名绑定的测试

### 4. CLI 和用户可见文案

如果 CLI 仍打印：

- `sigma_bottom`
- `sigma_top`
- `charged atoms near bot`
- `charged atoms near top`

请全部替换成语义一致的 `aligned/opposed`。

### 5. 文档同步是硬要求

本项目有明确要求：**代码变更必须同步 `context4agent`**。  
请不要遗漏：

- 接口暴露文档
- 实现准则文档
- 数据契约文档
- 如果术语变化触及“术语/单位说明”，也要同步 glossary

## 约束

- 不要修改与本次命名迁移无关的功能逻辑
- 不要顺手修复其他无关问题
- 不要引入新的语义混用，比如文档写 `aligned/opposed`，代码还叫 `bottom/top`
- 不要保留“新旧术语混搭”的中间状态
- 若必须保留兼容层，请明确说明兼容策略和理由；否则默认直接全量切换

## 验收标准

完成后请确认：

1. `charge` 相关代码中对外两列语义已统一为 `aligned/opposed`
2. 不再存在公开输出字段或用户可见文案中的 `bottom/top`
3. `context4agent` 中相关文档已同步
4. `charge` 相关测试全部通过
5. 如有必要，新增/更新测试以验证：
   - `layer` 方法 `PBC shift` 不换列
   - `counterion` 方法输出列也遵循 `aligned/opposed`

## 建议执行步骤

建议按以下顺序操作：

1. 全局搜索 `sigma_bottom|sigma_top|bottom|top`，限定在 `charge` 相关代码与文档
2. 先统一 `BaderAnalysis.py` 的变量名、CSV 字段、绘图标签、verbose 输出
3. 再同步 `CLI.py` / `main.py`
4. 再同步 `context4agent`
5. 最后更新测试并运行：
   - `pytest test/unit/charge/test_charge_analysis.py test/integration/charge/test_charge_trajectory.py -q`

## 最终输出要求

完成后请给出：

1. 修改摘要
2. 受影响的公开命名/CSV 列名变化
3. 运行了哪些测试
4. 是否存在潜在兼容性影响
