# 主线工作流（指针）

[← 回索引](README.md)

按"想做什么"组织。每条工作流详细页给：
- **预期输入目录结构**（你拿到/产出的数据该长什么样）
- **CP2K inp / VASP INCAR 要求打印的内容**（哪些 section / keyword 必须开）
- **内部数据流示意**（文件 → 解析模块 → dataclass → 分析函数 → 输出）
- **执行步骤**（CLI 终端会话）
- **输出位置 + 常见坑**

---

## 五条主线

| # | 工作流 | 详细页 |
|---|---|---|
| 01 | 算电极电势 U vs SHE（cSHE） | [workflows/01_potential_cshe.md](workflows/01_potential_cshe.md) |
| 02 | 表面电荷 σ + φ 标定 | [workflows/02_surface_charge.md](workflows/02_surface_charge.md) |
| 03 | Slow-Growth 粗扫 → TI 精算 | [workflows/03_sg_to_ti.md](workflows/03_sg_to_ti.md) |
| 04 | Bader 批量生成 VASP 工作目录 | [workflows/04_bader_batch.md](workflows/04_bader_batch.md) |
| 05 | SP 单点为 DeePMD 训练集抽帧 | [workflows/05_sp_for_dp.md](workflows/05_sp_for_dp.md) |

---

## 非交互场景出口

CLI 是给人看的。如果你要在脚本 / Jupyter / 集群作业里跑：

### 程序化入口（`main.py`）

```python
from md_analysis.main import (
    run_water_analysis, run_potential_analysis, run_charge_analysis,
    run_tracked_charge_analysis, run_counterion_charge_analysis, run_all,
)

# 跟 CLI 菜单一一对应；output_dir 是最终写入目录（不再前置 water/ 等）
run_water_analysis(
    xyz_path="md-pos-1.xyz",
    cell_abc=(10.22, 10.22, 26.42),
    output_dir="output/water/",
)
```

### Agent 入口（`agent.dispatch`）

JSON 序列化、JSON Schema 自描述，方便给 LLM / MCP / 批处理用：

```python
from md_analysis.agent import dispatch, list_tasks, get_task_schema

for t in list_tasks():
    print(t["name"], "—", t["description"])

schema = get_task_schema("ti_full_analysis")

result = dispatch("ti_full_analysis", {
    "root_dir": "./ti_runs/",
    "output_dir": "./output/",
    "epsilon_tol_ev": 0.01,
    "equilibration": 500,
})
print(result.success, result.summary["delta_A_eV"])
```

任务列表（14 个）见 [menu_reference.md 末尾](menu_reference.md#agent-任务列表)。

---

[← 回索引](README.md)
