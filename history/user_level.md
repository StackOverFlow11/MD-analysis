# 用户目前水平（持续更新）

用途：让后续协作方/Agent 快速知道当前背景与偏好，减少反复确认。

## 背景与经验

- **研究方向/体系**：金属-水界面、抗衡离子；当前重点：恒电荷框架下的 CP2K AIMD/后处理
- **CP2K 熟悉度**：常用输入模块（DFT/SCF/MD/PRINT）、常见输出文件的含义与基本排错能力
- **Python 熟悉度**：
  - 不熟悉：numpy/pandas/matplotlib 基本用法、Python 语法特性
- **数据处理习惯**：偏好命令行工具；希望自动画图与导出 CSV

## 协作偏好

- **输出语言**：中文为主；代码注释/变量命名倾向英文
- **节奏偏好**：先跑通最小可用流程（MVP）→ 再扩展功能
- **可接受的依赖**：允许 ASE，技术栈主要采用 Python
- **conda 环境**：`env_md_an`（位于 `/home/shaofl/miniconda3/envs/env_md_an`，Python 3.12，含 numpy/matplotlib/ase/pytest）

## 当前已掌握的样例数据（来自仓库）

- `data_example/potential/`：包含 `md.inp`、`md.out`、`md-1.ener`、`md-pos-*.xyz`、`CHARGE.mulliken`、`md-POTENTIAL-*.cube` 等
