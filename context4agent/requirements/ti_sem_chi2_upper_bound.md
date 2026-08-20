# 功能改进:约束 TI 的 SEM 卡方上界(替代 N_eff ≥ 50 硬门槛)

> 立项日期:2026-08-20。提出背景:group Ag1Cu4 disp 2k Volmer TI
> `ti_target_0.585260` 实测 —— 该点 SEM 贡献对 ΔA 仅 ~0.008 eV(远小于
> ε_tol=0.05 eV),却因 N_eff=16.5 < 50 被判 fail。硬门槛把"精度已富裕"的
> 点误杀,且 N_eff 判据与积分精度判据(sem_max)职能重叠。

## 目标

把"SEM 估计本身的可信度"从**硬门槛**(N_eff ≥ 50)改为**连续惩罚**:
用卡方分布给出 SEM 的单侧 95% 上界作为报告值,膨胀后的 SEM 同时用于
pass/fail 判定和 σ_A 积分传播。

## 设计(用户已拍板的四个决策,2026-08-20)

1. **公式**:
   $$\mathrm{SEM}_{\mathrm{report}} = \mathrm{SEM} \cdot \sqrt{\nu\,/\,\chi^2_{0.05,\nu}}$$
   即真实 SEM 有 95% 概率不超过报告值。置信水平**单侧 95% 硬编码为常量**
   (如 `DEFAULT_SEM_CONFIDENCE = 0.95`),不做用户可调参数。
2. **自由度 ν 取 block 平台处的块数**(ν = n_blocks_plateau − 1),
   **不用** ACF 的 N_eff —— sem_final 来自 F&P block average,其估计误差
   由平台处独立块数决定。block_average 需在平台判定时额外返回
   `n_blocks_plateau`。
3. **N_eff ≥ 50 硬门槛移除**,保留低地板:N_eff < 10 时直接 fail
   (τ_corr 估计本身失效,膨胀无意义)。地板值沿用
   `DEFAULT_GEWEKE_MIN_NEFF_SUBSERIES = 10` 量级,新增独立常量。
4. **膨胀后 SEM 进入 σ_A 传播**:`_integrate_free_energy` 的 `sems`
   改用膨胀值,使总自由能误差反映 SEM 的估计不确定度。

## 行为变化分类

**approved tightening(用户 2026-08-20 拍板)**:会改变旧分析结论
(原本 N_eff fail 的点在精度富裕时变 pass)。实施时:

- commit message 标 `BEHAVIOR CHANGE`;
- `ConstraintPointReport` 需新增字段(如 `sem_inflated`、`sem_inflation_factor`、
  `n_blocks_plateau`),CSV 列相应增加(同步 `data_contract.md`);
- 旧字段 `sem_final` 语义保持不变(= 未膨胀的 plateau SEM),避免静默破坏
  下游;判定与传播统一用膨胀值;
- failure_reasons 文案更新:SEM 超限时报告膨胀后数值与膨胀因子。

## 影响面

- `analysis/block_average.py`:返回平台块数
- `workflow.py`:`analyze_single_point` 的 pass 判定(4 条件中
  `passed_neff` 移除、`sem_ok` 用膨胀值)、`_auto_equilibrate` 同步
- `integration.py`:σ_A 用膨胀 SEM
- `plot.py` / CSV:诊断图与报告体现膨胀后 SEM
- 测试:卡方因子数值(ν=15→×1.44 量级)、边界(ν 极小/极大)、
  低地板 fail 路径、σ_A 传播一致性
- 文档:`context4agent/architecture/modules/` 下 TI 相关镜像文档同步

## 明确的非目标(不做)

- 不做 Student-t 均值置信区间(用户明确只要求卡方约束 SEM 上界);
- 不解决多盆地/双稳态导致的系统性偏差(该类问题超出平稳序列假设,
  任何 SEM 修正都无效,需轨迹几何层面诊断)。
