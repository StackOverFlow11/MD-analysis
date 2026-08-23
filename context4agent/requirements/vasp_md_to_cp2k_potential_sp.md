# 需求：md-analysis 增加「VASP MD 工作目录 → CP2K 单点目录生成」（→ 端到端 PZC 工具包）

> 记录日期 2026-08-10，触发场景：同门 sjy 的 Pt(111)/水 PZC 计算（chem `~/temp/bader-5000md1-3/`）。
> 原始草稿 + 案例分析在 03 项目 `/home/shaofl/Projects/research_projects/03_CuAg_1_SAA_CO2RR/tmp/sjy/`（`pzc_estimate.md` / `plot_phi_z.py` 原型 / 50 帧 PLANAR_AVERAGE 数据）。

## 产品定位（2026-08-10 user 补充）

**目标用户 = 只熟悉 VASP、不懂 CP2K、也不熟悉程俊组电势方法（Le 2017 cSHE alignment）的人。**
用户（本仓库主人）要用 md-analysis 给同门封装一个**开箱即用的工具包**：给 VASP MD 目录，出 PZC vs SHE，中间不要求任何 CP2K / 方法学知识。

这决定了设计取向：

- **模板机制（SpGen 范式）**：内置经过验证的 CP2K 电势 SP 模板做兜底（PBE + GTH + DZVP-MOLOPT、金属 smearing、`&V_HARTREE_CUBE`），同时支持显式传模板 / config 配置——维护者调好模板丢给同门即可，用户不需要会写 CP2K inp；元素 → 基组/赝势的映射工具内做（GTH_POTENTIALS + BASIS_MOLOPT_UZH 已确认覆盖 Pt/O/H/Cu/Ag 等本项目常用元素）
- **方法学内置**：Le 2017 中心方程 + 常数（15.35 / 15.81 / 0.35 eV）直接封装进分析端，用户拿到的是 U_PZC 数字和图，不是一堆 cube
- **护栏代替知识**：平台窗口自动检测 + 平台质量告警、赝势约定检查、正交 cell 校验——把"sjy 案例踩过的坑"全部变成自动检查/报错

## 背景与动机

sjy 手头是 **VASP MD** 轨迹（Pt(111)/水，5000 步，每 100 fs 抽帧）。要用 Le 2017（PRL 119, 016801）cSHE alignment 算 PZC，公式常数 Δ_DP A^(w)=15.35 eV **绑定 CP2K-GTH 赝势约定**（HPS 依赖赝势，原文明确警告；VASP-PAW 实测偏差 ~1.2 V，详 03 项目 `tmp/sjy/pzc_estimate.md`）。

因此对 VASP MD 的帧，需要跑 **CP2K 单点**（PBE + GTH + MOLOPT，打 `V_HARTREE` cube + smearing 出 E_F），让 φ_wat/E_F 落在 CP2K 约定里，常数才能直接用。同理也覆盖一般性需求：VASP 轨迹帧 → CP2K 电势/普通 SP。

## 现状缺口（2026-08-10 核实源码）

`md_analysis.scripts.PotentialGen`（CLI 431/432）本身已经是「轨迹帧 → CP2K SP 目录」，但**输入只接 CP2K xyz 轨迹**：

- `generate_potential_workdir(atoms, ...)` / batch 版吃的是 CP2K `md-pos-1.xyz` + CP2K `sp.inp` 模板 + CP2K 提交脚本
- 没有「VASP MD 工作目录」这个输入入口：XDATCAR / vasprun.xml 的解析、POSCAR 元素顺序、VASP 晶胞（可能非正交）都没接
- `BaderGen` 有现成的「帧 → VASP 目录」和 `IndexMapper`（CP2K↔VASP 原子序映射），但方向是 VASP，不是 CP2K

**缺口一句话**：帧源层只认 CP2K xyz，缺 VASP MD 目录 → ase 帧序列 → 现有 CP2K 生成器的桥。

## 需求规格

### A. VASP MD 帧源（新）

- 输入：VASP MD 工作目录，解析优先级 `vasprun.xml` > `XDATCAR`（ase.io 都能读；XDATCAR 需配 POSCAR 定元素顺序）
- 抽帧参数复用现有 `FrameSelection`（start/end/step）
- 晶胞校验：md-analysis 边界是**正交晶胞**；VASP 非正交晶胞要显式报错或警告，别静默写出错 cell 的 sp.inp
- 元素列表从 POSCAR/vasprun 读，校验 CP2K sp.inp 模板的 `&KIND` 覆盖所有元素（缺 KIND 直接报清哪个元素缺）

### B. CP2K SP 目录生成（复用 + 小改）

- 主流程直接复用 `PotentialGen`：帧 → xyz 坐标文件 + 模板替换 + script.sh
- **模板三级解析（照 SpGen 范式，保灵活度）**：
  1. 显式参数 `inp_template_path`（CLI/workflows 传入）——**优先**
  2. 持久化 config（Settings 设置，参考 `KEY_DP_SP_INP_TEMPLATE_PATH` 的做法新增一个 potential 模板 key）
  3. 包内置默认电势模板（PBE + GTH + DZVP-MOLOPT + smearing + `&V_HARTREE_CUBE`，随包发行）——兜底，零知识开箱即用
- 这样维护者（本仓库主人）只需把调好的模板文件交给同门 / 或写进 config，同门不用懂模板内容；模板本身完全可改（换泛函、换基组、加 print 项都自由）
- 模板硬要求（写进校验或文档）：`&PRINT &V_HARTREE_CUBE`、`&SCF &SMEAR`（金属出 Fermi energy）
- 入口：CLI 菜单编号（如 433/434）+ `workflows` facade（如 `run_potential_gen_from_vasp_md`）；agent dispatch 入口按 `development_agent` 分支策略后接

### C. 端到端封装（产品定位要求的核心交付）

一条命令完成：**VASP MD 目录 →（抽帧 + CP2K SP 目录生成）→ 用户提交计算 →（收集 cube/out → φ(z) 图 + E_F 累积平均图 + 平台自动检测 → U_PZC vs SHE）**

- 下游分析复用现有 `potential_full` 分布式模式（cube + sp.out Fermi），只需确认对 VASP 来源帧端到端通
- **Le 2017 方程内置**：`U_PZC = (−E_F − eφ_wat − 0.81)/e`，常数固化并注明出处与适用条件（CP2K-GTH only）
- 平台自动检测 + 质量告警：宏观平均后中段漂移 > ~0.2 eV 时明确警告"水膜可能无 bulk 平台"（sjy 案例的教训）
- 输出物：φ(z) 叠加图、E_F 累积平均图、U_PZC 数值 + 不确定度声明

## ⚠️ 设计约束

- 这个功能存在的理由就是赝势约定：文档写死「Le 2017 常数只能配 CP2K-GTH 的 φ/E_F」；VASP LOCPOT 路线（LVHAR）如需支持另开需求，且必须走 $V_0$ 水-真空自洽公式
- 正交晶胞约束沿用包级边界条件，不做任意法向量

## 验收标准（建议）

- 给定 sjy 的 VASP MD 目录，一条 CLI / workflows 调用生成 50 帧 CP2K SP 目录（xyz + sp.inp + script.sh），sp.inp 的 ABC 与 POSCAR 晶胞一致、KIND 覆盖 Pt/O/H
- 跑完后现有 `potential_full` 分布式模式直接出图出数，无需手工整理文件名
- 同帧 CP2K vs VASP 的 (φ_wat, E_F) 对比可**直接测出 PAW↔GTH 零点偏移 δ**（预期 ~1.2 V），作为正确性检查
- 单元测试：XDATCAR 帧源解析（帧数/元素顺序/cell）、非正交 cell 报错、KIND 缺失报错、平台质量告警触发
