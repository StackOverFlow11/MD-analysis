# 待修正问题说明（API 语义 / 单位一致性）

## 目的

本文用于明确当前项目中两个需要优先修正的问题，并给出可实践、可验收的修改思路，确保分析结果在**可复用性**与**科学表达一致性**上更可靠。

---

## 问题 1：API 语义与行为存在潜在偏差风险

### 现状描述

在 `scripts/structure/utils/LayerParser.py` 中，`detect_interface_layers()` 暴露了参数 `n_interface_layers`，但当前实现行为与该参数语义并不完全一致（实现更偏向固定策略，而非按参数控制的可变策略）。

### 风险

- 调用方可能误以为 `n_interface_layers` 会严格控制接口层数量，导致使用预期与真实行为不一致。
- 后续二次开发者很难从函数签名直接判断“哪些参数真实生效、哪些参数仅保留兼容”。
- 测试通过并不代表语义正确，容易产生“隐性回归”。

### 可实践修改思路

#### 方案 A（推荐，低风险渐进）

保留参数以兼容历史调用，但明确“当前版本策略”并增加迁移提示：

1. 在 `detect_interface_layers()` 的 docstring 中明确说明：
   - 当前默认策略是什么；
   - `n_interface_layers` 的实际作用范围；
   - 未来版本计划（若有）。
2. 当传入“当前不支持的语义值”时，给出显式 `warnings.warn(...)`，提示调用方行为可能与参数名直觉不一致。
3. 在 `scripts/structure/utils/__init__.py` 与相关分析层 docstring 中同步更新说明，避免入口文档与底层实现描述不一致。
4. 增加测试覆盖：
   - 参数取值与输出行为关系测试；
   - warning 触发测试（若采用 warning 方案）。

#### 方案 B（语义彻底统一，改动更大）

让 `n_interface_layers` 真正驱动结果：

1. 明确定义“每侧取 N 层”或“总计取 N 层”的规则（必须二选一，写进文档）。
2. 修改接口层选择逻辑，使结果严格遵循规则。
3. 对外发布时标注行为变更，并补充迁移说明。

### 验收标准

- 函数签名、文档、实际行为三者一致。
- 调用方不再需要“读源码猜语义”。
- 测试可稳定验证参数语义，不依赖隐含实现细节。

---

## 问题 2：图示单位与数据单位存在不一致风险

### 现状描述

在 `scripts/structure/Analysis/Water.py` 中，图中 y 轴单位文本与 CSV 数据列含义存在不一致风险（例如密度单位书写格式与常规定义、取向加权量纲表达与数据列说明不一致）。

### 风险

- 图与数据不能直接互证，影响报告可信度。
- 外部读者可能基于图注做错误物理解释。
- 后续论文/汇报复用图件时，单位错误会被放大。

### 可实践修改思路

1. 先建立“单位唯一事实源”（建议写在分析层文档中），明确每个输出量的物理量纲：
   - 质量密度：`g/cm^3`（或规范记法 `g·cm^-3`）；
   - 取向加权密度：明确是 `1/Å^3` 还是无量纲缩放量（需二选一）。
2. 统一修改以下位置，确保文本一致：
   - `scripts/structure/Analysis/Water.py` 的坐标轴标签；
   - `scripts/structure/Analysis/WaterAnalysis/*.py` 的 CSV header；
   - 相关 README / 使用文档中的单位说明。
3. 若第二幅图保留 `a.u.`，必须在文档中写清“a.u. 的构成与归一化方式”；否则建议直接改为物理单位，避免歧义。
4. 增加回归测试：
   - 检查关键 CSV header 是否包含约定单位；
   - 检查绘图标签文本是否与约定一致。

### 验收标准

- 同一物理量在图、CSV、文档中的单位表达完全一致。
- 任意用户仅通过输出文件即可正确解释量纲，不依赖口头说明。

---

## 建议执行顺序（可直接落地）

1. 先修问题 2（单位一致性）：改动小、见效快、风险低。
2. 再修问题 1（API 语义）：先做方案 A（文档+warning+测试），稳定后再评估是否推进方案 B。
3. 每一步都附带测试更新，避免“修了说明、没锁住行为”。

---

## 建议改动文件清单

- `scripts/structure/utils/LayerParser.py`
- `scripts/structure/utils/__init__.py`
- `scripts/structure/Analysis/Water.py`
- `scripts/structure/Analysis/WaterAnalysis/WaterDensity.py`
- `scripts/structure/Analysis/WaterAnalysis/WaterOrientation.py`
- `scripts/structure/Analysis/WaterAnalysis/AdWaterOrientation.py`
- 相关测试文件（`test/integration/structure/Analysis/...`）

