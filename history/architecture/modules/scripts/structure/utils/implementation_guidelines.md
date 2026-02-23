# `scripts/structure/utils/` 内部实现准则（当前实现口径）

> 适用范围：`scripts/structure/utils/`（`config.py`、`LayerParser.py`、`WaterParser.py`）。
>
> 目标：在明确物理口径与输入输出契约的前提下，提供可复用、可测试、可维护的底层实现。

## 1. 职责分层与边界

### `config.py`

- 只承载默认参数与常量：
  - 元素集合
  - 分箱步长
  - 几何阈值
  - 单位换算相关常量
- 禁止写任何业务计算逻辑。

### `LayerParser.py`

- 负责金属层识别、界面层标记、法向符号判定。
- 负责层级数据结构与摘要输出。
- 不负责水分子拓扑识别、角度 PDF 统计。

### `WaterParser.py`

- 负责水分子标记、氧索引提取、质量密度与取向统计。
- 负责 c 轴窗口角度 PDF 统计。
- 不负责界面层判别与金属层聚类。

## 2. 坐标与几何口径（必须一致）

### 2.1 z/c 方向统计统一口径

- z 分箱基于分数坐标第三轴与 `|c|`：
  - `z = f_c * |c|`
- 取向参考方向使用晶胞 `c_unit = c / |c|`，不使用笛卡尔固定 z 轴替代。

### 2.2 周期性与区间规则

- 分数坐标统一使用 wrap 后值（`[0, 1)`）。
- c 轴窗口区间采用左闭右开：
  - `start < end`：`[start, end)`
  - `start > end`：跨边界窗口 `[start, 1) U [0, end)`
  - `start == end (mod 1)`：全区间

### 2.3 分箱边界规则

- z 分箱数：`nbins = ceil(Lz / dz)`
- 最后一段边界直接落在 `Lz`，避免累计误差
- bin 归属使用 `searchsorted(..., side="right") - 1` 再 clip，确保边界可控

## 3. `LayerParser.py` 实现准则

### 3.1 输入与默认值

- `metal_symbols=None` 时，必须回退到 `DEFAULT_METAL_SYMBOLS`
- `normal` 支持 `"a" / "b" / "c"` 与显式向量
- 界面层策略固定为每侧 1 层（共 2 层），不提供 `n_interface_layers` 配置参数

### 3.2 分层与界面判定

- 先按法向投影做 1D 聚类
- 界面层仅保留直接面向非金属环境的一层（low/high 两侧各一层）
- 对 `"a" / "b" / "c"` 轴采用分数坐标 MIC 判法向；无法判定时允许几何回退

### 3.3 数据结构约束

- `Layer` 与 `SurfaceDetectionResult` 保持 dataclass 不可变语义
- 对外返回应优先使用 tuple，避免可变容器泄漏

## 4. `WaterParser.py` 实现准则

### 4.1 水分子标记

- O-H 连通判据基于 MIC 距离与 `oh_cutoff_A`
- 分配策略为"按距离排序的贪心匹配"
- 每个 H 在最终分配中最多归属一个 O

### 4.2 质量密度 z 分布

- 一个 O 计作一个 H2O 分子
- 质量密度单位固定为 `g/cm^3`
- 体积来自 `A_xy * Δz`，其中 `A_xy = ||a × b||`

### 4.3 取向加权 z 分布

- 角度向量为 H-O-H 角平分线（O 为起点）
- `cos(theta)` 必须按 `dot(bisector, c_unit) / |bisector|`
- 结果单位固定为 `g/cm^3`（公式：$\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin,cm^3}}$）

### 4.4 c 窗口角度 PDF

- 输入窗口为分数坐标 c 区间
- 输出长度为 `180 / ndeg`
- `ndeg` 允许浮点；要求 `180 / ndeg` 为整数（容差校验）
- PDF 单位固定为 `degree^-1`
- 若窗口内无样本，返回全零数组（不抛错）

## 5. 参数管理准则

- 默认参数只能从 `config.py` 读取，不得在函数内部重复硬编码
- 新增默认参数时必须：
  1. 写入 `config.py`
  2. 在调用链中显式接入
  3. 更新 `data_contract.md` 对应条目

## 6. 异常与错误信息准则

### 6.1 异常分类

- 几何/界面类异常优先使用 `SurfaceGeometryError`
- 水拓扑/取向类异常优先使用 `WaterTopologyError`
- 参数合法性问题使用 `ValueError`

### 6.2 错误信息要求

- 报错信息必须包含最小定位信息（参数名、索引、符号名或条件）
- 禁止吞错；捕获异常时需保留原问题语义

## 7. 数值稳定性准则

- 所有归一化前必须检查向量范数是否为 0
- `cos(theta)` 在反三角函数前需 clip 到 `[-1, 1]`
- 浮点整除判定（如 `180/ndeg`）使用容差，不使用纯整数 `%`

## 8. 性能与可扩展性准则

- 优先 vectorized 处理批量数据，循环仅用于必要几何步骤
- 大轨迹统计应采用逐帧累积策略，避免一次性持有全部帧
- 输出结构优先使用 numpy 数组，保持与后续统计/绘图链路兼容

## 9. 文档与测试同步准则

发生以下改动时，必须同步测试与文档：

- 新增/修改公开函数
- 默认参数变化（如 `dz`、`ndeg`、阈值）
- 统计口径变化（单位、方向、分箱规则、窗口规则）

同步范围至少包括：

- `test/unit/structure/utils/`
- `test/integration/structure/utils/`
- `history/architecture/modules/data_contract.md`
- `history/architecture/modules/glossary_units.md`（若涉及术语/单位）

## 10. 提交前检查清单

- [ ] 模块职责是否越界（facade/实现/config 是否混杂）
- [ ] c 方向口径是否全链路一致
- [ ] 默认参数是否全部来自 `config.py`
- [ ] 异常类型与报错信息是否清晰
- [ ] 单元与集成测试是否通过
- [ ] 契约文档是否同步更新
