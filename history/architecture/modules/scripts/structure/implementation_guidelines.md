# `scripts/structure/` 内部实现准则（当前实现口径）

> 适用范围：`scripts/structure/` 层（当前核心为 `scripts/structure/__init__.py`）。
>
> 目标：作为"结构分析"统一入口层，向外提供稳定 API，向内屏蔽实现细节。

## 1. 角色定位与边界

- `scripts/structure/` 是**领域入口层（facade）**：
  - 面向调用方提供稳定导入路径
  - 向下聚合 `scripts/structure/utils/` 与 `scripts/structure/Analysis/` 子包能力
- 本层不应承载复杂算法实现，算法必须留在 `utils` 子层。
- 面向"多帧统计 + 文件输出"的分析流程实现应放在 `Analysis` 子包，不应放在 facade 的 `__init__.py`。

## 2. 导出分层模型

当前导出建议分三类维护：

1. **核心数据结构与异常**
   - 如 `Layer`、`SurfaceDetectionResult`、`SurfaceGeometryError`、`WaterTopologyError`
2. **高频业务函数**
   - 如层识别、水分子标记、密度/取向统计函数
3. **默认参数常量**
   - 如 `DEFAULT_Z_BIN_WIDTH_A`、`DEFAULT_THETA_BIN_DEG`

要求：

- 仅导出"跨调用方复用"的稳定符号
- 仅内部使用的工具函数不得上浮到本层

## 3. `scripts/structure/__init__.py` 实现准则

- 必须显式维护 `__all__`，使公开 API 可审计
- 导出命名应与实现语义一致，避免别名漂移
- 禁止在 `__init__.py` 内：
  - 写统计计算逻辑
  - 触发 I/O 或绘图副作用
  - 引入与入口聚合无关的依赖

## 4. 依赖方向约束

- 允许方向：`scripts.structure` -> `scripts.structure.utils`
- 允许方向：`scripts.structure.Analysis` -> `scripts.structure.utils`
- 禁止方向：`scripts.structure.utils` 反向依赖 `scripts.structure` 入口
- 目标：避免导入环，保持"入口层薄、实现层厚"

## 5. 命名与语义一致性

- 分布类接口命名必须体现统计方向/对象（如 `*_z_distribution`）
- 角度/窗口类接口命名必须体现统计域（如 `*_theta_pdf_*`）
- 相同物理量在不同函数中的命名含义必须保持一致（例如 `oxygen_indices`、`c_fraction_range`）

## 6. 契约同步要求

凡是以下变化，均视为"契约变化"：

- 输入输出 shape 改变
- 单位或物理口径改变
- 默认参数改变
- 异常行为（抛错条件/类型）改变

发生契约变化时，必须同步更新：

- `history/architecture/modules/data_contract.md`
- `history/architecture/modules/glossary_units.md`（若涉及术语/单位）
- 对应 `interface_exposure.md`

## 7. 兼容性策略

- 默认保持导入路径稳定（`from scripts.structure import ...`）
- 若需替换接口：
  - 优先新增而非直接覆盖
  - 先保留旧接口并提供迁移窗口
- 破坏性移除前，必须明确迁移目标与时间点

## 8. 测试与验证要求

每次入口层导出变更，至少完成：

1. 导入烟雾测试：
   - `import scripts.structure`
   - 关键导出 `hasattr(...)` 检查
2. 相关单元/集成测试回归
3. 文档一致性检查（导出列表与文档一致）

## 9. 反模式清单（禁止）

- 把实现层私有 helper 直接导出到入口层
- 在未更新 `__all__` 的情况下新增导出
- 接口改名但未保留兼容路径
- 文档与测试未同步即合并导出改动

## 10. 维护检查清单（提交前）

- [ ] 导出符号是否仍是"面向外部调用"的稳定集合
- [ ] `__all__` 是否与实际导出一致
- [ ] 是否引入了不必要的跨层耦合
- [ ] 契约文档是否同步更新
- [ ] 导入烟雾测试与回归测试是否通过
