# `scripts/` 内部实现准则（当前实现口径）

> 适用范围：`scripts/` 顶层包（当前核心为 `scripts/__init__.py`）。
>
> 目标：保持一个稳定、轻量、低耦合的顶层入口，不在该层承载业务算法。

## 1. 角色定位与边界

- `scripts/` 是**顶层命名空间与聚合入口**，不是业务计算层。
- 顶层负责“把子包暴露出来”，不负责“执行分析逻辑”。
- 任何与结构识别、水统计、分布计算相关的实现，必须下沉到 `scripts/structure/` 及其子层。
- 目录文档治理采用硬约束：`history/architecture/modules/scripts/` 必须镜像 `scripts/`，且每个子目录都需维护
  `interface_exposure.md` 与 `implementation_guidelines.md`。

## 2. 导出稳定性等级

- 顶层导出按“稳定接口”管理：
  - 已进入 `scripts/__init__.py` 且在 `__all__` 中声明的符号，视为对外契约。
  - 未进入 `__all__` 的符号，不承诺稳定性。
- 当前建议顶层仅暴露子包（如 `structure`），避免暴露细粒度函数。

## 3. `__init__.py` 实现准则

- `scripts/__init__.py` 只允许：
  - 导出聚合（`from . import xxx`）
  - `__all__` 显式声明
  - 最小注释/模块说明
- `scripts/__init__.py` 禁止：
  - 数值计算、文件 I/O、图像生成、CLI 逻辑
  - 运行时副作用（例如导入即执行分析）
  - 重型三方库导入（如 `numpy`、`matplotlib`、`ase`）

## 4. 依赖方向约束

- 允许方向：`scripts` -> `scripts.structure`
- 禁止反向依赖：`scripts.structure` 不应依赖 `scripts` 顶层内部状态
- 禁止跨层耦合：
  - `scripts` 不直接依赖 `test/`、`history/`、`data_example/`
  - 顶层不得感知具体实现文件（如 `WaterParser.py`）

## 5. 导出变更规则

当新增/删除/重命名顶层导出时，必须同步完成：

1. 更新 `scripts/__init__.py` 中的导入与 `__all__`
2. 更新接口文档 `history/architecture/modules/scripts/interface_exposure.md`
3. 运行导入烟雾测试（至少覆盖 `import scripts` 与 `from scripts import ...`）
4. 运行回归测试，确认无破坏性影响

## 6. 向后兼容策略

- 默认策略：**不破坏现有导入路径**。
- 若必须调整顶层导出：
  - 优先采用“先兼容、后迁移”的两阶段方式
  - 在迁移窗口内保留旧导入路径的可用性
- 只有在明确需要时，才进行破坏性移除。

## 7. 反模式清单（禁止）

- 在顶层 `__init__.py` 写业务函数
- 为了“省导入”把下层所有符号全部扁平导出到顶层
- 未更新 `__all__` 就修改导出行为
- 未更新文档与测试就变更公开接口

## 8. 维护检查清单（提交前）

- [ ] 顶层导出是否仍保持最小且清晰
- [ ] 是否引入了不必要的新依赖
- [ ] `__all__` 与实际导出是否一致
- [ ] 导入烟雾测试是否通过
- [ ] 对应接口文档是否同步
