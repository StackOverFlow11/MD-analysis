# `md_analysis` 内部实现准则（当前实现口径）

> 适用范围：`src/md_analysis/` 顶层包（当前核心为 `src/md_analysis/__init__.py`）。
>
> 目标：保持一个稳定、轻量、低耦合的顶层入口，不在该层承载业务算法。

## 1. 角色定位与边界

- `md_analysis` 是**顶层命名空间与聚合入口**，不是业务计算层。
- 顶层负责"把子包暴露出来"，不负责"执行分析逻辑"。
- 业务实现下沉到三个子包：
  - `md_analysis.utils`：单帧底层工具
  - `md_analysis.water`：水分析多帧工作流
  - `md_analysis.potential`：电势分析多帧工作流
- 编程入口 `main.py` 和 CLI 入口 `CLI.py` 位于顶层，负责协调子包调用。
- 目录文档治理采用硬约束：`context4agent/architecture/modules/src/` 必须镜像 `src/md_analysis/`，且每个子目录都需维护
  `interface_exposure.md` 与 `implementation_guidelines.md`。

## 2. 导出稳定性等级

- 顶层导出按"稳定接口"管理：
  - 已进入 `src/md_analysis/__init__.py` 且在 `__all__` 中声明的符号，视为对外契约。
  - 未进入 `__all__` 的符号，不承诺稳定性。
- 当前顶层仅暴露三个子包：`utils`、`water`、`potential`。

## 3. `__init__.py` 实现准则

- `src/md_analysis/__init__.py` 只允许：
  - 导出聚合（`from . import xxx`）
  - `__all__` 显式声明
  - 最小注释/模块说明
- `src/md_analysis/__init__.py` 禁止：
  - 数值计算、文件 I/O、图像生成、CLI 逻辑
  - 运行时副作用（例如导入即执行分析）
  - 重型三方库导入（如 `numpy`、`matplotlib`、`ase`）

## 4. 依赖方向约束

- 允许方向：`md_analysis` -> `md_analysis.utils` / `md_analysis.water` / `md_analysis.potential`
- 允许方向：`md_analysis.water` -> `md_analysis.utils`
- 允许方向：`md_analysis.potential` -> `md_analysis.utils`
- 禁止反向依赖：子包不应依赖 `md_analysis` 顶层内部状态
- 禁止跨层耦合：
  - `md_analysis.water` 与 `md_analysis.potential` 之间不互相依赖
  - 顶层不得感知具体实现文件（如 `WaterParser.py`）

## 5. 导出变更规则

当新增/删除/重命名顶层导出时，必须同步完成：

1. 更新 `src/md_analysis/__init__.py` 中的导入与 `__all__`
2. 更新接口文档 `context4agent/architecture/modules/src/interface_exposure.md`
3. 运行导入烟雾测试（至少覆盖 `import md_analysis` 与 `from md_analysis import ...`）
4. 运行回归测试，确认无破坏性影响

## 6. 向后兼容策略

- 默认策略：**不破坏现有导入路径**。
- 若必须调整顶层导出：
  - 优先采用"先兼容、后迁移"的两阶段方式
  - 在迁移窗口内保留旧导入路径的可用性
- 只有在明确需要时，才进行破坏性移除。

## 7. 反模式清单（禁止）

- 在顶层 `__init__.py` 写业务函数
- 为了"省导入"把下层所有符号全部扁平导出到顶层
- 未更新 `__all__` 就修改导出行为
- 未更新文档与测试就变更公开接口

## 8. 维护检查清单（提交前）

- [ ] 顶层导出是否仍保持最小且清晰
- [ ] 是否引入了不必要的新依赖
- [ ] `__all__` 与实际导出是否一致
- [ ] 导入烟雾测试是否通过
- [ ] 对应接口文档是否同步
