# `scripts/` 接口暴露约定（当前实现）

> 对应代码：`scripts/__init__.py`
>
> 本文档定义 `scripts` 顶层包的公开接口边界与兼容规则。

## 1. 接口角色定义

- `scripts/` 是顶层命名空间入口，只负责暴露子包。
- `scripts/` 不直接暴露业务算法函数、数据结构或配置常量。
- 业务能力统一通过 `scripts.structure` 向下访问。
- 目录治理硬约束见：`history/architecture/modules/README.md`（镜像对齐 + 双文档）。

## 2. 当前公开接口清单

### 2.1 稳定导出（Stable）

- `structure`

与代码一致性基线：

- `scripts/__init__.py` 中已导入并写入 `__all__` 的符号，视为顶层稳定接口。

### 2.2 非公开符号（Non-public）

- 未在 `__all__` 中声明的任何名字，均不视为公开契约。
- 调用方不得依赖非公开符号名称或存在性。

## 3. 推荐导入方式

- `import scripts`
- `from scripts import structure`

不推荐：

- 从 `scripts` 顶层直接期待细粒度函数（应转到 `scripts.structure`）。

## 4. 稳定性级别与兼容承诺

### 4.1 顶层稳定性级别

- `scripts.structure`：**稳定（Stable）**
  - 正常迭代中默认保持导入路径不变。

### 4.2 兼容承诺

- 在无明确迁移说明前，不移除顶层稳定导出。
- 若确需调整，必须提供迁移路径与窗口说明。

## 5. 接口变更触发条件

以下任一情况视为“顶层接口变更”：

- 新增顶层导出符号
- 删除顶层导出符号
- 重命名顶层导出符号
- 修改 `__all__` 与实际导出的一致性

## 6. 顶层接口变更流程（必须执行）

1. 更新 `scripts/__init__.py` 的导出与 `__all__`
2. 更新本文档公开接口清单
3. 更新 `history/architecture/modules/scripts/implementation_guidelines.md`（若规则被触发）
4. 执行导入烟雾测试：
   - `import scripts`
   - `from scripts import <symbol>`
5. 运行回归测试确认无导入层回归

## 7. 反模式（禁止）

- 在顶层暴露实现层私有函数
- 把 `scripts` 当作业务函数大杂烩入口
- 修改导出但不更新 `__all__` 与文档
- 未做导入验证即提交顶层导出变更

## 8. 提交前检查清单

- [ ] 顶层导出是否仍保持最小集合
- [ ] `__all__` 是否与公开接口清单一致
- [ ] 兼容导入路径是否保持稳定
- [ ] 文档与测试是否同步通过
