# modules 维护硬约束（强制）

> 适用范围：`history/architecture/modules/` 全目录。

## 1. 目录镜像硬约束

- `history/architecture/modules/scripts/` 的目录结构必须与代码目录 `scripts/` 保持一致。
- 对 `scripts/` 中每一个业务子目录（忽略 `__pycache__`、二进制缓存目录）：
  - 必须在 `history/architecture/modules/scripts/` 下存在同名镜像子目录。

## 2. 双文档硬约束

- 每一个镜像子目录都必须包含两份文档：
  - `interface_exposure.md`
  - `implementation_guidelines.md`
- 以上两份文档缺一不可，不允许只写其中一份。

## 3. 颗粒度硬约束

- 文档细致程度至少对齐 `scripts/structure/utils/` 的现有粒度。
- 最低要求：
  - **接口文档**：公开符号清单、导入方式、非公开边界、变更流程
  - **实现文档**：职责边界、依赖方向、口径约束、同步要求、检查清单

## 4. 变更同步硬约束

- 若 `scripts/` 目录发生新增/删除/重命名子目录：
  1. 同步调整 `history/architecture/modules/scripts/` 目录镜像
  2. 为受影响目录同步维护两份文档
  3. 更新 `history/architecture/README.md` 的模块索引

## 5. 例外范围

- 以下目录不纳入镜像硬约束：
  - `__pycache__/`
  - 运行期临时目录与缓存目录

## 6. 记录落位硬约束

- 全局性约束仅允许写入：
  - `history/architecture/modules/data_contract.md`
  - `history/architecture/modules/glossary_units.md`
- 除上述两份全局文档外，其他所有记录（模型上下文补充、目录级约定、实现细节）必须写入
  与项目架构镜像一致的：
  - `**/interface_exposure.md`
  - `**/implementation_guidelines.md`

## 7. 落位自检清单模板（每次变更后执行）

> 用法：复制下面模板到当次变更记录中，逐项勾选；未满足项需在提交前补齐。

```md
### 落位自检清单

- [ ] 本次是否涉及“全局契约/单位术语”变更？
  - 若是：仅更新 `data_contract.md` 与/或 `glossary_units.md`
- [ ] 本次是否涉及某个 `scripts/**` 子目录的接口或实现行为变化？
  - 若是：已同步更新该镜像目录下
    - `interface_exposure.md`
    - `implementation_guidelines.md`
- [ ] 本次是否新增/删除/重命名了 `scripts/**` 业务子目录？
  - 若是：`history/architecture/modules/scripts/**` 是否完成镜像同步
- [ ] 是否错误地把目录级细节写进了全局文档？
  - 若是：已迁移到对应镜像目录文档
- [ ] 目录粒度是否达到 `scripts/structure/utils/` 参考水平？
  - 至少包含：公开清单、边界、依赖方向、变更流程、检查项
- [ ] `history/architecture/README.md` 索引是否需要同步更新？
- [ ] 本次约定是否均为“已协商确认”内容（无未确认硬写）？
```

## 8. 可移植性提示（工具包发布）

- 本地 Python 解释器路径（如 `D:\...`）属于机器私有配置，不应作为项目契约写入仓库。
- 面向 `git clone` 的可扩展工具包，应优先依赖：
  - 环境创建说明（README/requirements）
  - 相对路径与参数化配置
  - 避免在仓库文档中固化个人机器绝对路径
