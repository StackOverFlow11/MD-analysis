# 总体破坏性重构流程

> 用途：记录本轮大重构的总体顺序。具体技术方案由后续计划文件展开。

## Phase 0：备份当前上下文

重构前先备份当前状态，方便后续理解旧业务与回查测试数据。

建议备份内容：

- `src/`
- `context4agent/`
- `docs/`
- `data_example/`
- `test/`

推荐放在仓库外部，避免污染 `rg`、pytest 和 import 扫描：

```bash
cd /home/shaofl/Projects/git_projects
mkdir -p MD-analysis-backup/pre-reconstruction-YYYYMMDD
rsync -a MD-analysis/src MD-analysis/context4agent MD-analysis/docs \
  MD-analysis/data_example MD-analysis/test \
  MD-analysis-backup/pre-reconstruction-YYYYMMDD/
```

如果必须放仓库内 `backup/`，需加入本地 ignore，并在后续扫描中排除。

## Phase 1：重写 utils 边界

第一优先级是清理底层 `utils`：

- `utils.formats` 按 `cp2k` / `vasp` / `common` / `bader` 分包。
- `utils.io` 瘦身为真正通用 I/O helper。
- `utils.structure` 只保留结构/几何/拓扑职责，不读 CP2K/VASP 文件。
- 保留必要旧 import shim，避免一次性打断所有调用。

## Phase 2：梳理业务层数据需求

从顶层业务反推所需数据结构，而不是从文件格式往上堆接口。

- `water` 需要哪些结构/轨迹数据。
- `electrochemical` 需要哪些电势帧、电荷轨迹、结构字段。
- `enhanced_sampling` 需要哪些约束 MD metadata、lambda series、point collection。

## Phase 3：设计 engines 数据接口

先设计，给用户审阅，确认后再实现。

- 确定 `models`。
- 确定 `protocols`。
- 判断是否需要泛型 `ParserRegistry[T]`。
- 判断哪些数据接口暂不抽象，只保留 facade。

## Phase 4：实现 engines

在 `utils` 边界稳定后实现 engine adapter。

- 拆分/重写 `engines.cp2k`。
- 维持稳定 facade 导入路径。
- VASP 可保留 placeholder，不在本轮强行实现。

## Phase 5：迁移业务层

业务模块迁到新的 `engines` 数据接口。

- 减少业务层对 `utils.formats.*` 的直接依赖。
- 保持科学行为不变；需要测试保护关键数值结果。

## Phase 6：迁移入口层

最后处理入口层。

- `workflows` 先稳定。
- `agent` 跟随 workflows。
- `cli` 最后迁移，只保留参数收集和 workflow 调用。

## Phase 7：文档与测试收尾

主体重构完成后再统一维护：

- `context4agent`
- `docs`
- 全量测试结构
- 旧 shim 清理计划

测试可以阶段性保留最小回归集；完整文档同步可等主体结构稳定后再做。
