# 当前重构计划：Breaking Refactor for engines + utils

> 维护要求：执行中每个 phase 完成后**勾掉对应 checkbox**；遇到偏离原计划的架构决策，直接更新本文档。
>
> 已确认取向：可以接受大范围破坏式重构，只要最终代码架构合理。`agent/` 层本轮不作为必须维护的公开接口；如果它阻碍主架构整理，可以先删除或停用，后续再重写。

---

## 核心判断

1. **不要为旧路径支付长期兼容成本**
   `md_analysis.utils.CubeParser`、`BaderParser`、`RestartParser`、`StructureParser` 这些旧 import 路径可以破坏。重构完成后，全仓统一到新路径，不保留一堆 wrapper。

2. **agent 层不要半维护**
   当前 agent contract 会把重构锁死在旧 FQN、旧 dataclass 名和旧模块布局上。本轮默认把 `agent/` 从验收范围移除；更推荐直接 `git rm src/md_analysis/agent` 以及相关 agent contract 测试和镜像文档，后续按新架构重写。

3. **核心分析工作流必须保住**
   本轮验收重点是 CLI、`md_analysis.main`、`water`、`electrochemical`、`enhanced_sampling`、`scripts`。这些路径要在每个 phase 后保持可 import、可测试。

4. **边界比兼容更重要**
   `utils/formats` 只做单文件解析；`utils/structure` 做几何/化学语义；`utils/io` 做路径发现和通用 IO；`engines` 做 CP2K/VASP 等计算引擎门面。

5. **移动、改名、行为改造分开做**
   先移动模块，再抽 dataclass，再改 CP2K engine 门面。不要把路径移动、类型重命名和业务逻辑重写混在一个 commit。

6. **外部副作用需要确认**
   `git push` 不作为自动步骤。最终本地验证完成后，等待用户确认再推远端。

---

## 非目标

- 不在本轮维护旧 agent contract。
- 不实现 VASP 解析，只放清晰的占位接口。
- 不保留旧 `utils.*Parser` import 兼容层，除非某个短期兼容点被明确证明必要。
- 不在 `engines` 里反向依赖 `electrochemical`、`water` 等上层分析模块。

---

## 目标终态

```
src/md_analysis/
├── engines/
│   ├── __init__.py
│   ├── protocols.py
│   ├── models.py
│   ├── cp2k.py
│   └── vasp.py
│
└── utils/
    ├── constants.py
    │
    ├── formats/
    │   ├── __init__.py
    │   ├── cube.py
    │   ├── bader.py
    │   ├── cp2k_colvar.py
    │   ├── cp2k_cell.py
    │   ├── cp2k_stdout.py
    │   ├── cp2k_xyz.py
    │   ├── vasp_report.py
    │   ├── vasp_outcar.py
    │   └── vasp_locpot.py
    │
    ├── structure/
    │   ├── __init__.py
    │   ├── layer.py
    │   ├── water.py
    │   ├── cluster.py
    │   └── cell.py
    │
    └── io/
        ├── __init__.py
        ├── frame_discovery.py
        ├── cell_resolver.py
        ├── csv_writer.py
        └── path_helpers.py
```

### 分层说明

| 层 | 职责 | 不应该做的事 |
|---|---|---|
| `utils/formats` | 解析单个文件或单类文件格式 | 不判断 workflow，不找目录，不持有 engine registry |
| `utils/structure` | 几何、层、水、团簇、cell 等结构语义 | 不读 CP2K/VASP 文件 |
| `utils/io` | 路径发现、通用读写、cell 文件选择调度 | 不解析具体文件内容 |
| `engines` | CP2K/VASP 门面、协议、engine-neutral model | 不依赖 `electrochemical` 等上层业务 |
| 上层业务模块 | 调用 engine 和 utils 完成分析工作流 | 不直接散落解析 CP2K stdout/cube/xyz 的细节 |

---

## 设计基线

| # | 决策 | 选择 |
|---|---|---|
| 1 | 重构策略 | Breaking refactor，统一新路径，不维护旧 import |
| 2 | agent 层 | 本轮删除或停用，不作为验收范围 |
| 3 | engines 范围 | CP2K 做全文件门面；VASP 只留占位 |
| 4 | dataclass 位置 | 放在 `engines/models.py`，作为 engine-neutral 类型 |
| 5 | dataclass rename | `ColvarRestart` → `ConstraintMetadata`；`LagrangeMultLog` → `LambdaSeries` |
| 6 | potential 相关类型 | `PotentialFrame`、`FermiRecord` 放到 `engines.models` 或中立层，避免 `engines` 依赖 `electrochemical` |
| 7 | cell parser | 建 `utils/formats/cp2k_cell.py`；不要把 cell 逻辑塞进 `cp2k_colvar.py` |
| 8 | VASP parser | 不默认注册到 auto discovery；显式调用时清楚地 `NotImplementedError` |
| 9 | 提交粒度 | 每个 phase 一个可 review commit，测试随 phase 跑 |
| 10 | push | 最终由用户确认后再执行 |

---

## Phase 0：基线确认与 agent 清场决策

**目的**：先把验收面收窄，避免后续每一步都被 agent contract 和旧 FQN 拖住。

- [ ] 记录当前工作树状态：`git status --short`
- [ ] 识别当前 modified / untracked 文件里哪些属于 agent contract 或 agent 镜像文档
- [ ] 明确本轮不维护 `agent/` 的公开接口
- [ ] 决定 agent 执行方式：
  - A. 删除整个 `src/md_analysis/agent`，同步删除相关 tests 和 context mirror
  - B. 只保留已经 contract-backed 且不依赖旧 `*_with_report` 的 task，删除 legacy agent task
  - C. 全保留 agent，但所有 handler 后续都改调 `workflows`
- [ ] 当前推荐 B：保留成熟 contract-backed 资产，同时不让旧 agent wrapper 阻碍底层重构
- [ ] 跑当前最小基线测试：`pytest test/unit -q`
- [ ] 记录当前失败项，如果失败来自 agent contract，归入本轮清场范围

**Acceptance**：本轮验收范围清楚；后续 phase 不再为了维护 agent 旧路径而扭曲设计。

---

## Phase 1：移除或收缩 agent 层

**目的**：用一次明确变更切断旧 agent contract 对架构重构的约束。

- [ ] 按 Phase 0 的 A/B/C 决策处理 `src/md_analysis/agent`
- [ ] 如果选 B：保留成熟 contract-backed task，删除依赖旧 `*_with_report` 的 legacy task
- [ ] 删除或停用已废弃 task 对应的 `test_agent_contracts.py`
- [ ] 删除或更新 `context4agent/architecture/modules/src/md_analysis/agent*` 镜像
- [ ] 检查 `src/md_analysis/CLAUDE.md`、`context4agent`、`docs` 中对 agent 的引用，改成“后续重写”
- [ ] 跑非 agent 单元测试

**Acceptance**：仓库不再宣称 legacy agent 层在本轮可用；失败测试不再来自已废弃的 agent contract。

---

## Phase 2：平移 `utils/formats`

**目的**：先完成单文件解析器的物理移位，不做行为改造。

- [ ] `CubeParser.py` → `utils/formats/cube.py`
- [ ] `BaderParser.py` → `utils/formats/bader.py`
- [ ] `RestartParser/ColvarParser.py` → `utils/formats/cp2k_colvar.py`
- [ ] `RestartParser/CellParser.py` 中的 cell 解析逻辑 → `utils/formats/cp2k_cell.py`
- [ ] 建 `utils/formats/__init__.py`
- [ ] 删除空的 `RestartParser/`
- [ ] 全仓更新 import 和 `lazy_import` 字符串
- [ ] 跑 `pytest test/unit -q`

**Acceptance**：所有单文件解析器进入 `utils/formats`；旧 `RestartParser` 目录不存在；业务行为不变。

---

## Phase 3：平移 `utils/structure`

**目的**：把几何/化学语义工具从旧 `StructureParser` 目录移到清晰的结构层。

- [ ] `StructureParser/LayerParser.py` → `utils/structure/layer.py`
- [ ] `StructureParser/WaterParser.py` → `utils/structure/water.py`
- [ ] `StructureParser/ClusterUtils.py` → `utils/structure/cluster.py`
- [ ] 建 `utils/structure/cell.py`，只放 cell 数据结构和无 IO helper
- [ ] 建 `utils/structure/__init__.py`
- [ ] 删除空的 `StructureParser/`
- [ ] 全仓更新 import 和 `lazy_import` 字符串
- [ ] 跑 `pytest test/unit -q`

**Acceptance**：结构语义集中在 `utils/structure`；不再从 `StructureParser` import。

---

## Phase 4：建立 `utils/io`

**目的**：把路径发现、通用写出、cell 文件选择调度从解析器里分出来。

- [ ] `_frame_discovery.py` → `utils/io/frame_discovery.py`
- [ ] `_io_helpers.py` → `utils/io/csv_writer.py`
- [ ] `cell_resolver.py` → `utils/io/cell_resolver.py`
- [ ] 抽 `utils/io/path_helpers.py`，放目录扫描、候选文件选择等通用逻辑
- [ ] 让 `cell_resolver.py` 调用 `utils/formats/cp2k_cell.py`，自己只负责调度
- [ ] 建 `utils/io/__init__.py`
- [ ] 全仓更新 import 和 `lazy_import` 字符串
- [ ] 跑 `pytest test/unit -q`

**Acceptance**：`utils/io` 只处理路径和通用 IO；具体 CP2K cell 解析不留在 resolver 里。

---

## Phase 5：建立 `engines` 骨架与中立 models

**目的**：先把 engine 边界立起来，再接业务模块。

- [ ] 建 `src/md_analysis/engines/{__init__.py,protocols.py,models.py,cp2k.py,vasp.py}`
- [ ] `engines/protocols.py`：迁入 `ConstraintMDParser`、`ParserInferenceError`
- [ ] `engines/models.py`：定义 engine-neutral dataclass
  - `ConstraintMetadata`
  - `LambdaSeries`
  - `ConstraintPoint`
  - `ConstraintSet`
  - `ConstraintRun`
  - `PotentialFrame`
  - `FermiRecord`
- [ ] `utils/formats/cp2k_colvar.py` 改为只解析文件并返回 `engines.models` 类型
- [ ] `engines/cp2k.py` 初步提供 constraint MD 读取门面
- [ ] `engines/vasp.py` 只放占位类和清晰错误，不进入自动推断
- [ ] 全仓更新 dataclass import 和类型名
- [ ] 跑 `pytest test/unit -q`

**Acceptance**：dataclass 单一真相源在 `engines.models`；`utils/formats` 不再定义跨 engine 的业务类型。

---

## Phase 6：改造 `enhanced_sampling` 使用 `engines`

**目的**：让 TI/SG 等约束 MD workflow 依赖 engine 协议，而不是私有 parser 模块。

- [ ] `enhanced_sampling/_parsers.py` 删除或变成短期内部迁移层
- [ ] `enhanced_sampling` 全部改为从 `md_analysis.engines` import
- [ ] `infer_parser` / `resolve_parser` 的逻辑移到 `engines/__init__.py` 或专门 registry 模块
- [ ] 保证 CP2K parser 行为不变
- [ ] 不把 VASP parser 放进自动推断路径
- [ ] 跑 enhanced sampling 相关 unit tests

**Acceptance**：约束 MD 工作流通过 `engines` 访问 CP2K；没有业务模块继续依赖 `enhanced_sampling._parsers`。

---

## Phase 7：扩展 `engines/cp2k.py` 为 CP2K 全文件门面

**目的**：把 21x potential 里散落的 CP2K stdout / xyz / cube 解析入口收敛到 CP2K engine。

- [ ] `utils/formats/cp2k_stdout.py`：迁入 Fermi、step、time 等 stdout 行解析
- [ ] `utils/formats/cp2k_xyz.py`：迁入 CP2K xyz 注释中的 step/frame 解析
- [ ] `engines/cp2k.py` 增加门面 API：
  - `read_constraint_metadata(directory) -> ConstraintMetadata`
  - `read_lambda_series(directory) -> LambdaSeries`
  - `read_cube_frames(directory, mode=...) -> list[PotentialFrame]`
  - `read_fermi_series(md_out_path) -> list[FermiRecord]`
- [ ] `electrochemical/potential/_frame_source.py` 改成调度 engine API，不直接持有 CP2K 解析细节
- [ ] 确认 `engines` 不 import `electrochemical`
- [ ] 跑 potential 相关 unit / integration tests

**Acceptance**：CP2K 文件格式细节集中在 `utils/formats`，CP2K 工作流入口集中在 `engines/cp2k.py`。

---

## Phase 8：VASP 占位接口

**目的**：给后续 VASP 协作者留清楚扩展点，但不制造“看起来能用”的假实现。

- [ ] `utils/formats/vasp_report.py`：占位函数，显式 `NotImplementedError`
- [ ] `utils/formats/vasp_outcar.py`：占位函数，显式 `NotImplementedError`
- [ ] `utils/formats/vasp_locpot.py`：占位函数，显式 `NotImplementedError`
- [ ] `engines/vasp.py`：实现协议形状，但所有读取方法显式报未实现
- [ ] 不在 `infer_parser` 中自动返回 VASP parser
- [ ] 加测试：显式构造 VASP parser 时抛出清楚错误

**Acceptance**：VASP 扩展点存在；不会被自动发现误用。

---

## Phase 9：文档同步

**目的**：让文档反映新的架构，而不是继续描述已删除目录。

- [ ] `src/md_analysis/CLAUDE.md`：更新目录树和模块边界
- [ ] `src/md_analysis/utils/CLAUDE.md`：重写为 formats / structure / io 三层结构
- [ ] 新建 `src/md_analysis/engines/CLAUDE.md`
- [ ] `src/md_analysis/enhanced_sampling/CLAUDE.md`：改为引用 `engines`
- [ ] 删除旧 `RestartParser/CLAUDE.md`、`StructureParser/CLAUDE.md`
- [ ] 更新 `context4agent/architecture/modules/src/` 镜像
- [ ] 扫描 `docs/workflows/*.md` 中的旧路径引用
- [ ] 标注 agent 层为后续重写，不再描述为当前可用接口

**Acceptance**：代码路径、文档路径、context mirror 三者一致。

---

## Phase 10：最终回归与提交整理

- [ ] `pytest test/unit -q`
- [ ] 可行时跑 `pytest test/integration -q`
- [ ] `rg "CubeParser|BaderParser|RestartParser|StructureParser|enhanced_sampling._parsers" src test docs context4agent`
- [ ] `rg "md_analysis.agent|test_agent_contracts" src test docs context4agent`
- [ ] 检查 `engines` 依赖方向：不能 import `electrochemical`、`water`、`enhanced_sampling`
- [ ] 整理 commit：每个 phase 独立、可 review
- [ ] 用户确认后再 `git push`

**Acceptance**：本地测试通过；旧路径残留清楚处理；远端推送由用户确认。

---

## 重点风险

| 风险 | 处理方式 |
|---|---|
| agent tests 阻塞重构 | Phase 0-1 明确删除或停用，不在后续 phase 中继续维护 |
| CLI `lazy_import` 字符串漏改 | 每个移动 phase 后用 `rg "lazy_import|CubeParser|RestartParser|StructureParser"` 扫描 |
| dataclass rename 影响面大 | 先把 models 搬到 `engines.models`，字段名尽量保持不变，只改类名 |
| `engines` 反向依赖上层业务 | 把 `PotentialFrame`、`FermiRecord` 放到中立 model；用 import 检查守住方向 |
| VASP 占位被误自动发现 | 不默认注册到 `infer_parser`；只允许显式调用，并明确抛错 |
| 21x potential 改造风险高 | 放到 Phase 7；前面先把 utils 和 models 稳住 |
| 文档与代码不同步 | Phase 9 集中扫 `src`、`test`、`docs`、`context4agent` |

---

## 建议验收命令

```bash
pytest test/unit -q
pytest test/integration -q
rg "CubeParser|BaderParser|RestartParser|StructureParser|enhanced_sampling._parsers" src test docs context4agent
rg "md_analysis.agent|test_agent_contracts" src test docs context4agent
```

如果 integration 当前本来就不稳定，以 Phase 0 记录的基线为准；不能把本轮重构新引入的失败混进已有失败里。

---

## 预计工作量

| Phase | 估时 |
|---|---|
| 0：基线确认与 agent 决策 | 30-45 min |
| 1：agent 清场 | 30-60 min |
| 2：formats 平移 | 1 h |
| 3：structure 平移 | 30-45 min |
| 4：io 拆分 | 1 h |
| 5：engines + models | 1.5-2 h |
| 6：enhanced_sampling 接入 engines | 1 h |
| 7：CP2K 全文件门面 | 2-3 h |
| 8：VASP 占位 | 30 min |
| 9：文档同步 | 1-2 h |
| 10：最终回归与提交整理 | 30-60 min |

**合计：约 10-14 小时**。建议拆成 2 个工作日。

---

## 建议暂停点

- Phase 1 后：agent 清场完成，后续可以专心改核心架构
- Phase 4 后：`utils` 三层结构基本完成，`engines` 还没大规模接入
- Phase 6 后：`enhanced_sampling` 已接入 `engines`，CP2K 全文件门面尚未动 21x
- Phase 7 后：CP2K 全文件门面落地，适合做一次完整 review

---

[← 回 README](README.md) | [→ short_term.md](short_term.md)
