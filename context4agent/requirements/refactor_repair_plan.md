# 重构收尾修补计划：utils/engines + main/workflows

> 本文件取代并归档原 `current_reconstructions.md` 与
> `entrance_reconstruction.md`。这两个重构主计划已经完成主体实现，不再作为
> 当前执行依据；历史细节从 git 历史查看。
>
> 当前目标不是再做一次大迁移，而是修补完成后留下的测试、文档、CLI 表面和验收命令
> 不一致问题。

---

## 当前结论

两个重构的主体架构已经落地：

- `src/md_analysis/main.py` 已变成 import-only facade。
- `src/md_analysis/workflows/` 已成为程序化入口的 canonical 层，统一返回
  `WorkflowResult`。
- `src/md_analysis/utils/formats/`、`utils/structure/`、`utils/io/` 已替代旧的
  `CubeParser` / `BaderParser` / `RestartParser` / `StructureParser` 路径。
- `src/md_analysis/engines/` 已建立 CP2K facade、parser registry、engine-neutral
  dataclass re-export 和 VASP placeholder。
- `agent/` 已收缩为保留的 contract-backed 任务集合，不再承载 legacy
  `main.py` wrapper。

当前不应按旧重构计划重复移动文件。后续工作应只处理下面列出的收尾缺口。

---

## 已接受的实现偏离

这些点与旧计划文字不同，但按当前代码判断可以接受；后续不要误当成未完成任务。

| 主题 | 旧计划写法 | 当前实现 | 处理 |
|---|---|---|---|
| `utils/io` 文件名 | `frame_discovery.py` / `csv_writer.py` / `path_helpers.py` | `_frame_discovery.py` / `_io_helpers.py` / `cell_resolver.py` | 暂时接受当前命名；只有真正需要公开稳定 API 时再改名 |
| `ConstraintMetadata` / `LambdaSeries` 定义位置 | `engines.models` 单一真相源 | 定义在 `utils.formats.cp2k_colvar`，由 `engines.models` re-export | 接受，避免 `utils -> engines -> utils` 循环 |
| CP2K potential facade 名称 | `read_cube_frames(directory, mode=...)` | `read_continuous_potential_frames` / `read_distributed_potential_frames` | 接受，名称更贴近 continuous / distributed 两种输入模式 |
| `ConstraintPoint` / `ConstraintSet` / `ConstraintRun` | 预留在 `engines.models` | 未实现 | 接受，等真实消费者出现再加 |
| CLI 全部走 workflows | 原 Phase 6 全量目标 | 仅 water three-panel 与 full potential 已走 workflows，其余仍有底层直调 | 作为本修补计划的可选收尾项追踪 |

---

## 当前已知缺口

> 状态摘要（2026-05-12）：Phase 0–3 ✅ 完成（commits `ea30387`、`886a013`）；
> Phase 4 / Phase 5 仍按下文计划保留为后置项。

1. ~~**Integration fixture 路径未完全同步**~~ — ✅ Phase 1 已修
   （commit `ea30387`，11 个测试文件从 `data_example/potential/` 改为
   `data_example/potential/dense/`，integration 全套 44 passed）。

2. ~~**验收命令容易 import 到已安装旧包**~~ — ✅ Phase 2 已修
   （commit `886a013`，`README.md` / `test/README.md` / `docs/quickstart.md`
   推荐 `pip install -e .` 或显式 `PYTHONPATH=src`）。

3. **CLI 表面仍未完全收敛到 workflows**（Phase 4 后置项）
   - `cli/_charge.py`、`cli/_calibration.py`、`cli/_enhanced_sampling.py`、
     `cli/_constrained_ti.py`、`cli/_scripts.py` 仍有直接调用底层业务模块的路径。
   - 这不是 `main.py` facade 成功与否的阻塞项，但会让 CLI 与 notebook/API 的
     行为面继续分叉。后续视真实需求决定是否做。

4. ~~**旧计划路径引用需要清理**~~ — ✅ Phase 3 已修
   （commit `886a013`，源码 docstring / CLAUDE / context4agent / 测试说明全部
   指向本文件；`rg "current_reconstructions|entrance_reconstruction"` 在 src /
   test / docs / context4agent 范围零命中本文件除外）。

---

## Phase 0：冻结现状与验收基线 ✅ 完成

**目的**：确认“主体重构已完成”，避免后续 agent 重复执行旧迁移。

- [x] 记录 `git status --short`
- [x] 确认旧计划文件已删除：
  - `context4agent/requirements/current_reconstructions.md`
  - `context4agent/requirements/entrance_reconstruction.md`
- [x] 确认新计划为唯一重构收尾入口：
  - `context4agent/requirements/refactor_repair_plan.md`
- [x] 跑最小健康检查：
  - `PYTHONPATH=src pytest test/unit/workflows test/unit/engines -q` → 105 passed
  - `PYTHONPATH=src pytest test/integration/test_main.py -q` → 8 passed

**Acceptance**：小范围 workflows / engines / main facade 验收通过；没有 agent 再引用旧计划作为执行入口。

---

## Phase 1：修复 integration fixture 路径 ✅ 完成

**目的**：让全量 integration 测试反映当前 `data_example/potential/dense/` 布局。

实现：commit `ea30387`（11 个测试文件，34+/24-）。

- [x] 更新 `test/integration/potential/test_center_potential.py`
  - `_DATA_DIR` 改为 `data_example/potential/dense`
  - docstring / skip reason 同步
- [x] 更新 `test/integration/potential/test_phi_z_profile.py`
  - `_DATA_DIR` 改为 `data_example/potential/dense`
  - docstring / skip reason 同步
- [x] 更新 water / utils integration 测试里的旧 fixture 路径：
  - `test/integration/water/test_*`（5 个测试文件）
  - `test/integration/utils/test_water_layer_pipeline.py`
  - 3 个 preview helper 脚本顺手同步
- [x] 扫描确认：
  - `rg 'data_example.*potential.*md-(pos|POTENTIAL)|data_example.*potential.*md\\.inp|data_example.*potential.*md\\.out' test/integration` → 仅 dense/ 路径
- [x] 跑：
  - `PYTHONPATH=src pytest test/integration -q` → **44 passed**（之前 11 failed / 33 passed）

**Acceptance**：✅ 全量 integration 不再因缺失 `data_example/potential/{md.inp,md.out,md-pos-1.xyz,cube}` 失败。

---

## Phase 2：统一验收命令与开发文档 ✅ 完成

**目的**：避免测试误用 site-packages 旧版本。

实现：commit `886a013`。

- [x] 更新 `README.md` 的测试命令：
  - 推荐 `pip install -e .`（"Install" section + "All tests" block）
  - 显式 `PYTHONPATH=src pytest ...` 作为非 editable 安装时的前缀
- [x] 更新 `test/README.md` 中的测试命令说明（两选项 + 显式警告非 editable + 无前缀的坑）
- [x] 更新 `context4agent/requirements/README.md`，说明本文件是重构收尾计划
- [x] 更新 `docs/quickstart.md` 中与测试/本地运行相关的命令

**Acceptance**：✅ 新用户或 agent 按文档执行时，会 import 当前工作区源码，而不是环境中的旧安装包。

---

## Phase 3：清理旧计划引用 ✅ 完成

**目的**：删除旧计划后不留下坏链接。

实现：commit `886a013`（绝大部分由用户手动改完，CC 整合 + 提交）。

- [x] 扫描：
  - `rg 'current_reconstructions|entrance_reconstruction' src test docs context4agent README.md -g '!context4agent/requirements/refactor_repair_plan.md'` → 0 命中
- [x] 将仍有意义的引用改为：
  - `context4agent/requirements/refactor_repair_plan.md`
- [x] 对“重构进行中”的文字改成“主体已完成，剩余收尾见修补计划”
- [x] VASP placeholder 的错误消息改为指向本文件中的 VASP 扩展说明（`engines/vasp.py` + 3 个 `utils/formats/vasp_*.py`）

**Acceptance**：✅ 上述 `rg` 无结果（本文件归档说明除外）。

---

## Phase 4：CLI facade 收敛评估 ✅ 完成

**目的**：决定是否继续把 CLI 命令统一改走 `md_analysis.workflows`。

执行结论：本轮迁移 11 个 CLI 命令到 workflows facade（合并 Phase 6 Step A
的 2 条共 13 条已走 facade）；其余 20 条因 workflow 签名 gap 或交互流程
复杂保留底层直调，分析命令总数 33。完整对照表与原因详见
`src/md_analysis/cli/CLAUDE.md` 的 "Workflows facade 迁移状态" 段。

### 本轮已迁移到 workflows（11 条新增 + 2 条 Phase 6 Step A 已迁，共 13 条）

| CLI 编号 | 命令 | workflow facade |
|---|---|---|
| 225 | `TrackedChargeCmd` | `workflows.charge.run_tracked_charge` |
| 226 | `CounterionChargeCmd` | `workflows.charge.run_counterion_charge` |
| 301 | `SGQuickPlotCmd` | `workflows.enhanced_sampling.run_slowgrowth_quick_plot` |
| 302 | `SGPublicationPlotCmd` | `workflows.enhanced_sampling.run_slowgrowth_publication_plot` |
| 411 | `BaderSingleCmd` | `workflows.scripts.run_bader_single` |
| 412 | `BaderBatchCmd` | `workflows.scripts.run_bader_batch` |
| 421 | `TISingleCmd` | `workflows.scripts.run_ti_single` |
| 431 | `PotentialSingleCmd` | `workflows.scripts.run_potential_single` |
| 432 | `PotentialBatchCmd` | `workflows.scripts.run_potential_batch` |
| 441 | `SpGenSingleCmd` | `workflows.scripts.run_sp_single` |
| 442 | `SpGenBatchCmd` | `workflows.scripts.run_sp_batch` |
| 105 | `WaterThreePanelCmd` | `workflows.water.run_water_three_panel`（Phase 6 Step A） |
| 216 | `FullPotentialCmd` | `workflows.potential.run_potential_full`（Phase 6 Step A） |

### 保留底层直调（workflow facade 当前 gap）

| CLI 编号 | 命令 | gap |
|---|---|---|
| 101–104 | water 单步 | workflows.water 只 export `run_water_three_panel`（composite），无 density / orientation / adsorbed / theta 单步 facade |
| 211–215 | potential 单步 | workflows.potential 只 export `run_potential_full`（composite），无 center / fermi / electrode / phi_z / thickness_sensitivity 单步 facade |
| 221 / 222 / 223 | `SurfaceChargeCmd` | `run_surface_charge` 不接受 `potential_reference` / `potential_pH` / `potential_temperature_K` / `potential_phi_pzc`（底层 `surface_charge_analysis` 已支持） |
| 224 | `SingleSideChargeCmd` | `run_surface_charge` 不接受 `target_side`（底层支持） |
| 231 / 232 / 233 | calibration fit / predict | workflow 要求 `calibration_json_path` 必填；CLI 当前允许传 None 走全局默认 `~/.config/md_analysis/calibration.json` |
| 311 / 312 / 313 | constrained TI | Python-slice 切片 UI、逐点 equilibration override、约束点交互列表等交互流程；workflow `run_ti_full_analysis` / `run_ti_constant_potential_correction` 是端到端单次调用，目前不能完整覆盖 |
| 422 | `TIBatchCmd` | `generate_ti_batch_with_report`（即 workflow `run_ti_batch`）不接受 `colvar_id`（MVP 限定 primary CV） |

### 后续工作（不在本 Phase 范围）

要继续往 workflows 收敛，前提是先扩展 workflow 签名覆盖以上 gap，
而不是在 CLI 端 hack 绕过。建议优先级（如未来要做）：

1. `run_surface_charge` 加 `target_side` + 4 个 `potential_*` kwargs → 覆盖 charge 221–224
2. `run_calibration_fit` / `run_calibration_predict` 允许 `calibration_json_path=None` 时使用 config 默认
3. 给 `workflows.{water,potential}` 加单步 facade（5 + 4 个 leaf workflow）
4. 给 `run_ti_batch` 加 `colvar_id`
5. 给 TI workflow 暴露 discover_only / select-points 流程后再迁 311–313

**Acceptance**：CLI 主要低风险命令统一走 workflows（13/33 已迁）；剩余 20 个
命令的保留原因在 `cli/CLAUDE.md` 的迁移状态表中显式记录，未来扩展 workflow
时按此 gap 列表对照修补。

---

## Phase 5：VASP placeholder 后续扩展记录

**目的**：保留 VASP 扩展点，但不让 placeholder 看起来可用。

- [ ] 保持 `engines.vasp.VASPParser` 不自动注册
- [ ] 保持 `utils/formats/vasp_report.py` / `vasp_outcar.py` /
  `vasp_locpot.py` 显式抛 `NotImplementedError`
- [ ] 后续真正实现 VASP 时，需要同时补：
  - parser 文件格式测试
  - engine facade 测试
  - 自动发现策略测试
  - docs/workflows 中的 VASP 数据契约说明

**Acceptance**：显式构造 VASP placeholder 仍报清楚错误；`infer_parser` 不会误返回 VASP stub。

---

## 最终验收命令

```bash
PYTHONPATH=src pytest test/unit -q
PYTHONPATH=src pytest test/integration -q
PYTHONPATH=src pytest test/unit/workflows test/unit/engines -q
PYTHONPATH=src pytest test/integration/test_main.py -q
rg "current_reconstructions|entrance_reconstruction" src test docs context4agent README.md -g '!context4agent/requirements/refactor_repair_plan.md'
rg "md_analysis.utils.(CubeParser|BaderParser|RestartParser|StructureParser)" src test docs context4agent
rg "enhanced_sampling._parsers" src test docs context4agent -g '!context4agent/requirements/refactor_repair_plan.md'
```

> 最后一条扫描 `enhanced_sampling._parsers` 时只允许“历史说明”类命中
> （例如 `"该 shim 已经在 Phase 6 删除"`、`"早期由 _parsers.py 适配"`），
> 不允许当前架构描述里把它当作活跃模块。当前 canonical 入口为
> `md_analysis.engines`（`engines.protocols.ConstraintMDParser` + `engines.cp2k.CP2KParser`）。

如果不使用 `PYTHONPATH=src`，必须先运行：

```bash
pip install -e .
```

---

## 当前审查基线

最近一次审查记录（Phase 0–3 收尾完成后，2026-05-12）：

- `PYTHONPATH=src pytest test/unit -q`：**719 passed**
- `PYTHONPATH=src pytest test/unit/workflows test/unit/engines -q`：**105 passed**
- `PYTHONPATH=src pytest test/integration/test_main.py -q`：**8 passed**
- `PYTHONPATH=src pytest test/integration -q`：**44 passed**（Phase 1 修复了之前的 11 failed）

Phase 4（CLI facade 收敛）已在 commit `58ce053` 完成；唯一仍开放的工作是
Phase 5（VASP placeholder 后续扩展记录），作为 backlog 按需启动。
Phase 0–4 不应再被任何 agent 重新执行。
