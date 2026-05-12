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

1. **Integration fixture 路径未完全同步**
   - `data_example/potential/` 根目录下不再有 continuous fixture 文件。
   - 实际 continuous fixture 位于 `data_example/potential/dense/`。
   - 当前全量 integration 仍有若干测试指向旧路径，导致 `pytest test/integration`
     失败。

2. **验收命令容易 import 到已安装旧包**
   - 当前环境中直接执行 `pytest` 可能 import site-packages 里的旧
     `md_analysis`。
   - 本仓库内验收应使用 `PYTHONPATH=src pytest ...`，或先执行
     `pip install -e .`。

3. **CLI 表面仍未完全收敛到 workflows**
   - `cli/_charge.py`、`cli/_calibration.py`、`cli/_enhanced_sampling.py`、
     `cli/_constrained_ti.py`、`cli/_scripts.py` 仍有直接调用底层业务模块的路径。
   - 这不是 `main.py` facade 成功与否的阻塞项，但会让 CLI 与 notebook/API 的
     行为面继续分叉。

4. **旧计划路径引用需要清理**
   - 删除旧计划文件后，源码 docstring、CLAUDE/context 文档、测试说明不能再指向
     已删除的 `current_reconstructions.md` 或 `entrance_reconstruction.md`。

---

## Phase 0：冻结现状与验收基线

**目的**：确认“主体重构已完成”，避免后续 agent 重复执行旧迁移。

- [ ] 记录 `git status --short`
- [ ] 确认旧计划文件已删除：
  - `context4agent/requirements/current_reconstructions.md`
  - `context4agent/requirements/entrance_reconstruction.md`
- [ ] 确认新计划为唯一重构收尾入口：
  - `context4agent/requirements/refactor_repair_plan.md`
- [ ] 跑最小健康检查：
  - `PYTHONPATH=src pytest test/unit/workflows test/unit/engines -q`
  - `PYTHONPATH=src pytest test/integration/test_main.py -q`

**Acceptance**：小范围 workflows / engines / main facade 验收通过；没有 agent 再引用旧计划作为执行入口。

---

## Phase 1：修复 integration fixture 路径

**目的**：让全量 integration 测试反映当前 `data_example/potential/dense/` 布局。

- [ ] 更新 `test/integration/potential/test_center_potential.py`
  - `_DATA_DIR` 改为 `data_example/potential/dense`
  - docstring / skip reason 同步
- [ ] 更新 `test/integration/potential/test_phi_z_profile.py`
  - `_DATA_DIR` 改为 `data_example/potential/dense`
  - docstring / skip reason 同步
- [ ] 更新 water / utils integration 测试里的旧 fixture 路径：
  - `test/integration/water/test_*`
  - `test/integration/utils/test_water_layer_pipeline.py`
  - preview helper 脚本可顺手同步，但不作为 pytest 阻塞项
- [ ] 扫描确认：
  - `rg 'data_example.*potential.*md-(pos|POTENTIAL)|data_example.*potential.*md\\.inp|data_example.*potential.*md\\.out' test/integration`
- [ ] 跑：
  - `PYTHONPATH=src pytest test/integration -q`

**Acceptance**：全量 integration 不再因缺失 `data_example/potential/{md.inp,md.out,md-pos-1.xyz,cube}` 失败。

---

## Phase 2：统一验收命令与开发文档

**目的**：避免测试误用 site-packages 旧版本。

- [ ] 更新 `README.md` 的测试命令：
  - 推荐 `pip install -e .`
  - 或明确用 `PYTHONPATH=src pytest ...`
- [ ] 更新 `test/README.md` 中的测试命令说明
- [ ] 更新 `context4agent/requirements/README.md`，说明本文件是重构收尾计划
- [ ] 如有必要，更新 `docs/quickstart.md` 中与测试/本地运行相关的命令

**Acceptance**：新用户或 agent 按文档执行时，会 import 当前工作区源码，而不是环境中的旧安装包。

---

## Phase 3：清理旧计划引用

**目的**：删除旧计划后不留下坏链接。

- [ ] 扫描：
  - `rg 'current_reconstructions|entrance_reconstruction' src test docs context4agent README.md -g '!context4agent/requirements/refactor_repair_plan.md'`
- [ ] 将仍有意义的引用改为：
  - `context4agent/requirements/refactor_repair_plan.md`
- [ ] 对“重构进行中”的文字改成“主体已完成，剩余收尾见修补计划”
- [ ] VASP placeholder 的错误消息改为指向本文件中的 VASP 扩展说明

**Acceptance**：上述 `rg` 无结果，或只剩本文件中对旧文件的归档说明。

---

## Phase 4：CLI facade 收敛评估

**目的**：决定是否继续把 CLI 命令统一改走 `md_analysis.workflows`。

这一步不是主重构完成的阻塞项，可以单独安排。

- [ ] 列出 CLI 仍直接调用底层业务的命令：
  - charge
  - calibration
  - enhanced_sampling / constrained_ti
  - scripts
- [ ] 对每类命令判断：
  - 是否已有等价 `workflows.run_*`
  - CLI 是否需要打印 `WorkflowResult.artifacts`
  - 是否仍需要底层 report 的额外字段
- [ ] 优先迁移低风险命令：
  - charge 三个主入口
  - calibration fit / predict
  - scripts batch 入口
- [ ] 对复杂 TI / correction 命令，如果 workflow 入口还不能覆盖 CLI 细节，明确保留底层直调并写入 `cli/CLAUDE.md`

**Acceptance**：要么 CLI 主要命令统一走 workflows，要么文档明确说明哪些命令因交互细节暂时保留底层直调。

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
rg "enhanced_sampling._parsers" src test docs context4agent
```

如果不使用 `PYTHONPATH=src`，必须先运行：

```bash
pip install -e .
```

---

## 当前审查基线

最近一次审查记录：

- `PYTHONPATH=src pytest test/unit -q`：719 passed
- `PYTHONPATH=src pytest test/unit/workflows test/unit/engines -q`：105 passed
- `PYTHONPATH=src pytest test/integration/test_main.py -q`：8 passed
- `PYTHONPATH=src pytest test/integration -q`：33 passed, 11 failed

11 个 integration 失败集中在旧 `data_example/potential/` fixture 路径；这正是
Phase 1 的修补目标。
