# 远期诉求（按需更新）

> 维护要求：仅当明确要求"更新远期诉求/路线图/愿景"时更新本文件。

## 路线图（草案）

- **从恒电荷金属-水界面拓展到更多体系**
  - 引入反应物，分析电化学过程相关指标（按未来研究方向补充）
- **更丰富的电势/电荷分析**
  - 更稳健的电势对齐与参考（真空/体相/滑动平均/误差条）
  - 不同电荷分析方案对比（Mulliken / Bader / Hirshfeld，如未来需要）
- **可复现与工程化**
  - CLI 工作流、配置文件（yaml/toml）、批处理整个项目目录
  - 单元测试覆盖率提升、文档站点、发布到 PyPI（如需要）

## 多引擎适配策略

- 当前主线：CP2K MD 产出结构 → VASP 单点/Bader 后处理；解析层（`BaderParser` 等）已按 VASP 格式实现
- 适配原则：保持解析层与分析层分离，解析层按引擎区分（如 `_read_acf` 对应 Bader/VASP），分析层（表面电荷密度、电荷转移统计）仅依赖 ASE Atoms + arrays，不绑定引擎
- 未来若需支持其他引擎（Gaussian、QE、ORCA 等）的电荷输出，只需新增对应解析函数，分析层无需修改

## 长期架构方向：数据接口与 parser registry

- **`engines.models.ConstraintMetadata` / `LambdaSeries` 当前是 re-export 过渡方案**：
  - 现状：`src/md_analysis/engines/models.py` 从 `utils.formats.cp2k_colvar` 导入 `ConstraintMetadata` / `LambdaSeries` 后重新暴露；运行时它们与 `utils.formats.cp2k_colvar.ConstraintMetadata` / `LambdaSeries` 是同一个类，不是继承关系。
  - 当前合理性：`utils.formats.cp2k_colvar` 是 CP2K 单文件解析层，直接构造这两个 dataclass；若立即把定义搬到 `engines.models`，会迫使 `utils` 运行时反向 import `engines`，破坏现有依赖方向（`engines -> utils`，禁止 `utils -> engines`）。
  - 架构风险：这两个名字语义上已经是 engine-neutral 类型，但物理定义仍在 CP2K-specific 文件解析模块里；未来实现 VASP parser 时，VASP 也要返回同一类 `ConstraintMetadata` / `LambdaSeries`，届时“neutral 类型住在 cp2k_colvar.py”会造成概念混淆。
  - 后续决策点：Phase 5 或真正启动 VASP 支持前，评估是否新建更底层的 neutral 类型模块（例如 `md_analysis.core.models` / `md_analysis.types`），让 `utils.formats.cp2k_colvar` 和 `engines.models` 都依赖该模块；在此之前，上层代码仍应统一从 `md_analysis.engines.models` 或 `md_analysis.engines` 导入 canonical 名。
- **后续 parser / protocol 抽象原则：按数据契约抽象，不按功能入口抽象**：
  - 不建议为每个 workflow / CLI 功能创建一个 parser（例如 `WaterDensityParser` / `TIFullAnalysisParser` 这类按功能命名的解析器），这会把业务功能和底层文件读取强绑定，并导致 CP2K/VASP 扩展时重复实现。
  - 推荐按上层业务共享的数据需求定义接口：例如 `ConstraintMDParser`（约束 MD 元数据 + λ(t)）、未来可评估的 `PotentialFrameReader`（电势帧 + Fermi 序列）、`ChargeTrajectoryReader`（电荷轨迹 / Bader / Mulliken 数据）等。
  - 抽象触发条件：同一类数据明确有多个后端（CP2K/VASP 等）、多个业务模块重复读取同一类原始数据、业务层被文件格式细节污染，或测试需要稳定 mock 数据源。否则优先保持 `engines.cp2k.read_*` 这类模块级 facade，避免过早抽象。
- **业务层依赖边界的长期取舍**：
  - 激进方向：业务层不应直接解析 engine/file-format 原始数据；CP2K/VASP/Bader/cube/stdout/restart 等文件读取与组合应先通过 `engines` 或未来的数据接口层结构化，再交给 `electrochemical` / `enhanced_sampling` / `water` 等业务模块。
  - 不宜绝对化为“业务层不能调用任何非 engines 内容”。业务层仍可直接依赖 engine-neutral 工具，例如 `utils.constants`（单位/物理常量）、`utils.structure`（界面层、水拓扑、周期聚类）和必要的通用 I/O helper。
  - 长期目标：减少业务层对 `utils.formats.*` 的直接依赖，尤其是 `cp2k_*` / `vasp_*` / `cube` / `bader` 等原始文件 parser；新增功能优先从 `engines` 的 dataclass / facade / protocol 获取结构化数据。
- **多数据接口 registry 的候选设计：泛型 `ParserRegistry[T]`（方案 B）**：
  - 当前 `engines.protocols` 只有 `_REGISTRY: dict[str, Callable[[], ConstraintMDParser]]`，本质上只服务 `ConstraintMDParser` 这一类数据接口。
  - 若未来新增 `PotentialFrameReader` / `ChargeTrajectoryReader` 等多个 protocol，优先评估一个可复用的泛型 registry 类，而不是为每类接口复制一套注册/查找/嗅探代码。
  - 目标形状示例：
    ```python
    class ParserRegistry(Generic[T]):
        def register(self, name: str, factory: Callable[[], T]) -> None: ...
        def get(self, name: str) -> T: ...
        def infer(self, path: Path) -> T: ...

    constraint_registry = ParserRegistry[ConstraintMDParser]()
    potential_registry = ParserRegistry[PotentialFrameReader]()
    ```
  - 采用前提：至少出现第二类稳定数据接口，且它同样需要按 engine 名称注册、按目录/文件自动嗅探、并被多个业务模块共享；否则继续保持当前单 registry + 模块级 facade，避免为当前代码引入过重抽象。

## 非目标（暂不考虑）

- 与 CP2K 之外的引擎做大规模统一 I/O 适配（仅在明确需求时逐步扩展）
