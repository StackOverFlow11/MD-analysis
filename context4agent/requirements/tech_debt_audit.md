# 技术债务审计报告（2026-03-09）

> 维护要求：每次解决一项 debt，将对应条目标记为 ✅ 并注明 commit hash 和日期。新发现的 debt 追加到对应分类末尾。

## 审计范围

对 `src/md_analysis/` 全部 ~45 源文件（~3,700 LOC）和 `test/` 全部 ~27 测试文件（~2,700 LOC）的完整代码审查，涵盖架构设计、代码质量、测试覆盖、性能、可观测性五个维度。

---

## 1. 测试覆盖缺口

### 1.1 零覆盖模块（Critical）

| 模块 | 文件 | 公共函数数 | 影响 |
|------|------|-----------|------|
| `cli/` | `__init__.py`, `_prompt.py`, `_water.py`, `_potential.py`, `_charge.py`, `_scripts.py`, `_settings.py` | ~30+ | 用户主交互界面完全无测试 |
| `main.py` | `main.py` | 4 (`run_water_analysis`, `run_potential_analysis`, `run_charge_analysis`, `run_all`) | 程序化 API 入口无验证 |
| `utils/CubeParser.py` | `CubeParser.py` | 5 (`read_cube_header_and_values`, `slab_average_potential_ev`, `plane_avg_phi_z_ev`, `z_coords_ang`, `extract_step_from_cube_filename`) | 底层解析器仅靠集成测试间接覆盖 |

### 1.2 覆盖不足的函数（High）

| 函数 | 位置 | 现状 |
|------|------|------|
| `thickness_sensitivity_analysis()` | `potential/CenterPotential.py` | 完全未测试 |
| `_savgol_smooth_window5()` | `water/Water.py` | 仅集成测试间接覆盖，无单元测试 |
| `detect_adsorbed_layer_range_from_density_profile()` | `water/WaterAnalysis/AdWaterOrientation.py` | 缺少边界情况单元测试 |
| `phi_z_planeavg_analysis()` | `potential/PhiZProfile.py` | 集成测试仅 65 行，仅验证文件存在 |

### 1.3 测试基础设施问题（Medium）

- **无 pytest markers**：缺少 `@pytest.mark.requires_data` 标记区分需要 `data_example/` 的测试
- **conftest.py 极简**：仅导出 1 个 helper，无共享 fixture（临时目录、mock Atoms 等）
- **测试数据构建器分散**：`_build_fake_trajectory()` 在 charge 测试中本地定义，未提取为共享 fixture
- **部分测试目录缺少 `__init__.py`**：`test/unit/`、`test/unit/utils/`、`test/unit/charge/`、`test/integration/` 等

---

## 2. 代码重复

### 2.1 Fortran 浮点解析 `_float()`（2 处）

- `utils/CubeParser.py:33-35` — `def _float(s): return float(s.replace("D","E").replace("d","e"))`
- `potential/CenterPotential.py:59-60` — 完全相同的函数

**建议**：提取到 `utils/config.py` 或新建 `utils/parsing.py`。

### 2.2 Circular Mean 计算（2 处）

- `utils/StructureParser/ClusterUtils.py:165-174` — `_circular_mean(values, period)`
- `utils/StructureParser/LayerParser.py:112-128` — `_circular_mean_fractional(f)`（period 硬编码为 1.0）

两者数学逻辑相同：`angles = 2π * f/period` → `atan2(mean(sin), mean(cos))`。

**建议**：统一为 `ClusterUtils._circular_mean(values, period)`，`LayerParser` 调用时传 `period=1.0`。

### 2.3 Cube 文件发现模式（4 处）

- `potential/CenterPotential.py` — 至少 3 处 `sorted(workdir.glob(cube_pattern))` + FileNotFoundError 检查
- `potential/PhiZProfile.py` — 1 处相同模式

**建议**：提取到 `utils/CubeParser.py` 新增 `discover_cube_files(workdir, pattern) -> list[Path]`。

### 2.4 界面检测调用模式（3 处）

- `charge/BaderAnalysis.py:118` 和 `:168`
- `water/WaterAnalysis/_common.py:65`

```python
det = detect_interface_layers(atoms, metal_symbols=metal_symbols, normal=normal)
iface_aligned = det.interface_normal_aligned()
iface_opposed = det.interface_normal_opposed()
```

三处代码块近乎相同。

**建议**：提取为 `utils/StructureParser/LayerParser.py` 中的便捷函数 `detect_both_interfaces(atoms, ...) -> tuple[Layer, Layer]`。

### 2.5 `_AXIS_MAP` 字典（2 处）

- `utils/StructureParser/LayerParser.py:41` — `_AXIS_MAP = {"a": 0, "b": 1, "c": 2}`
- `charge/BaderAnalysis.py:193` — 相同定义

**建议**：导出到 `utils/config.py` 作为 `AXIS_MAP`。

### 2.6 魔法字符串（散布多文件）

`"normal_aligned"`、`"normal_opposed"`、`"counterion"`、`"layer"` 等字符串在 10+ 处以字面量出现，无集中常量定义。

**建议**：在 `utils/config.py` 中定义：
```python
INTERFACE_NORMAL_ALIGNED = "normal_aligned"
INTERFACE_NORMAL_OPPOSED = "normal_opposed"
```

---

## 3. 错误处理与可观测性

### 3.1 异常体系不统一（Medium）

~~已修复~~（✅ 2026-03-09）。`exceptions.py` 定义 `MDAnalysisError(Exception)` 基类，9 个领域异常统一继承：

```
MDAnalysisError(Exception)          ← exceptions.py
├── ConfigError                     ← config.py
├── BaderParseError                 ← utils/BaderParser.py
├── WaterTopologyError              ← utils/StructureParser/WaterParser.py
├── SurfaceGeometryError            ← utils/StructureParser/LayerParser.py
├── CellParseError                  ← utils/RestartParser/CellParser.py
├── SlowGrowthParseError            ← utils/RestartParser/SlowgrowthParser.py
├── BaderGenError                   ← scripts/BaderGen.py
└── IndexMapError                   ← scripts/utils/IndexMapper.py
    └── IndexMapParseError          ← scripts/utils/IndexMapper.py
```

剩余问题：部分函数抛出通用 `ValueError`/`FileNotFoundError` 而非领域特定异常（更大改造，另行规划）。

### 3.2 CLI 无容错包装（High）

CLI 子菜单（`_water.py`、`_potential.py`、`_charge.py`、`_scripts.py`）的 `_cmd_*()` 函数无 try-except 包装。分析函数报错时直接 traceback 退出，无用户友好提示，无"重试/返回菜单"选项。

**建议**：在每个 `_cmd_*()` 函数外层添加统一的 try-except，捕获 `MDAnalysisError`、`FileNotFoundError`、`ValueError`，输出上下文提示后返回菜单。

### 3.3 零日志框架（High）

- 0 处 `import logging`，全部使用 `print()` 输出（CLI 中 11+ 处）
- 无结构化日志级别（DEBUG/INFO/WARNING/ERROR）
- 无日志文件输出
- 仅靠 `verbose: bool` 参数 + `tqdm` 进度条提供运行可见性

**建议**：引入 `logging` 模块，在 `__init__.py` 或 `main.py` 中配置根 logger `md_analysis`，替换所有 `print()` 为 `logger.info()`/`logger.warning()`。CLI 入口配置 log level。

---

## 4. 配置管理碎片化

### 4.1 现状

5 个独立的 `config.py` 文件：

| 文件 | 职责 | 条目数 |
|------|------|--------|
| `utils/config.py` | 物理常数 + 分析默认值 | ~20 |
| `water/config.py` | 输出文件名 + 界面默认值 | ~8 |
| `potential/config.py` | 输出文件名 + 厚度默认值 | ~10 |
| `charge/config.py` | 输出文件名 + 单位转换因子 | ~6 |
| 根 `config.py` | 持久化用户配置 | ~5 函数 |

### 4.2 问题

- 科学常数、分析默认值、输出文件名、用户设置混在不同层级
- 函数签名中硬编码默认值（如 `thickness_ang=7.5`）而非引用 `potential.config.DEFAULT_THICKNESS_ANG`
- 无配置值校验（负厚度、不存在的路径等不做检查）

### 4.3 建议（非紧急）

保持当前结构（各模块 config 独立），但：
- 确保函数签名默认值引用 config 常量而非字面量
- 在 `utils/config.py` 中集中定义跨模块共享的魔法字符串

---

## 5. 架构局限

### 5.1 正交晶胞假设深度嵌入（Architectural）

正交性检查和假设分布在 5 个核心模块中：

| 模块 | 硬编码位置 | 表现 |
|------|-----------|------|
| `RestartParser/CellParser.py` | `parse_abc_from_restart()` | off-diagonal > 1e-6 → raise CellParseError |
| `StructureParser/LayerParser.py` | `detect_interface_layers()` | 仅接受 `"a"/"b"/"c"` 字符串轴，拒绝自定义法向量 |
| `StructureParser/WaterParser.py` | 密度/取向分布函数 | `area_xy = |cross(a,b)|`，`lz = |c|` 假设矩形几何 |
| `water/WaterAnalysis/_common.py` | `_detect_interface_fractions()` | 硬编码 `normal="c"` |
| `utils/CubeParser.py` | `plane_avg_phi_z_ev()` | 隐含假设 cube grid 与坐标轴对齐 |

**影响**：hex、monoclinic、triclinic 体系完全不支持。

**改造估算**：需重构分数坐标投影、面积/体积计算、binning 逻辑，约 3-4 周工作量。

**当前决策**：根据 `short_term.md` 已确认的体系前提（"三基矢正交的周期性体系"），这是有意的设计约束，非 bug。仅在用户明确需要非正交体系时才启动改造。

### 5.2 无并行化基础设施（Performance）

- 0 处使用 `multiprocessing`/`ThreadPoolExecutor`/`concurrent.futures`
- 帧级分析天然可并行（每帧独立完成：界面检测 → 水分子拓扑 → binning）
- 水分析需要 2 次完整轨迹遍历（pass 1: 密度+取向, pass 2: θ 分布），可优化为单 pass
- Cube 文件完全加载到内存（`np.fromstring`），无 lazy/mmap 选项

**建议**：
1. 在 `_compute_density_orientation_ensemble()` 中添加可选 `n_workers` 参数
2. 评估单 pass 合并水分析的可行性（需在密度 pass 中同时累积 θ 数据）
3. 大型 cube 文件（>1GB）考虑 memory-mapped I/O

### 5.3 双 pass 水分析设计（Design Debt）

`water/Water.py:plot_water_three_panel_analysis()` 流程：
1. Pass 1：`_compute_density_orientation_ensemble()` → 密度 + 取向 CSV
2. 中间步骤：`detect_adsorbed_layer_range_from_density_profile()` → 确定吸附层范围
3. Pass 2：`compute_adsorbed_water_theta_distribution()` → θ 分布 CSV

Pass 2 依赖 Pass 1 的结果来确定吸附层范围，因此不能简单合并。但可以考虑：
- 在 Pass 1 中缓存每帧的原子坐标/取向数据
- 吸附层检测后直接从缓存中提取 θ 分布，避免重读轨迹

**权衡**：内存占用增加 vs I/O 时间减少。对典型体系（~300 原子，~2000 帧），轨迹文件 ~10MB，I/O 不是瓶颈；对大体系需评估。

---

## 6. 文档字符串缺口

### 公共函数缺少 docstring

| 函数 | 文件 | 状态 |
|------|------|------|
| `detect_water_molecule_indices()` | `utils/StructureParser/WaterParser.py` | ✅ 已有 docstring |
| `get_water_oxygen_indices_array()` | `utils/StructureParser/WaterParser.py` | ✅ 已有 docstring |
| `_read_xyz_metal_z_for_steps()` | `potential/CenterPotential.py` | ✅ 已有 docstring |
| CLI 各 `_cmd_*()` 函数（15 个） | `cli/_water.py`, `cli/_potential.py`, `cli/_scripts.py`, `cli/_settings.py` | ✅ 2026-03-09 已补全 |
| `_plot_series_with_cumavg()` | `potential/CenterPotential.py` | ✅ 2026-03-09 已补全 |
| `_parse_csv_symbols()` | `potential/CenterPotential.py` | ✅ 2026-03-09 已补全 |
| `_smooth_1d()` | `water/WaterAnalysis/AdWaterOrientation.py` | ✅ 2026-03-09 已补全 |
| `_plot_surface_charge()` | `charge/BaderAnalysis.py` | ✅ 2026-03-09 已补全 |

---

## 7. 待办清单

> 与第 1–6 节各条目一一对应，汇总完成状态。

### P0 — 立即修复（影响正确性或阻碍开发）

| # | 对应章节 | 待办项 | 状态 |
|---|---------|--------|------|
| 1 | 1.1 零覆盖模块 | 为 `CubeParser.py` 5 个公共函数编写单元测试 | ✅ 2026-03-09 `test/unit/utils/test_cube_parser.py` (22 tests) |
| 2 | 1.1 零覆盖模块 | 为 `main.py` 4 个入口函数编写集成测试 | ✅ 2026-03-09 `test/integration/test_main.py` (8 tests) |
| 3 | 2.1 Fortran 浮点解析 | 提取 `_float()` Fortran 浮点解析到共享位置 | ✅ 2026-03-09 `CenterPotential.py` 导入 `CubeParser._float` |

### P1 — 短期改善（提升可维护性和用户体验）

| # | 对应章节 | 待办项 | 状态 |
|---|---------|--------|------|
| 4 | 3.2 CLI 无容错包装 | CLI 子菜单 `_cmd_*()` 添加统一 try-except 容错 | ✅ 2026-03-09 `_handle_cmd_error` 装饰器统一 14 个 handler |
| 5 | 3.3 零日志框架 | 引入 `logging` 框架，替换 `print()` 调用 | ✅ 2026-03-09 NullHandler + per-module logger + CLI StreamHandler |
| 6 | 3.1 异常体系不统一 | 创建 `MDAnalysisError` 异常基类，统一继承关系 | ✅ 2026-03-09 `exceptions.py` + 9 个异常类统一继承 |
| 7 | 2.3 Cube 文件发现模式 | 提取 cube 文件发现逻辑到 `CubeParser.discover_cube_files()` | ✅ 2026-03-09 新增 `discover_cube_files()` + 4 个单测 |
| 8 | 1.2 覆盖不足的函数 | 为 `thickness_sensitivity_analysis()` 补充测试 | [ ] |

### P2 — 中期优化（提升性能和代码质量）

| # | 对应章节 | 待办项 | 状态 |
|---|---------|--------|------|
| 9 | 5.2 无并行化基础设施 | 添加帧级并行化选项（`ProcessPoolExecutor`） | [ ] |
| 10 | 5.3 双 pass 水分析设计 | 评估水分析单 pass 合并可行性 | [ ] |
| 11 | 1.1 零覆盖模块 | CLI 模块测试（mock stdin/stdout） | [ ] |
| 12 | 2.6 魔法字符串 | 集中管理界面标签等魔法字符串为常量 | ✅ 2026-03-09 `AXIS_MAP`/`AREA_VECTOR_INDICES`/`INTERFACE_*`/`CHARGE_METHOD_*` |
| 13 | 2.2 Circular Mean | 统一 `_circular_mean` 实现 | ✅ 2026-03-09 `LayerParser._circular_mean_fractional` 委托 `ClusterUtils._circular_mean` |

### P3 — 长期架构演进（需求驱动）

| # | 对应章节 | 待办项 | 状态 |
|---|---------|--------|------|
| 14 | 5.1 正交晶胞假设 | 非正交晶胞支持（仅在用户明确需求时启动） | [ ] |
| 15 | 4. 配置管理碎片化 | 配置体系重构（统一 schema + 校验） | [ ] |
| 16 | — | 分析模块插件化注册机制 | [ ] |
