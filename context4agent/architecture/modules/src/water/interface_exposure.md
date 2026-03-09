# `md_analysis.water` 接口暴露约定（当前实现）

> 对应代码（对外导出）：`src/md_analysis/water/__init__.py`
>
> 主要实现文件：
> - `src/md_analysis/water/Water.py`
> - `src/md_analysis/water/WaterAnalysis/*`
>
> 本文档定义 `md_analysis.water` 子包的公开接口边界。

## 1. 接口角色定义

- `md_analysis.water` 是高层水分析工作流子包。
- 本层负责多帧统计组织与结果导出（如 CSV/PNG）。
- 本层可组合 `md_analysis.utils` 的底层能力，但不重写其物理定义。
- 子目录文档必须镜像代码目录：`water/WaterAnalysis/` 由独立接口/实现文档维护。

## 2. 当前公开接口清单

### 2.1 配置常量（Stable）

- `DEFAULT_OUTPUT_DIR`（`None`，语义：`None` 表示使用 `Path.cwd()` 作为输出目录）
- `DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME`
- `DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME`
- `DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME`
- `DEFAULT_START_INTERFACE`（`"normal_aligned"`）
- `DEFAULT_WATER_MASS_DENSITY_CSV_NAME`
- `DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME`
- `DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME`

注：各分析函数的 `dz_A` 参数默认值取自 `md_analysis.utils.config.DEFAULT_Z_BIN_WIDTH_A`（当前为 `0.1` Angstrom），它是函数参数而非模块级常量。

### 2.2 Re-exported Utility Types（Stable）

以下符号从 `md_analysis.utils` 转发，包含在 `water/__init__.py` 的 `__all__` 中：

数据结构与异常：

- `Layer`（dataclass，来自 `utils.StructureParser.LayerParser`）
- `SurfaceDetectionResult`（来自 `utils.StructureParser.LayerParser`）
- `SurfaceGeometryError`（来自 `utils.StructureParser.LayerParser`）
- `WaterTopologyError`（来自 `utils.StructureParser.WaterParser`）

函数：

- `detect_interface_layers`（来自 `utils.StructureParser.LayerParser`）
- `format_detection_summary`（来自 `utils.StructureParser.LayerParser`）
- `detect_water_molecule_indices`（来自 `utils.StructureParser.WaterParser`）
- `get_water_oxygen_indices_array`（来自 `utils.StructureParser.WaterParser`）

常量：

- `DEFAULT_METAL_SYMBOLS`（来自 `utils.config`）
- `DEFAULT_Z_BIN_WIDTH_A`（来自 `utils.config`）
- `DEFAULT_THETA_BIN_DEG`（来自 `utils.config`）
- `DEFAULT_WATER_OH_CUTOFF_A`（来自 `utils.config`）
- `WATER_MOLAR_MASS_G_PER_MOL`（来自 `utils.config`）

### 2.3 分析函数（Stable）

- `water_mass_density_z_distribution_analysis(xyz_path, md_inp_path, *, output_dir, output_csv_name, start_interface, dz_A, metal_symbols, frame_start, frame_end, frame_step) -> Path`
  - 文件位置：`src/md_analysis/water/WaterAnalysis/WaterDensity.py`
  - 统计口径：A 口径（逐帧等权系综平均）
  - 返回：CSV 路径

- `water_orientation_weighted_density_z_distribution_analysis(xyz_path, md_inp_path, *, output_dir, output_csv_name, start_interface, dz_A, metal_symbols, frame_start, frame_end, frame_step) -> Path`
  - 文件位置：`src/md_analysis/water/WaterAnalysis/WaterOrientation.py`
  - 统计口径：A 口径（逐帧等权系综平均）
  - 返回：CSV 路径

- `detect_adsorbed_layer_range_from_density_profile(distance_A, rho_g_cm3, *, near_zero_ratio, smoothing_window_bins) -> tuple[float, float, float]`
  - 文件位置：`src/md_analysis/water/WaterAnalysis/AdWaterOrientation.py`
  - 主峰定义：直接取水密度分布最高 bin 位置
  - 返回：`(start_distance_A, end_distance_A, peak_distance_A)`

- `ad_water_orientation_analysis(xyz_path, md_inp_path, *, output_dir, output_profile_csv_name, output_range_txt_name, start_interface, dz_A, near_zero_ratio, smoothing_window_bins, frame_start, frame_end, frame_step) -> tuple[Path, Path]`
  - 文件位置：`src/md_analysis/water/WaterAnalysis/AdWaterOrientation.py`
  - 功能：自动识别吸附层区间并输出取向分析结果
  - 返回：`(profile_csv_path, range_txt_path)`

- `compute_adsorbed_water_theta_distribution(xyz_path, md_inp_path, *, adsorbed_range_A, output_dir, output_csv_name, start_interface, dz_A, ndeg, near_zero_ratio, smoothing_window_bins, frame_start, frame_end, frame_step, verbose) -> tuple[ndarray, ndarray, Path]`
  - 文件位置：`src/md_analysis/water/WaterAnalysis/AdWaterOrientation.py`
  - 功能：统计吸附层内 `0-180` 度取向分布（PDF）
  - 返回：`(theta_centers_deg, theta_pdf_degree_inv, csv_path)`

- `plot_water_three_panel_analysis(xyz_path, md_inp_path, *, output_dir, output_png_name, start_interface, dz_A, ndeg, frame_start, frame_end, frame_step, verbose) -> Path`
  - 文件位置：`src/md_analysis/water/Water.py`
  - 功能：集成三联图（密度、取向、吸附层 $\theta$ 分布）并输出 PNG，同时落盘相关 CSV/TXT
  - 轨迹读取：固定两次（第一次密度+取向；第二次吸附层 $\theta$ 分布）
  - 返回：PNG 路径

## 3. 推荐导入方式

- `from md_analysis.water import water_mass_density_z_distribution_analysis`
- `from md_analysis.water import DEFAULT_START_INTERFACE`
- `from md_analysis.water import plot_water_three_panel_analysis`

## 4. 非公开边界

- 未在对应 `__init__.py` 的 `__all__` 中声明的符号，不属于公开接口。
- 内部 helper（如 `_parse_*`、`_single_frame_*`）仅供模块内复用。
- `src/md_analysis/water/WaterAnalysis/_common.py` 为私有实现模块，仅限 `water` 包内部使用（`Water.py` 和 `WaterAnalysis/` 子模块可导入），不属于对外契约。

## 5. 接口变更流程（必须执行）

1. 更新 `src/md_analysis/water/**/__init__.py` 导出与 `__all__`
2. 更新本文档公开接口清单
3. 若涉及统计口径变化，同步更新：
   - `context4agent/architecture/modules/data_contract.md`
4. 运行导入烟雾测试与相关回归测试
