# calibration — 接口暴露

> 对应代码：`src/md_analysis/electrochemical/calibration/__init__.py`
>
> 本文档定义 `md_analysis.electrochemical.calibration` 的符号级公开接口与暴露边界。

## 1. 接口角色定义

- `calibration` 是表面电荷密度 σ → 电极电势 φ 的标定映射子包。
- 输入电势始终为 V vs SHE，输出可通过 `convert_reference()` 转换为 RHE/PZC。
- 公开范围由 `__all__` 严格定义（15 个符号）。

## 2. 当前公开接口清单（按模块）

### 2.1 `_data.py` 数据 I/O

数据结构：

- `CalibrationData`
  - 标定数据容器（dataclass，非 frozen）
  - 字段：
    - `potentials: np.ndarray` — 电势数组（V vs SHE）
    - `charge_densities: np.ndarray` — 表面电荷密度数组（μC/cm²）
    - `reference: str` — 电势参考标度（默认 `"SHE"`）
    - `metadata: dict` — 附加元数据
  - 属性：`n_points: int` — 数据点数

函数：

- `load_calibration_csv(csv_path: str | Path) -> CalibrationData`
  - 输入：CSV 文件路径（列 1 = φ V vs SHE，列 2 = σ μC/cm²）
  - 自动检测首行是否为表头（尝试 float 转换）
  - 输出：`CalibrationData` 实例
- `load_calibration_json(json_path: str | Path | None = None) -> tuple[CalibrationData, dict[str, Any]]`
  - 输入：JSON 文件路径（默认 `~/.config/md_analysis/calibration.json`）
  - 输出：`(CalibrationData, fit_params)` 元组
- `save_calibration_json(data: CalibrationData, fit_params: dict[str, Any], json_path: str | Path | None = None) -> Path`
  - 输入：标定数据 + 拟合参数字典
  - 输出：保存的 JSON 文件路径

### 2.2 `_mapper.py` 映射器

数据结构：

- `FittingInfo`
  - 拟合质量容器（dataclass，非 frozen）
  - 字段：
    - `method: str` — 拟合方法名称
    - `r_squared: float` — 决定系数 R²
    - `residuals: np.ndarray` — 逐点残差（φ_predicted − φ_actual）
    - `rmse: float` — 均方根误差
    - `equation_str: str` — 人类可读方程字符串
    - `params: dict[str, Any]` — 方法特定拟合参数

抽象基类：

- `ChargePotentialMapper` (ABC)
  - σ → φ 映射的抽象接口
  - 抽象方法：
    - `fit(sigma: np.ndarray, phi: np.ndarray) -> FittingInfo` — 从标定数据拟合
    - `predict(sigma: np.ndarray | float) -> np.ndarray` — 预测 φ
    - `to_dict() -> dict[str, Any]` — 序列化为 JSON 兼容字典
    - `from_dict(d: dict[str, Any]) -> ChargePotentialMapper` (classmethod) — 从字典反序列化

具体实现：

- `LinearMapper`
  - 线性回归：`φ = slope · σ + intercept`
  - 序列化参数：`slope`, `intercept`
- `PolynomialMapper`
  - 多项式回归：`φ = Σ c_i · σ^i`
  - 构造参数：`degree: int = 2`
  - 序列化参数：`degree`, `coefficients`
- `SplineMapper`
  - 三次样条插值（需 scipy）
  - 延迟导入 `scipy.interpolate.CubicSpline`
  - 序列化参数：`sigma_data`, `phi_data`（原始数据点，重建时重新构造样条）
- `DifferentialCapacitanceMapper`
  - 分段线性 σ→φ 映射，基于逐区间微分电容
  - 算法：
    1. 按 φ 升序排列标定点
    2. 计算相邻点间微分电容 C_i = Δσ_i / Δφ_i (μF/cm²)
    3. 分段线性插值预测 φ；超出范围时使用最近端点的电容外推
  - 约束：σ 必须随 φ 单调递增（正微分电容）
  - 异常：`ValueError` — 重复 φ 值或 σ 非单调
  - 序列化参数：`phi_nodes`, `sigma_nodes`, `capacitances`

工厂函数：

- `create_mapper(method: str, **kwargs: Any) -> ChargePotentialMapper`
  - 按方法名称创建 mapper 实例
  - 支持方法：`"linear"`, `"polynomial"`, `"spline"`, `"differential_capacitance"`
  - `**kwargs` 传递给构造函数（如 `degree=3`）
  - 未知方法抛出 `ValueError`
- `mapper_from_dict(d: dict[str, Any]) -> ChargePotentialMapper`
  - 从序列化字典恢复已拟合的 mapper
  - 根据 `d["method"]` 分发到对应 `from_dict()` 类方法
  - 未知方法抛出 `ValueError`

### 2.3 `_plot.py` 可视化

函数：

- `plot_calibration(sigma: np.ndarray, phi: np.ndarray, mapper: ChargePotentialMapper, fitting_info: FittingInfo, *, reference: str = "SHE", output_dir: Path | None = None) -> Path`
  - 散点 + 拟合曲线 + R²/RMSE 标注
  - 图表尺寸：11×4.8 inch, 160 DPI
  - `DifferentialCapacitanceMapper` 特殊处理：自动标注各段 C 值和区间分界线
  - 输出：PNG 文件路径

### 2.4 `CalibrationWorkflow.py` 工作流

函数：

- `calibrate(csv_path: str | Path | None = None, *, data_points: list[tuple[float, float]] | None = None, method: str = DEFAULT_FITTING_METHOD, poly_degree: int = DEFAULT_POLY_DEGREE, output_dir: Path | None = None, calibration_json_path: Path | None = None) -> Path`
  - 完整标定流程：加载 → 拟合 → 保存 JSON → 输出 CSV → 绘图
  - `csv_path` 与 `data_points` 互斥（恰好提供一个）
  - 默认拟合方法：`"linear"`（`DEFAULT_FITTING_METHOD`）
  - 返回：保存的 JSON 文件路径

- `predict_potential(sigma: float | np.ndarray, *, calibration_json_path: Path | None = None, target_reference: str | None = None, temperature_K: float = DEFAULT_TEMPERATURE_K, pH: float = DEFAULT_PH, phi_pzc: float | None = None) -> np.ndarray`
  - 加载标定 JSON → 重建 mapper → 预测 φ → 可选参考转换
  - `target_reference` 为 None 时使用标定文件中存储的参考标度
  - 返回：预测电势数组（V vs target_reference）

- `convert_reference(phi: float | np.ndarray, *, from_ref: str, to_ref: str, temperature_K: float = DEFAULT_TEMPERATURE_K, pH: float = DEFAULT_PH, phi_pzc: float | None = None) -> np.ndarray`
  - 电势参考标度转换，以 SHE 为中枢
  - 支持标度：`"SHE"`, `"RHE"`, `"PZC"`
  - SHE → RHE：`φ_RHE = φ_SHE + (RT/F)·ln(10)·pH`
  - SHE → PZC：`φ_PZC = φ_SHE − φ_pzc`
  - PZC 转换需提供 `phi_pzc`，否则抛出 `ValueError`

### 2.5 `config.py` 常量

| 常量 | 类型 | 值 | 描述 |
|------|------|-----|------|
| `DEFAULT_CALIBRATION_DIR` | `Path` | `~/.config/md_analysis/` | 标定文件目录 |
| `DEFAULT_CALIBRATION_FILE` | `Path` | `~/.config/md_analysis/calibration.json` | 默认标定 JSON 路径 |
| `DEFAULT_CALIBRATION_CSV_NAME` | `str` | `"calibration_data.csv"` | 输出 CSV 文件名 |
| `DEFAULT_CALIBRATION_PNG_NAME` | `str` | `"calibration_fit.png"` | 输出 PNG 文件名 |
| `FITTING_LINEAR` | `str` | `"linear"` | 线性拟合方法标识 |
| `FITTING_POLYNOMIAL` | `str` | `"polynomial"` | 多项式拟合方法标识 |
| `FITTING_SPLINE` | `str` | `"spline"` | 样条拟合方法标识 |
| `FITTING_DIFFERENTIAL_CAPACITANCE` | `str` | `"differential_capacitance"` | 微分电容拟合方法标识 |
| `DEFAULT_FITTING_METHOD` | `str` | `"linear"` | 默认拟合方法 |
| `DEFAULT_POLY_DEGREE` | `int` | `2` | 默认多项式阶数 |
| `DEFAULT_REFERENCE` | `str` | `"SHE"` | 默认电势参考标度 |
| `DEFAULT_TEMPERATURE_K` | `float` | `298.15` | 默认温度 (K) |
| `DEFAULT_PH` | `float` | `0.0` | 默认 pH |
| `R_J_PER_MOL_K` | `float` | `8.314462618` | 气体常数 (J/mol/K, CODATA 2018) |
| `F_C_PER_MOL` | `float` | `96485.33212` | 法拉第常数 (C/mol, CODATA 2018) |
| `LN10` | `float` | `2.302585092994046` | ln(10) |

## 3. Internal consumers

- `electrochemical.charge.Bader.SurfaceCharge.surface_charge_analysis()` 内部调用 `load_calibration_json()` + `mapper_from_dict()` + `predict()` 实现自动电势外推
- CLI 命令 231/232/233 分别调用 `calibrate()`/`calibrate()`/`predict_potential()`

## 4. 推荐导入方式

- `from md_analysis.electrochemical.calibration import calibrate, predict_potential, convert_reference`
- `from md_analysis.electrochemical.calibration import CalibrationData, FittingInfo`
- `from md_analysis.electrochemical.calibration import LinearMapper, create_mapper, mapper_from_dict`

## 5. 异常与错误条件

- `ValueError` — `csv_path` 与 `data_points` 同时或均未提供
- `ValueError` — 未知拟合方法（`create_mapper`/`mapper_from_dict`）
- `ValueError` — 未知参考标度（`convert_reference`）
- `ValueError` — PZC 转换缺少 `phi_pzc`
- `ValueError` — `DifferentialCapacitanceMapper`：重复 φ 值或 σ 非单调
- `RuntimeError` — mapper 未拟合即调用 `predict()`
- `ImportError` — `SplineMapper` 缺少 scipy
- `FileNotFoundError` — 标定 JSON 不存在

## 6. 接口变更触发条件

以下情况属于 calibration 公开接口变更：

- `__all__` 新增/删除/重命名
- 公开函数签名变化
- Mapper 序列化格式变化（影响 JSON 兼容性）
- 默认常量值变化
- 新增 Mapper 实现类
