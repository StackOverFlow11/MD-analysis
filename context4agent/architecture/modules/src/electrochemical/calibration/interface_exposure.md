# calibration — 接口暴露

## Public API (from `__init__.py`)

### Data
- `CalibrationData` — 标定数据容器（potentials, charge_densities, reference, metadata）
- `load_calibration_csv(csv_path)` → CalibrationData
- `load_calibration_json(json_path)` → (CalibrationData, fit_params)
- `save_calibration_json(data, fit_params, json_path)` → Path

### Mapper
- `ChargePotentialMapper` — ABC（fit, predict, to_dict, from_dict）
- `FittingInfo` — 拟合信息数据类（r_squared, rmse, residuals, equation_str, params）
- `LinearMapper`, `PolynomialMapper`, `SplineMapper`, `DifferentialCapacitanceMapper` — 具体实现
- `create_mapper(method, **kwargs)` — 工厂函数
- `mapper_from_dict(d)` — 从序列化字典恢复 mapper

### Workflow
- `calibrate(csv_path=, data_points=, method=, ...)` → Path（JSON）
- `predict_potential(sigma, calibration_json_path=, ...)` → ndarray
- `convert_reference(phi, from_ref=, to_ref=, ...)` → ndarray

### Visualization
- `plot_calibration(sigma, phi, mapper, fitting_info, ...)` → Path（PNG）

## Internal consumers

- `electrochemical.charge.Bader.SurfaceCharge.surface_charge_analysis()` 内部调用 `load_calibration_json()` + `mapper_from_dict()` + `predict()` 实现自动电势外推
