# calibration — 开发备忘

## 定位

表面电荷密度 σ → 电极电势 φ 的标定映射。输入电势始终为 V vs SHE，输出可转换为 RHE/PZC。

## 文件结构

| 文件 | 用途 |
|---|---|
| `config.py` | 默认文件名、Nernst 常量、拟合方法常量 |
| `_data.py` | `CalibrationData` 数据类、CSV/JSON I/O |
| `_mapper.py` | `ChargePotentialMapper` ABC + Linear/Polynomial/Spline 实现 |
| `_plot.py` | 标定拟合可视化 |
| `CalibrationWorkflow.py` | 公开 API：`calibrate()`、`predict_potential()`、`convert_reference()` |

## 约定

- **输入电势单位**：始终 V vs SHE
- **输出电势参考**：通过 `convert_reference()` 转换 SHE/RHE/PZC
- **RHE 转换**：φ_RHE = φ_SHE + (RT/F)·ln(10)·pH
- **标定文件**：独立 JSON（默认 `~/.config/md_analysis/calibration.json`），同时存储原始数据点和拟合参数
- **scipy 依赖**：仅 `SplineMapper` 需要 scipy（lazy import）
- **图表风格**：与 `PhiZProfile.py` 一致（11×4.8, 160 dpi, grid, tight_layout）

## CLI 命令

- 231：Calibrate from CSV File
- 232：Calibrate from Manual Input
- 233：Predict Potential from Charge

## 与 charge 模块的集成

`surface_charge_analysis()`（charge/Bader/SurfaceCharge.py）在分析完成后自动尝试加载 `~/.config/md_analysis/calibration.json`：
- **有标定**：通过 `mapper_from_dict()` + `predict(σ)` 得到 V vs SHE 电势，再根据 `potential_reference` 参数通过 `convert_reference()` 转换参考标度（SHE/RHE/PZC），追加 φ 列到 CSV（列名含参考标度标签），PNG 添加右轴
- **无标定**：log info 提示用户可通过菜单 23 标定，不影响 σ 输出
- **参考配置**：CLI 221-223 从用户配置（菜单 909 `SetPotentialReferenceCmd` 设置）读取 `potential_reference`/`pH`/`temperature_K`/`phi_pzc` 传递给分析函数
