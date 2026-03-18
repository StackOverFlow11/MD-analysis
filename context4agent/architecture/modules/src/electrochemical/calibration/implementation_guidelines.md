# calibration — 实现指南

## 定位

表面电荷密度 σ → 电极电势 φ 的标定映射模块。用于恒电荷 MD 中通过 Bader 电荷分析结果外推电极电势。

## 架构

```
calibration/
├── config.py              # 常量：Nernst 参数、默认文件名、拟合方法
├── _data.py               # CalibrationData 数据类、CSV/JSON I/O
├── _mapper.py             # ChargePotentialMapper ABC + Linear/Polynomial/Spline
├── _plot.py               # 标定拟合可视化
├── CalibrationWorkflow.py # 公开 API：calibrate(), predict_potential(), convert_reference()
└── __init__.py            # re-exports
```

## 数据接口

- **输入电势**：始终 V vs SHE
- **输出电势**：可通过 `convert_reference()` 转换为 RHE/PZC
- **电荷密度**：μC/cm²（与 charge 模块一致）
- **标定文件**：独立 JSON（默认 `~/.config/md_analysis/calibration.json`）

## 拟合方法

| 方法 | 类 | 依赖 |
|------|-----|------|
| linear | `LinearMapper` | numpy |
| polynomial | `PolynomialMapper` | numpy |
| spline | `SplineMapper` | scipy（lazy import） |

## CLI 命令

- 231: Calibrate from CSV File → `electrochemical/calibration/fit/`
- 232: Calibrate from Manual Input → `electrochemical/calibration/fit/`
- 233: Predict Potential from Charge → `electrochemical/calibration/predict/`
