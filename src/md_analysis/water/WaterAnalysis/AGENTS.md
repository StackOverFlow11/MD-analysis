# WaterAnalysis — 开发备忘

## 定位

水分析核心计算子包。`_common.py` 提供共享的帧处理和 ensemble 平均逻辑。

## 约定

### Ensemble 平均流程
1. 逐帧读取轨迹（`_iter_trajectory`），设置 cell + PBC
2. 每帧检测两侧界面分数坐标（`_detect_interface_fractions`）
3. 每帧计算半路径上的密度/取向剖面（`_single_frame_density_and_orientation`）
4. 所有帧的 path_length 取平均 → 公共 bin 数
5. 每帧剖面归一化到 [0,1] → 插值到公共网格 → 逐 bin 等权平均

### 半路径
- Profile 范围：选定界面 → 两界面中点（gap midpoint）
- **不是**全路径穿越 slab — 只取水层的一半
- `StartInterface = Literal["normal_aligned", "normal_opposed"]`

### 吸附层检测
- `detect_adsorbed_layer_range_from_density_profile(distance, rho)`
- 规则：主峰（argmax）→ 下界（最后一个 near-zero bin）→ 上界（smoothed local minimum）
- smoothing 用 moving average（window=5），不是 Savitzky-Golay

### 所有函数都是 `_` 前缀（私有）
- 外部应通过 `WaterDensity.py`/`WaterOrientation.py`/`AdWaterOrientation.py` 的公开函数调用

## 陷阱与历史 Bug

- 两界面分数坐标相同 → `SurfaceGeometryError`（无法构建水层间隙）
- 某帧无水分子 → `ValueError`（跳过或报错取决于调用层）
- 插值用 `np.interp`（1D 线性），要求源/目标 x 都是单调递增
- path_length ≤ 0 → bin geometry 计算出错
