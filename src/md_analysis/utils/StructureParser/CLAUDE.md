# StructureParser — 开发备忘

## 定位

结构解析子包：金属层检测与界面标记、水分子拓扑推断、1D 周期性聚类。

## 约定

### 层排序（LayerParser）
- `metal_layers_sorted` 始终排列为 `[normal_aligned, interior..., normal_opposed]`
- 从 +axis 界面穿过 slab 到 -axis 界面
- `center_frac`：分数坐标 [0,1)，不是绝对 Å 坐标
- `interface_label` ∈ `{"normal_aligned", "normal_opposed", None}`

### 轴支持
- 仅接受 `"a"`/`"b"`/`"c"` 字符串，传入自定义向量或数值会 `ValueError`
- `AXIS_MAP = {"a": 0, "b": 1, "c": 2}`，`AREA_VECTOR_INDICES` 给出垂直平面的两个轴索引

### 周期性聚类（ClusterUtils）
- `cluster_1d_periodic`：贪心聚类 + **wrap-around merge**（首尾 cluster 可能跨 0/period 边界合并）
- 中心用 circular mean 计算（复数指数法），不是算术平均
- `find_largest_gap_periodic`：找最大间隙来定位水层区域

### 水分子检测（WaterParser）
- 返回 `(n_water, 3)` 数组：`[O_index, H1_index, H2_index]`，H 按原子索引排序
- O-H cutoff 默认 1.25 Å，每个 H 只能属于一个 O（贪心最短距离）
- **bisector 方向**：`OH1_unit + OH2_unit`（两个 O→H 单位向量之和），不是 H1-H2 中点方向
- θ = bisector 与参考方向（如 +c 轴单位向量）的夹角

## 陷阱与历史 Bug

- `center_frac` 曾叫 `center_s`（commit df29467 重命名），旧注释/文档中可能残留
- `"low_c"`/`"high_c"` 已重命名为 `"normal_aligned"`/`"normal_opposed"`
- 单层 slab 情况：仅标记为 `normal_aligned`，没有 `normal_opposed`
- `gap_midpoint_periodic` 正确处理跨 0/period 边界的 gap（不要用简单算术中点）
- z-bin 边界情况：当 `lz_A / dz_A < 1` 时退化为单 bin
