# Data Contract（`scripts/structure/utils` 当前实现约定）

> 仅记录当前代码中已经实现并在协作中使用的约定；若后续改实现，需同步更新本文档。

## 通用约定

- 处理对象为**单帧** `ase.Atoms`
- 原子索引统一为 **0-based**
- 与 z 方向分布相关的方法，命名统一包含 `z_distribution`
- `oxygen_indices` 输入可为 `(n,)` 或 `(n, 1)`，内部会展平为一维索引数组

## 默认配置（`scripts/structure/utils/config.py`）

- `DEFAULT_METAL_SYMBOLS`：默认金属集合为所有过渡金属（当前实现同时包含 `La/Ac` 与 `Lu/Lr`）
- `DEFAULT_Z_BIN_WIDTH_A = 0.25`（Angstrom）
- `DEFAULT_THETA_BIN_DEG = 5.0`（degree）
- `DEFAULT_WATER_OH_CUTOFF_A = 1.25`（Angstrom）
- `WATER_MOLAR_MASS_G_PER_MOL = 18.01528`（g/mol）

## `LayerParser.py` 接口约定

### 输入与默认值

- 主入口：`detect_interface_layers(atoms, ...)`
- `metal_symbols=None` 时，使用 `DEFAULT_METAL_SYMBOLS`
- `normal` 默认为 `"c"`（沿晶胞 c 方向）
- `layer_tol_A` 默认为 `0.6`
- `n_interface_layers` 默认为 `2`
- `nonmetal_symbols_hint=None` 时，默认把“非金属”定义为“不在 `metal_symbols` 中的全部原子”

### 输出结构

- 返回 `SurfaceDetectionResult`
- 关键字段：
  - `axis_unit`：用于排序/判向的单位向量
  - `metal_indices`：金属原子索引元组
  - `metal_layers_sorted`：按法向投影坐标排序后的层列表
- 每层 `Layer` 含：
  - `atom_indices`
  - `center_s`
  - `is_interface`
  - `normal_unit`（仅 `is_interface=True` 时有值）

### 判层与界面标记规则

- 先按法向投影坐标做 1D 聚类得到金属层
- 默认对低端与高端两侧各取最外层 `n_interface_layers` 作为界面层
- 当 `normal` 是 `"a" / "b" / "c"` 时，用分数坐标最小像判定界面法向正负；否则回退到几何侧判定

## `WaterParser.py` 接口约定

### 水分子标记

- `detect_water_molecule_indices(atoms, ...) -> (n_water, 3)`
- 每行格式为 `[O_idx, H1_idx, H2_idx]`
- O-H 连通性判定使用 MIC 距离与阈值 `oh_cutoff_A`
- 分配策略是“按 O-H 距离从小到大贪心匹配”，且每个 H 最多归属一个 O

### 氧索引提取

- `get_water_oxygen_indices_array(water_molecule_indices) -> (n_water, 1)`

### 质量密度 z 分布

- `compute_water_mass_density_z_distribution(atoms, oxygen_indices, ...) -> (nbins, 1)`
- 单位：`g/cm^3`
- 约定：一个 O 代表一个 H2O 分子

分箱与体积口径：

$$
nbins = \left\lceil \frac{L_z}{dz} \right\rceil,\quad
V_{\mathrm{bin}} = A_{xy}\,\Delta z,\quad
A_{xy} = \lVert \mathbf{a}\times\mathbf{b}\rVert
$$

密度口径：

$$
\rho_{\mathrm{bin}} = \frac{N_{O,\mathrm{bin}}\cdot M_{\mathrm{H_2O}}/N_A}{V_{\mathrm{bin}}}
$$

其中：

- z 坐标先 wrap 到 `[0, L_z)`
- `L_z = |\mathbf{c}|`

### 取向加权 z 分布

- `compute_water_orientation_weighted_density_z_distribution(atoms, oxygen_indices, ...) -> (nbins, 1)`
- 单位：`1/Angstrom^3`
- $\theta$ 定义：以 O 为起点、沿 H-O-H 角平分线方向的向量，与晶胞 `+c` 方向（即分数坐标第三轴正向，`001`）的夹角
- 实现口径：`cos(theta)` 使用 `dot(bisector, c_unit) / |bisector|`，其中 `c_unit = c / |c|`

统计口径：

$$
P_{\mathrm{orient},\mathrm{bin}} = \frac{\sum_{i\in \mathrm{bin}}\cos\theta_i}{V_{\mathrm{bin}}}
$$

- 若输入 `oxygen_indices` 中存在未被水分子标记覆盖的 O，会抛错提示

### 指定 c 轴区间的取向角 PDF

- `compute_water_orientation_theta_pdf_in_c_fraction_window(atoms, oxygen_indices, c_fraction_range, ...) -> (180 / ndeg,)`
- `c_fraction_range = [start_c_fraction, end_c_fraction]`，使用分数坐标第三轴区间筛选水氧
- `ndeg` 默认取 `DEFAULT_THETA_BIN_DEG = 5.0`，并要求 `180 / ndeg` 为整数
- 角度范围固定为 `[0, 180]` 度，输出为角度概率密度（单位 `degree^-1`）
- 区间筛选口径为周期区间（左闭右开）：
  - `start < end`：`[start, end)`
  - `start > end`：跨越边界，按 `[start, 1) U [0, end)`
  - `start == end (mod 1)`：按全区间处理
