# Bader — 开发备忘

## 定位

Bader 电荷分析子包：表面电荷密度计算、指定原子电荷追踪、counterion 自动检测追踪。

## 模块布局

| 文件 | 用途 |
|------|------|
| `_frame_utils.py` | 帧目录发现、排序、step/time 提取（私有） |
| `BaderData.py` | `BaderTrajectoryData` frozen dataclass + `load_bader_trajectory()` — 加载轨迹并 remap 回 XYZ 序 |
| `SurfaceCharge.py` | 表面电荷密度：单帧 `compute_frame_surface_charge`、多帧 `trajectory_surface_charge`、端到端 `surface_charge_analysis` |
| `AtomCharges.py` | 原子电荷：单帧/多帧索引查询（POSCAR 序）、指定原子追踪（XYZ 序）、counterion 逐帧检测追踪（XYZ 序） |

## 约定

### 索引序注意
- `frame_indexed_atom_charges` / `trajectory_indexed_atom_charges`：接受 **POSCAR 序**索引（历史 API，向后兼容）
- `tracked_atom_charge_analysis` / `counterion_charge_analysis`：接受/输出 **XYZ 序**索引（新 API）
- XYZ 序还原通过 `remap_array(data, imap, "poscar_to_xyz")` 实现，依赖 POSCAR 注释行中的 IndexMap 编码

### 相对导入
从本包内文件到 `md_analysis.utils` 需要 4 个点：`from ....utils.config import ...`

路径层级：`Bader` → `charge` → `electrochemical` → `md_analysis`

### BaderTrajectoryData
- `net_charges` 形状 `(n_frames, n_atoms)`，XYZ 序
- 不包含 `bader_charge`（与赝势相关的伪电荷无物理意义）
- `steps` 和 `times` 从目录名 `bader_t{time}_i{step}` 提取

## 陷阱

- `remap_array` 要求 `data.shape[0] == imap.n_atoms`，即 POSCAR 中的原子数必须与 IndexMap 编码一致
- POSCAR 注释行必须包含 `md_analysis::v1` 格式的 IndexMap 编码（由 BaderGen 写入）
- counterion 逐帧检测结果中，不同帧的原子集合可能不同（CSV 中未检测到的原子对应空值）
