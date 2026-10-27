# scripts/utils — 开发备忘

## 定位

CP2K XYZ 与 VASP POSCAR 之间的双射索引映射（IndexMapper）。

## 约定

### POSCAR 注释行格式
```
md_analysis::v1 frame=<int> source=<urlencoded> n=<int> order=<csv> p2x=<base64>
```
- `p2x`：POSCAR→XYZ 方向的排列数组，base64 编码的 int16/int32（按原子数自动选择）
- `source`：URL 编码的源文件路径
- `order`：逗号分隔的元素排列顺序

### 双射不变量
- `poscar_to_xyz[xyz_to_poscar[i]] == i`（往返一致）
- `xyz_to_poscar[poscar_to_xyz[j]] == j`

### remap_array
- `direction="p2x"`：POSCAR 顺序 → XYZ 顺序
- `direction="x2p"`：XYZ 顺序 → POSCAR 顺序

## 陷阱

- `element_order` 决定 VASP 中元素分组顺序 — 与 ASE 的默认排序可能不同
- 注释行解码失败（格式不匹配）→ `IndexMapParseError`
- 该模块无 md_analysis 内部依赖（纯 ASE + numpy + urllib）
