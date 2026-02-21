# Glossary / Units（约定汇总）

> 仅记录**已与用户协商确认**的约定；未确认内容必须标注“待讨论/待确认”，不得当作既定事实写入。
>
> 记录落位硬约束：本文件仅承载“全局术语/单位约定”；非全局目录级补充记录必须写入
> 对应目录的 `interface_exposure.md` / `implementation_guidelines.md`。

## 文档公式格式（已确认）

- 行内公式统一使用 `$...$`
- 块级公式统一使用 `$$...$$`
- 不使用 `\(...\)` 与 `\[...\]` 作为公式定界符

## 术语与统计口径索引（当前实现）

- 跨模块统计数据形状、单位、公式口径：`data_contract.md`
- `scripts` 层级接口暴露约定：`modules/scripts/**/interface_exposure.md`
- `scripts` 层级内部实现准则：`modules/scripts/**/implementation_guidelines.md`

## 体系前提（已确认）

- 三基矢正交的周期性体系
- 存在金属/水界面；由于周期性边界条件，体系总会有两个表面（两个界面）

## 界面法向（已确认）

### normal_unit 的含义

- `normal_unit` 是一个**单位向量**（长度为 1），用于标记界面层的“朝向”（符号），并且只会在 `Layer.is_interface=True` 时出现。
- `normal_unit` 的方向定义为：**从金属界面层指向非金属环境（通常为水相）的一侧**。
- 在当前主流程默认使用 `normal="c"` 时，`normal_unit` 等价于 `± c_unit`（`c_unit = c / |c|`）。
  - 若晶胞与笛卡尔坐标轴对齐，则可简写理解为 **001 / 00-1**；但更通用的表述应以 `± c_unit` 为准。

### `±c_unit` 的判据（沿所选轴的分数坐标 MIC）

以分数坐标（fractional coordinate）为例，令：

- $f_z^{(m)}$：某个金属“界面层”（例如最外层金属层中心）的分数坐标 z
- $f_z^{(w)}$：与该界面层“最近的水”（水的参考点待确定，例如水氧/水分子COM）的分数坐标 z

定义 $\Delta f = f_z^{(w)} - f_z^{(m)}$，使用 z 方向的最小像Minimum Image Convention（MIC）：

$$
\Delta f_{\mathrm{MIC}} = \Delta f - \mathrm{round}(\Delta f)
$$

等价实现（避免 0.5 边界舍入歧义）：

$$
\Delta f_{\mathrm{MIC}} = ((\Delta f + 0.5)\bmod 1) - 0.5 \in [-0.5, 0.5)
$$

判据：

- 若 $\Delta f_{\mathrm{MIC}} > 0$，则环境在金属层的 +c 最近像一侧，取 `normal_unit = +c_unit`
- 若 $\Delta f_{\mathrm{MIC}} < 0$，则环境在金属层的 -c 最近像一侧，取 `normal_unit = -c_unit`
- 若 $\Delta f_{\mathrm{MIC}} = 0$，属于退化情况（几乎不会出现）；需要额外 tie-break 规则（待确认）

> 注：由于体系存在两个表面，该判据对两侧表面分别计算时，通常会得到一正一负（对应两侧界面的相反法向）。
