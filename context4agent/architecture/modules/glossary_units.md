# Glossary / Units（约定汇总）

> 仅记录**已与用户协商确认**的约定；未确认内容必须标注"待讨论/待确认"，不得当作既定事实写入。
>
> 记录落位硬约束：本文件仅承载"全局术语/单位约定"；非全局目录级补充记录必须写入
> 对应目录的 `interface_exposure.md` / `implementation_guidelines.md`。

## 文档公式格式（已确认）

- 行内公式统一使用 `$...$`
- 块级公式统一使用 `$$...$$`
- 不使用 `\(...\)` 与 `\[...\]` 作为公式定界符

## 术语与统计口径索引（当前实现）

- 跨模块统计数据形状、单位、公式口径：`data_contract.md`
- `src` 层级接口暴露约定：`modules/src/**/interface_exposure.md`
- `src` 层级内部实现准则：`modules/src/**/implementation_guidelines.md`

## 体系前提（已确认）

- 三基矢正交的周期性体系
- 存在金属/水界面；由于周期性边界条件，体系总会有两个表面（两个界面）

## 界面层策略（已确认）

- `detect_interface_layers()` 固定每侧标记 1 层界面层，共 2 层，不可配置。
- 原有 `n_interface_layers` 参数已移除（该参数从未影响实现逻辑）。

## 单位约定（已确认）

| 物理量 | 单位 | 说明 |
|---|---|---|
| 水质量密度 $\rho(z)$ | `g/cm^3` | 一个 O 计作一个 H2O 分子 |
| 取向加权密度 $\rho(z)\cdot\langle\cos\theta\rangle(z)$ | `g/cm^3` | $\sum_i \cos\theta_i \cdot m_{\mathrm{H_2O}} / V_{\mathrm{bin}}$ |
| 角度 $\theta$ | degree | H-O-H 角平分线与 `+c_unit` 夹角 |
| 角度 PDF | `degree^-1` | 归一化概率密度 |
| 距离 | Angstrom | 沿法向的实空间距离 |
| 路径坐标 | 无量纲 `[0,1]` | 界面 → 两界面中点路径的归一化值 |

## 界面法向（已确认）

### normal_unit 的含义

- `normal_unit` 是一个**单位向量**（长度为 1），用于标记界面层的"朝向"（符号），仅在 `Layer.is_interface=True` 时出现。
- `normal_unit` 的方向定义为：**从金属界面层指向非金属环境（通常为水相）的一侧**。
- 在当前主流程默认使用 `normal="c"` 时，`normal_unit` 等价于 $\pm c_{\mathrm{unit}}$（$c_{\mathrm{unit}} = \mathbf{c} / |\mathbf{c}|$）。

### $\pm c_{\mathrm{unit}}$ 的判据（沿所选轴的分数坐标 MIC）

令：

- $f_z^{(m)}$：某个金属界面层的分数坐标 z
- $f_z^{(w)}$：与该界面层最近的水（参考氧）的分数坐标 z

定义 $\Delta f = f_z^{(w)} - f_z^{(m)}$，使用 MIC：

$$
\Delta f_{\mathrm{MIC}} = ((\Delta f + 0.5)\bmod 1) - 0.5 \in [-0.5, 0.5)
$$

判据：

- 若 $\Delta f_{\mathrm{MIC}} > 0$，取 `normal_unit = +c_unit`
- 若 $\Delta f_{\mathrm{MIC}} < 0$，取 `normal_unit = -c_unit`
- 若 $\Delta f_{\mathrm{MIC}} = 0$，退化情况（几乎不出现），需额外 tie-break（待确认）

> 注：体系存在两个表面，该判据对两侧表面分别计算时，通常得到一正一负。
