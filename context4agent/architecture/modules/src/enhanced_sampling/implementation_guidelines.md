# md_analysis.enhanced_sampling — Implementation Guidelines

## 职责边界

增强抽样方法的数据处理与可视化。当前实现覆盖慢增长（slow-growth）热力学积分。

## 目录组织

```
enhanced_sampling/
  __init__.py          # 空包（re-export 预留）
  slowgrowth/          # 慢增长热力学积分子包
    __init__.py        # re-export 数据类 + 分析/绘图/CSV 函数
    config.py          # 输出文件名常量
    SlowGrowth.py      # 数据类 + 积分逻辑
    SlowGrowthPlot.py  # 绘图 + CSV 导出 + 统一入口
```

## 依赖方向

- `enhanced_sampling.slowgrowth` → `utils.RestartParser.ColvarParser`（`ColvarMDInfo`）
- `enhanced_sampling.slowgrowth` → `utils.config`（`HA_TO_EV`）
- `enhanced_sampling.slowgrowth` → `utils._io_helpers`（`_write_csv`）
- **无反向依赖**：其他模块不依赖 `enhanced_sampling`

## 扩展性

未来可在 `enhanced_sampling/` 下新增同级子包（如 `metadynamics/`），结构与 `slowgrowth/` 平行。顶层 `__init__.py` 可在需要时添加 re-export。
