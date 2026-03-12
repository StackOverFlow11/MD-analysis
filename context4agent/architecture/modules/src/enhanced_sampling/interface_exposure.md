# md_analysis.enhanced_sampling — Interface Exposure

## 模块角色

增强抽样分析工作流的顶层包。当前仅包含 `slowgrowth` 子包（慢增长热力学积分）。

## Public API

当前 `enhanced_sampling/__init__.py` 未 re-export 任何符号（空包），所有公开接口通过子包 `slowgrowth` 直接导入。

## 推荐导入方式

```python
# 推荐：从 slowgrowth 子包导入
from md_analysis.enhanced_sampling.slowgrowth import (
    Slowgrowth,
    SlowgrowthFull,
    SlowgrowthSegment,
    slowgrowth_analysis,
    plot_slowgrowth_quick,
    plot_slowgrowth_publication,
    write_slowgrowth_csv,
)
```

## 顶层 re-export

`enhanced_sampling` **未**从 `md_analysis.__init__` re-export，需显式导入完整路径。
