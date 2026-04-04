# water — 开发备忘

## 定位

水分析工作流：密度剖面、取向加权密度、吸附层检测与 θ 分布、三面板集成图。

## 约定

- `DEFAULT_START_INTERFACE = "normal_aligned"` — 从 +axis 界面开始的半路径分析
- 三面板图 (`Water.py`) 是 **two-pass**：第一遍算 density + orientation，第二遍算 theta PDF
- 输出 CSV/PNG 文件名常量在 `config.py` 中定义
- CSV 输出通过 `_write_csv_from_arrays()`��`utils/_io_helpers.py`），与 potential/charge 模块的 `_write_csv()` 格式一致
- Cell 参数获取：优先 `cell_abc` 参数，其次 `md_inp_path` 文件解析；两者均无则 `ValueError`
- matplotlib 延迟 import：`Water.py` 中 `matplotlib.ticker` 在 plot 函数体内导入，不在模块顶层

## 陷阱与历史 Bug

- `layer_tol_A` 必须穿透到所有内部调用（曾遗漏 — commit a96449b）
- Frame slicing 参数 (`frame_start/end/step`) 需要连接到所有分析函数（曾遗漏 — commit 937391c）

## 子目录

| 目录 | 用途 |
|---|---|
| `WaterAnalysis/` | 核心计算逻辑 → `WaterAnalysis/CLAUDE.md` |
