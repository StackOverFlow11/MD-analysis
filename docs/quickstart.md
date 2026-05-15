# Quickstart — 5 分钟跑通第一个分析

[← 回索引](README.md)

目标：在示例数据上跑一遍**水三联图分析**（菜单 105），看看 md-analysis 长什么样、输出落在哪里。

---

## 0. 安装

```bash
pip install numpy matplotlib ase pytest tqdm
pip install -e .   # editable 安装；开发者推荐
```

完成后会有一个 `md-analysis` 命令在你的 PATH 下。

> 如果不想用 editable 安装，可以用 `pip install .` 然后在测试 / 调用脚本时
> 加 `PYTHONPATH=src` 前缀，让 Python 优先加载工作区源码而不是 site-packages
> 里的旧版本。

---

## 1. 启动 CLI

```
$ md-analysis
================================================================================
            MD-Analysis: Metal-Water Interface Analysis Toolkit
                              Version 0.1.0
================================================================================

 ---------- MD-Analysis ----------

 1) Water Analysis
 2) Electrochemical Analysis
 3) Enhanced Sampling
 4) Scripts / Tools

 9) Settings

   0) Back / Exit

 Input:
```

---

## 2. 直达菜单 105（水三联图）

不需要逐级进菜单——**直接输 `105`** 就跳进去：

```
 Input: 105

 ---------- Full Water Three-Panel Analysis  (includes 101-104) ----------

 XYZ trajectory file: data_example/potential/dense/md-pos-1.xyz
 Cell parameters
   1) Auto from md.inp (../md.inp / ./md.inp)
   2) Manual
 Choice [1]: 1
   Detected ABC = (10.2239, 10.2239, 26.4220) A from md.inp

 Z-axis bin width (A) [0.1]:
 Modify advanced parameters? (y/n) [n]: n
```

按回车走默认值即可（方括号里是默认）。

---

## 3. 跑分析

按 Enter 后会看到进度条：

```
 Loading frames: 100%|██████████| 6/6 [00:01<00:00,  4.20it/s]
 Computing density: 100%|██████████| 6/6 [00:00<00:00, 13.50it/s]
 Computing orientation: 100%|██████████| 6/6 [00:00<00:00, 12.10it/s]

 Output:
   density:           ./output/water/water_density.csv
   orientation:       ./output/water/water_orientation.csv
   adsorbed-water:    ./output/water/ad_water_orientation.csv
   theta:             ./output/water/ad_water_theta.csv
   three-panel PNG:   ./output/water/water_three_panel.png

 Press Enter to continue...
```

---

## 4. 查看输出

```
$ tree output/
output/
└── water/
    ├── ad_water_orientation.csv
    ├── ad_water_theta.csv
    ├── water_density.csv
    ├── water_orientation.csv
    └── water_three_panel.png
```

PNG 是 4 联图（密度 / 取向加权密度 / 吸附水取向 / 吸附水 θ 分布），CSV 是各自的数据列。

---

## 5. 输出目录结构规则

`<outdir>` 是你输入的输出根目录（默认 `./output/`），下面**自动**按菜单路径分子目录：

```
<outdir>/
├── water/                              # 1xx
├── electrochemical/
│   ├── potential/                      # 21x
│   │   ├── center/  fermi/  electrode/
│   │   ├── phi_z/   thickness_sensitivity/
│   ├── charge/                         # 22x
│   │   ├── counterion/  layer/  full/
│   │   ├── counterion_aligned/
│   │   ├── tracked/  counterion_tracking/
│   └── calibration/                    # 23x
│       ├── fit/  predict/
└── enhanced_sampling/
    ├── slowgrowth/                     # 30x
    └── constrained_ti/                 # 31x
```

跑哪个菜单，**输出就只落在对应那一层**——你不会用 105 跑出来的水分析覆盖之前 213 跑出来的电极电势。

---

## 接下来

- 想跑真实分析？看 [workflows.md](workflows.md) 里的 5 条主线工作流
- 想知道每个菜单代码具体是什么？查 [menu_reference.md](menu_reference.md)
- CLI 的默认值（layer 容差 / bin 宽度等）想改？看 [settings.md](settings.md)
- 跑出来报错或结果异常？翻 [pitfalls.md](pitfalls.md)
