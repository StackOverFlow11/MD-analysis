# cli — 开发备忘

## 定位

VASPKIT 风格交互式编号菜单。无 argparse，所有输入通过 `input()` 提示。

## 约定

### 菜单编号方案
- `1xx`：Water (101-105)
- `21x`：Potential (211-216), `22x`：Charge (221-226), `23x`：Calibration (231-233)
- `30x`：Slow-Growth (301-302), `31x`：Constrained TI (311-313)
- `41x`：Bader (411-412), `42x`：TI (421-422)
- `9xx`：Settings (901-909)
- 编号前缀必须与父 MenuGroup 编号匹配

### 框架模式
- `MenuGroup`（非叶）/ `MenuCommand`（叶），在 `__init__.py` 的 `build_menu_tree()` 中组装
- `build_flat_index()` 使任意叶节点编号可从根直达（如直接输入 `421`）
- **flat index 只索引 `MenuCommand`**，不索引 `MenuGroup` — 输入子菜单编号（如 `42`）必须先进入父菜单

### 输出目录自动推导（output_name 机制）
- 每个 `MenuNode` 拥有 `output_name: str` 属性，贡献一段路径
- `MenuGroup` 通过构造参数 `output_name=` 设置（如 `MenuGroup("2", ..., output_name="electrochemical")`）
- `MenuCommand` 通过类属性或 `__init__` 中设置 `self.output_name`
- `MenuCommand.output_subdir` 是 `@property`，自动遍历父链拼接所有 `output_name` 段
- `MenuGroup.add()` 自动建立 `node.parent` 引用
- 移动命令到不同 `MenuGroup` 时输出路径自动更新
- 无 `output_name` 的节点（Scripts、Settings）不参与路径推导

### lazy_import
- 所有 `execute()` 方法通过 `lazy_import()` 延迟加载分析模块（36+ 处使用）
- 目的：CLI 启动只加载框架代码，不触发 numpy/matplotlib/ase

### 参数采集
- `K` 类：字符串键常量，防止拼写错误
- `ParamCollector` ABC：`collect(ctx)` 提示用户 + `apply_default(ctx)` 静默填充
- `params` 元组：总是提示；`advanced_params` 元组：用户选择"修改高级参数"时才提示
- `ConfigDefaultParam`：从 `~/.config/md_analysis/config.json` 读取用户覆盖值，fallback 到硬编码默认
- 所有 Potential 命令（211-216）的 `params` 元组首位为 `input_mode`（`ChoiceParam`："continuous"/"distributed"），后接模式相关参数（`sp_root_dir`、`sp_dir_pattern`、`sp_cube_filename`、`sp_out_filename`）。`execute()` 通过 `_is_distributed(ctx)` 分派调用

### 错误处理
- `MenuCommand.run()` 的 inline try-except 捕获 `MDAnalysisError`/`FileNotFoundError`/`ValueError`/`RuntimeError` → 打印简洁消息
- 未知异常 → `logger.error(..., exc_info=True)` 记录完整 traceback + 打印简洁消息到 stdout

### 测试钩子
- `_prompt.py` 的 `set_input_source(fn)` 可注入自定义输入函数

## 新增命令检查清单

1. 在 `_<module>.py` 中创建 `MenuCommand` 子类，定义 `params`/`advanced_params`/`output_name`/`execute()`
2. 在 `__init__.py` 中 import 并在 `build_menu_tree()` 中注册到对应 `MenuGroup`
3. 如需新参数键 → 在 `_params.py` 的 `K` 类中添加常量
4. 如需新参数类型 → 创建 `ParamCollector` 子类或使用现有泛型类

## 陷阱与历史 Bug

- 菜单码重编号（bace527）：旧代码中 401/402 已改为 411/412
- `CellAbcParam.collect()` 允许一次重试（.restart 失败 → 切 md.inp），第二次失败才报错
- `_discover_restart_file()` 排除 `_\d+.restart` 检查点文件（正则过滤）
- SG 命令会检测 LagrangeMultLog 中的 overflow（NaN 步），并在终端打印警告
- SG 命令 301/302 的 `output_name` 由父 `MenuGroup("30", output_name="slowgrowth")` 提供，`_SlowgrowthPlotCmd` 自身不定义 `output_name`（否则路径重复拼接为 `slowgrowth/slowgrowth`）
- `K.TI_DIR_PATTERN` 仅接受 `"ti_target"/"xi"/"auto"`，注意区分 `K.DIR_PATTERN`（Bader 用）
