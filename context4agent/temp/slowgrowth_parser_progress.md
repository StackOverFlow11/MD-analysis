# SlowgrowthParser 实现进展报告

## 完成状态：✅ 全部完成

## 新增文件

| 文件 | 说明 |
|---|---|
| `src/md_analysis/utils/RestartParser/__init__.py` | 子包接口，重导出 8 个公共符号 |
| `src/md_analysis/utils/RestartParser/SlowgrowthParser.py` | 核心解析模块（~300 行） |
| `test/unit/utils/test_slowgrowth_parser.py` | 14 个单元测试 |

## 修改文件

| 文件 | 变更 |
|---|---|
| `src/md_analysis/utils/__init__.py` | 添加 RestartParser 的 8 个符号导入 + `__all__` |
| `context4agent/.../utils/interface_exposure.md` | 添加 §2.8 RestartParser 接口文档 |
| `context4agent/.../utils/implementation_guidelines.md` | 添加 RestartParser 职责描述 |
| `context4agent/requirements/short_term.md` | 标记 SlowgrowthParser 为已完成 |

## 公共 API

### 数据类（frozen dataclass）

- `ColvarDef` — 集合变量定义（支持 ANGLE/DISTANCE/COMBINE_COLVAR 递归嵌套）
- `ConstraintInfo` — COLLECTIVE 约束参数
- `SlowGrowthRestart` — restart 文件元数据
- `LagrangeMultLog` — 拉格朗日乘子时序（含 `collective_shake`/`collective_rattle` 属性）

### 公共函数

- `parse_slowgrowth_restart(path)` → `SlowGrowthRestart`
- `parse_lagrange_mult_log(path)` → `LagrangeMultLog`
- `compute_target_series(restart, n_steps)` → `np.ndarray`

### 异常

- `SlowGrowthParseError`

## 测试结果

- 14/14 新测试通过
- 150/150 全套测试通过（无回归）
- 耗时 ~142s

### 测试覆盖

| 场景 | restart 解析 | LagrangeMultLog | 说明 |
|---|---|---|---|
| `angle/` | ✅ ANGLE CV | ✅ 单约束 | 3 原子角度约束 |
| `distance/` | ✅ DISTANCE CV (.bak-1) | ✅ 单约束 | 2 原子距离约束 |
| `distance_combinedCV/` | ✅ COMBINE_COLVAR (D1-D2) | ✅ 单约束 | 嵌套双距离组合 CV |
| `more_constrain/` | ✅ COMBINE_COLVAR + FIXED_ATOMS | ✅ 多约束 (97 constraints) | 范围展开 `7..10`、`\` 续行 |

- `compute_target_series`：验证 step_start 处 ξ==target_au，线性增长
- 边界：空文件报错、缺 &CONSTRAINT 报错、范围展开

## 设计要点

- **复用** `CellParser.parse_abc_from_restart()` 获取 cell 参数
- **递归** 解析 COMBINE_COLVAR 中嵌套的 &COLVAR 块
- **自动检测** LagrangeMultLog 格式（第 2 行是否以 Rattle 开头）
- **容错** 末尾不完整 Shake 块静默丢弃

## Git

- Commit: `5cda8d1` on `development`
- 已 push 到 origin
