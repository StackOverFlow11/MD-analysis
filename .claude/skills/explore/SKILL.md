---
name: explore
description: Parallel exploration of the entire project for context building
user-invocable: true
---

# Parallel Project Exploration

When this skill is invoked, launch **all of the following agents in parallel** (in a single message with multiple Agent tool calls) to explore the project comprehensively. Each agent should be of type `Explore` with thoroughness "very thorough".

## Agent 1: Project Structure & Git State

Explore the overall project structure and recent git activity:
- Run `git log --oneline -20` to see recent commits
- Run `git status` to see current working state
- Run `git branch -a` to see all branches
- Read `pyproject.toml` for project configuration
- Read `CLAUDE.md` for project conventions
- Glob `src/md_analysis/**/__init__.py` to map the package tree
- Summarize: project structure, recent changes, current branch state

## Agent 2: Core Utils Module

Explore `src/md_analysis/utils/` in depth:
- Read all files in `utils/` including sub-packages (`StructureParser/`, `RestartParser/`)
- Identify all public functions, classes, and constants
- Note any recent changes or TODO comments
- Summarize: public API, key data structures, dependencies between sub-modules

## Agent 3: Water Analysis Module

Explore `src/md_analysis/water/` in depth:
- Read all files in `water/` and `water/WaterAnalysis/`
- Identify the full public API and internal helpers
- Trace the data flow from input (xyz + md.inp) to output (CSV + PNG)
- Summarize: analysis pipeline, key functions, configuration options

## Agent 4: Potential & Charge Modules

Explore `src/md_analysis/potential/` and `src/md_analysis/charge/` in depth:
- Read all files in both modules
- Identify public APIs, configuration constants, and output formats
- Understand the cSHE formula implementation and Bader charge methods
- Summarize: analysis pipelines, key functions, output structure

## Agent 5: CLI, Scripts & Tests

Explore `src/md_analysis/cli/`, `src/md_analysis/scripts/`, and `test/`:
- Read all CLI menu files to understand the user-facing interface
- Read scripts module (BaderGen, IndexMapper)
- Glob `test/**/*.py` and read test files to understand coverage
- Check `context4agent/` for any outstanding requirements or tech debt
- Summarize: CLI menu structure, script capabilities, test coverage gaps

## After All Agents Complete

Synthesize the results from all 5 agents into a **concise project status report** in Chinese (中文), covering:
1. 项目当前状态（分支、最近改动）
2. 各模块概览（utils / water / potential / charge / scripts / cli）
3. 测试覆盖情况
4. 待办事项或技术债务
5. 值得注意的代码模式或问题
