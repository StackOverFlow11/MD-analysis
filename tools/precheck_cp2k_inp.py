#!/usr/bin/env python3
"""CP2K input pre-check: one-shot diagnostic for refactor step 4.

Runs the cit-based parser/generator (``md_analysis.utils.formats.cp2k.inp``)
over the four input faces defined in ``temp/cp2k_inp_refactor_plan.md`` §5
step 4 and writes a quantified report to ``temp/precheck_report.md``.

Input faces
-----------
1. ``data_example/``   every ``*.inp`` / ``*restart*`` file (incl. ``.bak``);
   the 121/121 baseline doubles as the script's self-check.
2. User real samples   inp/restart paths from the user config
   (``~/.config/md_analysis/config.json``; keys ``sp_inp_template_path`` /
   ``dp_sp_inp_template_path`` plus any string value naming an inp/restart).
   Unconfigured or missing-on-disk entries are reported but do NOT count
   as failures.
3. Packaged templates  ``importlib.resources`` enumeration of
   ``md_analysis.scripts:template/*.inp``.  ``sp_potential.inp`` still
   carries the ``{{KIND_BLOCKS}}`` placeholder at this stage and is a
   KNOWN EXPECTED FAILURE (re-run this face after step 5 legalises the
   template).  The placeholder is NOT a preprocessing target; the W1-W3
   pipeline stays untouched.
4. 03_CuAg corpus      every ``*.inp`` parsed; every COMBINE_COLVAR file
   gets a full ``generate_inp`` roundtrip (built-in re-parse self-check);
   5-file CRLF sampling verifies ``\\r\\n`` normalisation.

Faces 2/4 whose paths are absent (e.g. on another machine) are skipped
and noted in the report.  Exit codes: 0 = zero unexpected failures,
1 = unexpected failures found, 2 = setup error.
"""
from __future__ import annotations

import argparse
import importlib.resources
import json
import sys
import time
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
REPORT_PATH = REPO_ROOT / "temp" / "precheck_report.md"
DATA_EXAMPLE_DIR = REPO_ROOT / "data_example"
CORPUS_ROOT = Path("/home/shaofl/Projects/research_projects/03_CuAg_1_SAA_CO2RR")

DATA_EXAMPLE_BASELINE = 121  # plan §0: 27 inp + 94 restart = 121 (incl. .bak)
EXPECTED_COMBINE_COLVAR = 66  # plan §5.4 / §0 corpus measurement
PLACEHOLDER_TEMPLATE = "sp_potential.inp"
PLACEHOLDER_TOKEN = "{{KIND_BLOCKS}}"
CRLF_SAMPLE_N = 5


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _err_class(exc: Exception) -> str:
    """Group failures by exception type + first line of the message."""
    msg = str(exc).splitlines()
    first = msg[0] if msg else ""
    return f"{type(exc).__name__}: {first[:120]}"


def _group(failures: list[tuple[Path | str, Exception]]) -> dict[str, list[Path | str]]:
    groups: dict[str, list[Path | str]] = {}
    for src, exc in failures:
        groups.setdefault(_err_class(exc), []).append(src)
    return groups


def _parse_safe(src):
    """parse_inp; returns (tree, None) or (None, exc)."""
    from md_analysis.utils.formats.cp2k.inp import parse_inp

    try:
        return parse_inp(src), None
    except Exception as exc:  # diagnostic tool: surface everything
        return None, exc


def _roundtrip_safe(tree: dict):
    """generate_inp (self-checking); returns None or exc."""
    from md_analysis.utils.formats.cp2k.inp import generate_inp

    try:
        generate_inp(tree)
        return None
    except Exception as exc:
        return exc


def _git_info() -> str:
    import subprocess

    try:
        branch = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        sha = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return f"{branch}@{sha}"
    except Exception:
        return "<not a git checkout>"


def _dep_versions() -> str:
    bits = []
    try:
        import cp2k_input_tools  # noqa: F401

        bits.append(f"cp2k-input-tools {getattr(cp2k_input_tools, '__version__', '?')}")
    except Exception as exc:
        bits.append(f"cp2k-input-tools unavailable ({exc})")
    try:
        import pint

        bits.append(f"pint {pint.__version__}")
    except Exception as exc:
        bits.append(f"pint unavailable ({exc})")
    return ", ".join(bits)


# ---------------------------------------------------------------------------
# face 1: data_example (baseline self-check)
# ---------------------------------------------------------------------------


def collect_data_example() -> list[Path]:
    files = [
        p for p in DATA_EXAMPLE_DIR.rglob("*")
        if p.is_file()
        and (p.name.lower().endswith(".inp") or "restart" in p.name.lower())
    ]
    return sorted(files)


def run_face1(unexpected: list[tuple[str, Path | str, Exception]]) -> dict:
    files = collect_data_example()
    ok, failures = 0, []
    for path in files:
        _, exc = _parse_safe(path)
        if exc is None:
            ok += 1
        else:
            failures.append((path, exc))
    for src, exc in failures:
        unexpected.append(("①data_example", src, exc))
    baseline_note = (
        "与 121 基线一致"
        if len(files) == DATA_EXAMPLE_BASELINE
        else f"基线漂移：期望 {DATA_EXAMPLE_BASELINE}，实际 {len(files)}（数据目录变动，不算失败）"
    )
    return {
        "n": len(files), "ok": ok, "failures": failures,
        "baseline_note": baseline_note,
    }


# ---------------------------------------------------------------------------
# face 2: user real samples from config
# ---------------------------------------------------------------------------


def run_face2(unexpected: list[tuple[str, Path | str, Exception]]) -> dict:
    from md_analysis import config as mdconfig

    entries: list[tuple[str, str, Path]] = []  # (key, status, path)
    notes: list[str] = []
    parsed: list[tuple[Path, Exception]] = []
    ok = 0

    if not mdconfig.CONFIG_FILE.exists():
        notes.append(
            f"config 文件不存在（{mdconfig.CONFIG_FILE}）→ 本机未配置任何输入路径；"
            "不算失败，仅记录"
        )
        return {"entries": entries, "notes": notes, "n": 0, "ok": 0,
                "failures": [], "skipped": True}

    cfg = mdconfig.load_config()
    priority = [mdconfig.KEY_SP_INP_TEMPLATE_PATH, mdconfig.KEY_DP_SP_INP_TEMPLATE_PATH]
    scan_keys = priority + sorted(
        k for k, v in cfg.items()
        if k not in priority and isinstance(v, str)
    )
    seen_paths: set[Path] = set()
    for key in scan_keys:
        if key not in cfg:
            notes.append(f"`{key}`：未配置（不算失败）")
            continue
        val = cfg[key]
        if not isinstance(val, str) or not val:
            notes.append(f"`{key}`：值非路径字符串（{type(val).__name__}），跳过")
            continue
        name = Path(val).name.lower()
        known_key = key in priority
        if not known_key and not (name.endswith(".inp") or "restart" in name):
            continue  # unrelated config entry (vasp script path, numbers, ...)
        path = Path(val).expanduser()
        if path in seen_paths:
            entries.append((key, "duplicate（已由其他键检查）", path))
            continue
        seen_paths.add(path)
        if not path.is_file():
            entries.append((key, "缺失（文件不存在，不算失败）", path))
            continue
        _, exc = _parse_safe(path)
        if exc is None:
            ok += 1
            entries.append((key, "parse OK", path))
        else:
            entries.append((key, "parse 失败（意外）", path))
            parsed.append((path, exc))
            unexpected.append(("②user-config", path, exc))
    return {"entries": entries, "notes": notes, "n": len(seen_paths),
            "ok": ok, "failures": parsed, "skipped": False}


# ---------------------------------------------------------------------------
# face 3: packaged CP2K templates (importlib.resources)
# ---------------------------------------------------------------------------


def _diagnose_placeholder(text: str) -> str:
    """Supplementary diagnostic: parse with the placeholder line removed.

    This does NOT touch the preprocessing pipeline; it only confirms the
    plan §3.3 claim that deleting the placeholder line legalises the
    template (to be done in step 5).
    """
    lines = text.split("\n")
    hits = [i for i, line in enumerate(lines) if PLACEHOLDER_TOKEN in line]
    if len(hits) != 1:
        return f"占位符 `{PLACEHOLDER_TOKEN}` 出现 {len(hits)} 次（非 1），未做删行诊断"
    reduced = "\n".join(line for i, line in enumerate(lines) if i != hits[0])
    tree, exc = _parse_safe(reduced)
    if exc is None:
        return (
            f"删除第 {hits[0] + 1} 行占位符后 parse 成功"
            "（印证步骤 5『删行合法化』方案；本诊断未改任何管线）"
        )
    return (
        f"删除第 {hits[0] + 1} 行占位符后仍失败：{_err_class(exc)}"
        "（步骤 5 合法化时需人工复核）"
    )


def run_face3(unexpected: list[tuple[str, Path | str, Exception]]) -> dict:
    base = importlib.resources.files("md_analysis.scripts") / "template"
    templates = sorted(
        (p for p in base.iterdir() if p.name.endswith(".inp")),
        key=lambda p: p.name,
    )
    entries: list[tuple[str, str, str]] = []  # (name, status, detail)
    expected_failures = 0
    unexpected_failures: list[tuple[str, Exception]] = []
    diagnostics: list[str] = []
    for tpl in templates:
        try:
            text = tpl.read_text(encoding="utf-8")
        except Exception as exc:
            entries.append((tpl.name, "读取失败（意外）", _err_class(exc)))
            unexpected_failures.append((tpl.name, exc))
            unexpected.append(("③package-template", tpl.name, exc))
            continue
        tree, exc = _parse_safe(text)
        if exc is None:
            rt = _roundtrip_safe(tree)
            if rt is None:
                entries.append((tpl.name, "parse + generate roundtrip OK", ""))
            else:
                entries.append((tpl.name, "roundtrip 失败（意外）", _err_class(rt)))
                unexpected_failures.append((tpl.name, rt))
                unexpected.append(("③package-template", tpl.name, rt))
            continue
        if tpl.name == PLACEHOLDER_TEMPLATE and PLACEHOLDER_TOKEN in text:
            entries.append(
                (tpl.name, "已知预期失败（步骤 5 合法化后复跑）", _err_class(exc))
            )
            expected_failures += 1
            diagnostics.append(_diagnose_placeholder(text))
        else:
            entries.append((tpl.name, "parse 失败（意外）", _err_class(exc)))
            unexpected_failures.append((tpl.name, exc))
            unexpected.append(("③package-template", tpl.name, exc))
    return {"entries": entries, "n": len(templates),
            "expected_failures": expected_failures,
            "unexpected_failures": unexpected_failures,
            "diagnostics": diagnostics}


# ---------------------------------------------------------------------------
# face 4: 03_CuAg corpus
# ---------------------------------------------------------------------------


def run_face4(corpus_root: Path, unexpected: list[tuple[str, Path | str, Exception]]) -> dict:
    if not corpus_root.is_dir():
        return {"skipped": True,
                "note": f"语料根目录不存在：{corpus_root}（换机器场景，跳过不算失败）"}

    files = sorted(
        p for p in corpus_root.rglob("*")
        if p.is_file() and p.suffix.lower() == ".inp"
    )
    ok, failures = 0, []
    ok_trees: dict[Path, dict] = {}
    for path in files:
        tree, exc = _parse_safe(path)
        if exc is None:
            ok += 1
            ok_trees[path] = tree
        else:
            failures.append((path, exc))
    for src, exc in failures:
        unexpected.append(("④03_CuAg-parse", src, exc))

    # COMBINE_COLVAR files: full generate_inp roundtrip (built-in self-check)
    combine_files: list[Path] = []
    for path in files:
        try:
            text = path.read_text(encoding="utf-8")
        except Exception:
            continue  # unreadable files already failed parse above
        if "combine_colvar" in text.lower():
            combine_files.append(path)
    rt_failures: list[tuple[Path, Exception]] = []
    for path in combine_files:
        tree = ok_trees.get(path)
        if tree is None:
            rt_failures.append((path, RuntimeError("parse 阶段已失败，跳过 roundtrip")))
            continue
        exc = _roundtrip_safe(tree)
        if exc is not None:
            rt_failures.append((path, exc))
    for src, exc in rt_failures:
        unexpected.append(("④03_CuAg-roundtrip", src, exc))

    # CRLF sampling: deterministic spread over successfully parsed files
    crlf_results: list[tuple[Path, str]] = []
    pool = sorted(ok_trees)
    if pool:
        step = max(1, len(pool) // CRLF_SAMPLE_N)
        sample = pool[::step][:CRLF_SAMPLE_N]
        for path in sample:
            lf = path.read_text(encoding="utf-8").replace("\r\n", "\n")
            crlf = lf.replace("\n", "\r\n")
            tree_a, exc_a = _parse_safe(lf)
            tree_b, exc_b = _parse_safe(crlf)
            if exc_a is not None or exc_b is not None:
                crlf_results.append((path, f"失败：{_err_class(exc_a or exc_b)}"))
                unexpected.append(("④03_CuAg-crlf", path, exc_a or exc_b))
            elif json.dumps(tree_a, sort_keys=True) != json.dumps(tree_b, sort_keys=True):
                crlf_results.append((path, "LF/CRLF 树不一致（意外）"))
                unexpected.append(("④03_CuAg-crlf", path, RuntimeError("tree mismatch")))
            else:
                crlf_results.append((path, "LF 与 CRLF 树全等"))
    combine_note = (
        "与计划登记的 66 个一致"
        if len(combine_files) == EXPECTED_COMBINE_COLVAR
        else f"计划登记 {EXPECTED_COMBINE_COLVAR}，实测 {len(combine_files)}（语料变动）"
    )
    return {
        "skipped": False, "n": len(files), "ok": ok, "failures": failures,
        "combine_n": len(combine_files), "combine_note": combine_note,
        "rt_failures": rt_failures, "crlf_results": crlf_results,
    }


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def _failure_table(groups: dict[str, list[Path | str]], max_examples: int = 3) -> list[str]:
    lines = []
    for cls, paths in groups.items():
        lines.append(f"- `{cls}` — {len(paths)} 个")
        for p in paths[:max_examples]:
            lines.append(f"  - {p}")
        if len(paths) > max_examples:
            lines.append(f"  - …其余 {len(paths) - max_examples} 个见失败清单")
    return lines


def write_report(report_path: Path, r1, r2, r3, r4, unexpected, elapsed: float) -> None:
    import md_analysis

    lines: list[str] = []
    lines.append("# CP2K inp 输入预检报告（重构步骤 4）")
    lines.append("")
    lines.append(f"- 运行时间：{datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"- 分支/commit：{_git_info()}")
    lines.append(f"- Python：{sys.executable}")
    lines.append(f"- md_analysis：{md_analysis.__file__}")
    lines.append(f"- 依赖：{_dep_versions()}")
    lines.append(f"- 工具：`tools/precheck_cp2k_inp.py`（一次性诊断，不进 CI）")
    lines.append(f"- 总耗时：{elapsed:.1f}s")
    lines.append("")

    verdict = "**零意外失败**" if not unexpected else f"**{len(unexpected)} 个意外失败**（明细见下）"
    lines.append("## 总结论")
    lines.append("")
    lines.append(verdict)
    lines.append("")
    lines.append("| 输入面 | N | parse OK | 失败 | 已知预期失败 | 意外失败 |")
    lines.append("|---|---|---|---|---|---|")
    n2 = r2["n"]
    f2 = len(r2["failures"])
    lines.append(f"| ① data_example | {r1['n']} | {r1['ok']} | {len(r1['failures'])} | 0 | {len(r1['failures'])} |")
    if r2.get("skipped"):
        lines.append(f"| ② 用户真实样本 | 未配置 | 0 | 0 | 0 | 0 |")
    else:
        lines.append(f"| ② 用户真实样本（config） | {n2} | {r2['ok']} | {f2} | 0 | {f2} |")
    u3 = len(r3["unexpected_failures"])
    lines.append(
        f"| ③ 包内 CP2K 模板 | {r3['n']} | {r3['n'] - r3['expected_failures'] - u3} "
        f"| {r3['expected_failures'] + u3} | {r3['expected_failures']} | {u3} |"
    )
    if r4.get("skipped"):
        lines.append("| ④ 03_CuAg 语料 | 跳过 | - | - | - | - |")
    else:
        u4 = len(r4["failures"]) + len(r4["rt_failures"]) + sum(
            1 for _, s in r4["crlf_results"] if "全等" not in s
        )
        lines.append(
            f"| ④ 03_CuAg 语料 | {r4['n']} | {r4['ok']} | {u4} | 0 | {u4} |"
        )
    lines.append("")

    lines.append("## ① data_example（基线自校验）")
    lines.append("")
    lines.append(f"- 文件数 {r1['n']}（{r1['baseline_note']}）：parse {r1['ok']} 通过 / {len(r1['failures'])} 失败")
    if r1["failures"]:
        lines.extend(_failure_table(_group(r1["failures"])))
    lines.append("")

    lines.append("## ② 用户真实样本（config 指向清单）")
    lines.append("")
    for note in r2["notes"]:
        lines.append(f"- {note}")
    if r2["entries"]:
        lines.append("")
        lines.append("| config 键 | 状态 | 路径 |")
        lines.append("|---|---|---|")
        for key, status, path in r2["entries"]:
            lines.append(f"| `{key}` | {status} | `{path}` |")
    lines.append("")

    lines.append("## ③ 包内 CP2K 模板（importlib.resources 枚举）")
    lines.append("")
    lines.append(f"- 模板数：{r3['n']}")
    if r3["entries"]:
        lines.append("")
        lines.append("| 模板 | 状态 | 错误类别/说明 |")
        lines.append("|---|---|---|")
        for name, status, detail in r3["entries"]:
            lines.append(f"| `{name}` | {status} | {detail or '—'} |")
    for diag in r3["diagnostics"]:
        lines.append(f"- 补充诊断：{diag}")
    lines.append(
        "- 边界确认：`{{...}}` 占位符不是预处理对象，未为消化它扩展 W1-W3 管线；"
        "合法化（删占位符行）属步骤 5，届时复跑本面验收"
    )
    lines.append("")

    lines.append("## ④ 03_CuAg 语料")
    lines.append("")
    if r4.get("skipped"):
        lines.append(f"- {r4['note']}")
    else:
        lines.append(
            f"- parse 复跑：{r4['n']} 个文件，{r4['ok']} 通过 / {len(r4['failures'])} 失败"
        )
        if r4["failures"]:
            lines.extend(_failure_table(_group(r4["failures"])))
        lines.append("")
        lines.append(
            f"- COMBINE_COLVAR 全量 generate_inp roundtrip（含 inp.py 内建自校验）："
            f"{r4['combine_n']} 个（{r4['combine_note']}），"
            f"{r4['combine_n'] - len(r4['rt_failures'])} 通过 / {len(r4['rt_failures'])} 失败"
        )
        if r4["rt_failures"]:
            lines.extend(_failure_table(_group(r4["rt_failures"])))
        lines.append("")
        lines.append(f"- CRLF 抽样（{CRLF_SAMPLE_N} 个文件，LF→CRLF 后树全等性）：")
        for path, status in r4["crlf_results"]:
            lines.append(f"  - `{path}` — {status}")
    lines.append("")

    lines.append("## 意外失败明细")
    lines.append("")
    if not unexpected:
        lines.append("无——全部结果均在预期内（③ 的模板失败为已登记的预期失败）。")
    else:
        grouped: dict[str, list] = {}
        for face, src, exc in unexpected:
            grouped.setdefault(face, []).append((src, exc))
        for face, items in grouped.items():
            lines.append(f"### {face}（{len(items)}）")
            lines.extend(_failure_table(_group(items)))
            lines.append("")
    lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description="CP2K input pre-check (refactor step 4)")
    parser.add_argument(
        "--corpus", type=Path, default=CORPUS_ROOT,
        help=f"03_CuAg corpus root (default: {CORPUS_ROOT})",
    )
    parser.add_argument(
        "--report", type=Path, default=REPORT_PATH,
        help=f"report output path (default: {REPORT_PATH})",
    )
    args = parser.parse_args()

    t0 = time.perf_counter()
    unexpected: list[tuple[str, Path | str, Exception]] = []

    print("[1/4] data_example (baseline self-check) ...", flush=True)
    r1 = run_face1(unexpected)
    print(f"      {r1['ok']}/{r1['n']} parse OK, {len(r1['failures'])} failed", flush=True)

    print("[2/4] user real samples (config) ...", flush=True)
    r2 = run_face2(unexpected)
    if r2.get("skipped"):
        print("      config not present -> reported, not a failure", flush=True)
    else:
        print(f"      {r2['ok']}/{r2['n']} parse OK, {len(r2['failures'])} failed", flush=True)

    print("[3/4] packaged CP2K templates ...", flush=True)
    r3 = run_face3(unexpected)
    print(
        f"      {r3['n']} templates, {r3['expected_failures']} known-expected "
        f"failures, {len(r3['unexpected_failures'])} unexpected", flush=True,
    )

    print(f"[4/4] 03_CuAg corpus ({args.corpus}) ... (directory walk may take a while)", flush=True)
    t4 = time.perf_counter()
    r4 = run_face4(args.corpus, unexpected)
    if r4.get("skipped"):
        print(f"      skipped: {r4['note']}", flush=True)
    else:
        print(
            f"      parse {r4['ok']}/{r4['n']}, roundtrip "
            f"{r4['combine_n'] - len(r4['rt_failures'])}/{r4['combine_n']}, "
            f"CRLF sample {sum(1 for _, s in r4['crlf_results'] if '全等' in s)}"
            f"/{len(r4['crlf_results'])}  [{time.perf_counter() - t4:.1f}s]", flush=True,
        )

    elapsed = time.perf_counter() - t0
    write_report(args.report, r1, r2, r3, r4, unexpected, elapsed)
    print(f"report written to {args.report} [{elapsed:.1f}s total]", flush=True)

    if unexpected:
        print(f"UNEXPECTED FAILURES: {len(unexpected)} (see report)", file=sys.stderr)
        return 1
    print("PRECHECK OK: zero unexpected failures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
