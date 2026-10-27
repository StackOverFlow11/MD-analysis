"""Parse and generate CP2K input (``.inp``) files via cp2k-input-tools.

This is the single place in the package that touches ``cp2k-input-tools``
(cit) and ``pint``; both imports are deferred to first use so that package
startup does not pay their import cost (~0.4 s + ~0.2 s).

Tree conventions (cit ``CP2KInputParserSimplified`` output)
-----------------------------------------------------------
- Section and keyword keys are the **lowercase canonical XML names**
  (``NAME[@type='default']``).  Input aliases are normalised on parse:
  ``PROJECT`` becomes ``project_name``.  Always read/write canonical names;
  using an alias as a key when generating is not rejected but produces a
  duplicate keyword which fails re-parse (``NameRepetitionError``).
- A repeated section whose instances carry unique string parameters
  (e.g. ``&KIND O``) becomes a **dict keyed by the parameter**; a repeated
  section without parameters (e.g. nested ``&COLVAR`` in ``&COMBINE_COLVAR``,
  multiple ``&COLLECTIVE`` blocks) becomes a **list of dicts**.
- ``"_"`` is the section parameter key; ``"*"`` holds default-keyword lines
  (e.g. coordinates in ``&COORD``).
- Multi-value keywords are tuples; repeated keywords are lists of values;
  keyword values are native ``int``/``float``/``bool``/``str``.

Units
-----
Keyword values land on the XML *default unit* of each keyword.  For
``internal_cp2k`` keywords (``TARGET`` / ``TARGET_GROWTH`` in
``&COLLECTIVE``) that is CP2K's internal unit system: length bohr, angle
radian, time atomic unit of time (1 fs = 41.341374 aut), energy hartree --
the same convention as ``.restart`` files, so inp trees and restart trees
can be compared directly.  Explicit ``[unit]`` tokens on those keywords are
rejected by cit; the preprocessing pipeline below converts them to internal
units and strips the token before parsing.

Preprocessing pipeline (parse side, in order)
---------------------------------------------
``CRLF normalisation -> G1 -> W3 -> W1 -> W2 -> comment cleaning``

- CRLF: ``\\r\\n`` -> ``\\n`` (first step).
- G1: ``@`` preprocessor directives at line start are **rejected** with
  :class:`InpParseError` (cit would silently mis-handle e.g. bracketed
  ``@IF`` conditions); run ``cp2k --preprocess`` and feed the expanded file.
- W3: comma-separated values on ``LIST`` lines become space-separated
  (case-insensitive keyword match, only the value part is touched).
- W1: integer ranges ``2..4`` expand to ``2 3 4`` (cit issue #72); works
  across ``\\`` line continuations.
- W2: explicit units on ``TARGET`` / ``TARGET_GROWTH`` are converted to CP2K
  internal units and the ``[unit]`` token is stripped (atom-wise two-way
  mapping through pint; ``fs``/``ps`` map to ``atomic_unit_of_time``).
- Comment cleaning: quote-aware truncation of inline comments that contain
  quote characters (prevents cit tokenizer crash, upstream issue #111).
  Whole-line comments and comment-free lines are left byte-identical; quote
  characters inside quoted values are preserved.

Generation
----------
``generate_inp`` serialises a (possibly edited) tree back to CP2K syntax.
cit's generator emits same-level sections in reverse (stack/LIFO); to keep
the *semantic* order of repeated-section lists, ``_stabilize_tree``
pre-reverses every list-of-dicts (the list form of a repeated section)
before generation, and a fail-fast self-check re-parses the generated text
and compares it against the caller's original tree (dict key order
insensitive, list order sensitive, tuple/list normalised, int/float
distinct).  Any mismatch raises :class:`InpGenerateError` with the first
differing path.  Comments, blank lines, keyword letter case and the textual
order of *non-repeating* sections are not preserved (accepted readability
loss); repeated-section containers must be lists (tuples are rejected by
the generator).
"""

from __future__ import annotations

import io
import json
import logging
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

from ....exceptions import MDAnalysisError

__all__ = [
    "InpParseError",
    "InpGenerateError",
    "parse_inp",
    "generate_inp",
    "get_section",
    "get_keyword",
    "set_keyword",
    "iter_collectives",
    "set_cell_abc",
]

logger = logging.getLogger(__name__)


class InpParseError(MDAnalysisError):
    """Raised when a CP2K input text cannot be parsed (syntax or content)."""


class InpGenerateError(MDAnalysisError):
    """Raised when a CP2K input tree cannot be generated faithfully."""


# ---------------------------------------------------------------------------
# Lazy singletons (cit and pint are imported on first use only)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _get_parser(xmlspec: str | os.PathLike | None = None):
    """Return a cached cit simplified parser (keyed by *xmlspec*)."""
    from cp2k_input_tools.parser import CP2KInputParserSimplified

    return CP2KInputParserSimplified(xmlspec=xmlspec)


@lru_cache(maxsize=None)
def _get_generator(xmlspec: str | os.PathLike | None = None):
    """Return a cached cit generator (keyed by *xmlspec*).

    Unlike the parser, the generator does not fall back to the bundled XML
    when handed ``xmlspec=None``, so the default is applied here.
    """
    from cp2k_input_tools.generator import CP2KInputGenerator

    if xmlspec is None:
        return CP2KInputGenerator()
    return CP2KInputGenerator(xmlspec=xmlspec)


@lru_cache(maxsize=None)
def _get_ureg():
    """Return a cached pint unit registry."""
    from pint import UnitRegistry

    return UnitRegistry()


# ---------------------------------------------------------------------------
# Preprocessing pipeline: CRLF -> G1 -> W3 -> W1 -> W2 -> comment cleaning
# ---------------------------------------------------------------------------

_COMMENT_CHARS = ("#", "!")

# W1: integer range expansion (cit issue #72); look-around guards avoid
# matching inside longer tokens (e.g. floating point numbers, words).
_RANGE_RE = re.compile(r"(?<![\w.])(\d+)\.\.(\d+)(?![\w.])")

# W3: `LIST a,b,c` lines; only the value part after the keyword is rewritten.
_LIST_LINE_RE = re.compile(r"(?im)^(\s*list\s+)(\S.*)$")

# W2: explicit units on internal_cp2k keywords (cit rejects the unit token).
_W2_KEYWORDS = ("TARGET_GROWTH", "TARGET")
_UNIT_KW_RE = re.compile(
    r"(?im)^(\s*)(%s)\s+\[([^\]\n]+)\]\s+(.+)$" % "|".join(_W2_KEYWORDS)
)

# Atom-wise two-way unit mapping for W2:
#   source token -> unambiguous pint full name (fs must not become
#   femtosiemens, which pint resolves `fs` to with case_sensitive=False)
#   source atom  -> CP2K internal unit atom (fs/ps -> atomic_unit_of_time,
#   v7 plan section 0: 1 fs = 41.34137333518245 aut)
_SRC_ATOM_MAP = {
    "angstrom": "angstrom",
    "bohr": "bohr",
    "nm": "nanometer",
    "deg": "degree",
    "rad": "radian",
    "fs": "femtosecond",
    "ps": "picosecond",
    "hartree": "hartree",
    "ev": "electron_volt",
    "k": "kelvin",
}
_UNIT_ATOM_MAP = {
    "angstrom": "bohr",
    "bohr": "bohr",
    "nm": "bohr",
    "deg": "radian",
    "rad": "radian",
    "fs": "atomic_unit_of_time",
    "ps": "atomic_unit_of_time",
    "hartree": "hartree",
    "ev": "hartree",
    "k": "kelvin",
}
_UNIT_ATOM_RE = re.compile(r"[a-zA-Z]+")


def _gate_at_directives(text: str) -> None:
    """G1: reject CP2K ``@`` preprocessor directives instead of rewriting.

    A line whose first non-blank character is ``@`` is a directive for any
    CP2K preprocessor (this mirrors cit's own dispatch criterion, so there
    are no structural false positives; continuation lines are tracked).
    Whole-line comments and blank lines are skipped.
    """
    continuation = False
    for lineno, line in enumerate(text.split("\n"), start=1):
        if not continuation:
            stripped = line.lstrip()
            if not stripped or stripped.startswith(_COMMENT_CHARS):
                continuation = line.rstrip().endswith("\\")
                continue
            if stripped.startswith("@"):
                raise InpParseError(
                    f"line {lineno}: CP2K preprocessor directive "
                    f"{stripped[:60]!r} is not supported; run "
                    "`cp2k --preprocess` on the file and parse the "
                    "expanded output instead"
                )
        continuation = line.rstrip().endswith("\\")


def _normalize_list_commas(text: str) -> str:
    """W3: `LIST 1,2,7..10` -> `LIST 1 2 7..10` (value area only)."""

    def _sub(match: re.Match) -> str:
        return match.group(1) + match.group(2).replace(",", " ")

    return _LIST_LINE_RE.sub(_sub, text)


def _expand_ranges(text: str) -> str:
    """W1: `2..4` -> `2 3 4` (cit issue #72)."""

    def _sub(match: re.Match) -> str:
        start, end = int(match.group(1)), int(match.group(2))
        if end >= start and end - start < 10_000:
            return " ".join(str(i) for i in range(start, end + 1))
        return match.group(0)  # not a range; let cit judge the line

    return _RANGE_RE.sub(_sub, text)


def _convert_internal_units(text: str) -> str:
    """W2: convert explicit units on TARGET/TARGET_GROWTH to internal units.

    ``TARGET [angstrom] 1.0`` becomes ``TARGET 1.8897261246222794``.  Unit
    atoms are translated both to unambiguous pint names (source) and to the
    CP2K internal unit (destination); pint does the numeric conversion.
    Lines that cannot be converted are left untouched -- cit then rejects
    the unit token with its own line-numbered error (fail fast).
    """
    ureg = _get_ureg()

    def _atoms(unit_str: str, mapping: dict[str, str]) -> str:
        def _sub(match: re.Match) -> str:
            atom = match.group(0).lower()
            if atom not in mapping:
                raise ValueError(f"unknown unit atom: {atom!r}")
            return mapping[atom]

        return _UNIT_ATOM_RE.sub(_sub, unit_str.replace("^", "**"))

    def _sub(match: re.Match) -> str:
        indent, kw, unit, values = match.groups()
        try:
            src_unit = ureg.parse_expression(_atoms(unit, _SRC_ATOM_MAP))
            dst_unit = ureg.parse_expression(_atoms(unit, _UNIT_ATOM_MAP))
            converted = []
            for token in values.split():
                try:
                    value = float(token.replace("D", "E").replace("d", "e"))
                except ValueError:
                    converted.append(token)  # e.g. trailing comment words
                    continue
                converted.append(repr((value * src_unit).to(dst_unit).magnitude))
            return f"{indent}{kw.upper()} {' '.join(converted)}"
        except Exception:
            return match.group(0)

    return _UNIT_KW_RE.sub(_sub, text)


def _clean_comments(text: str) -> str:
    """Truncate inline comments containing quotes (cit issue #111 guard).

    Quote-aware line scan: a comment starting outside quotes is truncated
    at the comment character only when the comment part itself contains a
    quote character (otherwise cit handles the comment natively).  Whole
    line comments are never touched (cit skips them before tokenizing);
    quotes inside quoted values are preserved.  Lines that continue a
    ``\\`` continuation are treated as inline (cit joins them first).
    """
    lines = text.split("\n")
    out: list[str] = []
    continuation = False
    for line in lines:
        if not continuation:
            stripped = line.lstrip()
            if stripped.startswith(_COMMENT_CHARS):
                out.append(line)
                continuation = line.rstrip().endswith("\\")
                continue
        quote: str | None = None
        prev = ""
        cut: int | None = None
        for i, ch in enumerate(line):
            if quote is not None:
                if ch == quote and prev != "\\":
                    quote = None
            elif ch in "'\"":
                quote = ch
            elif ch in _COMMENT_CHARS:
                cut = i
                break
            prev = ch
        if cut is not None and ("'" in line[cut:] or '"' in line[cut:]):
            new_line = line[:cut].rstrip()
            if line.rstrip().endswith("\\"):
                new_line += " \\"
            line = new_line
        out.append(line)
        continuation = line.rstrip().endswith("\\")
    return "\n".join(out)


def _preprocess(text: str) -> str:
    """Apply the full parse-side pipeline (order pinned, see module doc)."""
    text = text.replace("\r\n", "\n")
    _gate_at_directives(text)
    text = _normalize_list_commas(text)
    text = _expand_ranges(text)
    text = _convert_internal_units(text)
    text = _clean_comments(text)
    return text


# ---------------------------------------------------------------------------
# parse_inp / generate_inp
# ---------------------------------------------------------------------------


def _load_input(text_or_path: Any) -> tuple[str, str]:
    """Dual dispatch: PathLike -> path; str with newline -> text; else path."""
    if isinstance(text_or_path, os.PathLike):
        path = Path(text_or_path)
        return path.read_text(encoding="utf-8"), str(path)
    if isinstance(text_or_path, str):
        if "\n" in text_or_path:
            return text_or_path, "<input text>"
        path = Path(text_or_path)
        try:
            return path.read_text(encoding="utf-8"), str(path)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"parse_inp: input {text_or_path[:80]!r} contains no newline, "
                "so it was treated as a file path, but no such file exists; "
                "pass the file contents (text) or an existing path"
            ) from None
    raise TypeError(
        "parse_inp expects input text (str) or a path (str/PathLike), "
        f"got {type(text_or_path).__name__}"
    )


def _format_parse_error(src: str, exc: Exception) -> str:
    """Attach cit's line context (if any) to the wrapped error message."""
    ctx = exc.args[-1] if exc.args else None
    if hasattr(ctx, "linenr") and ctx.linenr is not None:
        info = f"line {ctx.linenr}"
        if getattr(ctx, "line", None):
            info += f" ({str(ctx.line).strip()[:60]!r})"
        return f"cannot parse CP2K input {src}: {info}: {exc}"
    return f"cannot parse CP2K input {src}: {exc}"


def parse_inp(text_or_path: Any, *, xmlspec: str | os.PathLike | None = None) -> dict:
    """Parse a CP2K input file (or text) into a cit simplified tree.

    Parameters
    ----------
    text_or_path : str, os.PathLike
        ``os.PathLike`` is read as a file path (utf-8).  A ``str`` holding
        a newline is parsed as input text; a ``str`` without newlines is
        treated as a path (``FileNotFoundError`` with a hint if missing).
    xmlspec : str or os.PathLike, optional
        Path to a ``cp2k_input.xml`` produced by ``cp2k --xml``; defaults
        to the CP2K version bundled with cp2k-input-tools.

    Returns
    -------
    dict
        The simplified option tree (see module docstring for conventions).

    Raises
    ------
    InpParseError
        On invalid syntax/content (cit errors are wrapped with line info)
        and on unsupported ``@`` preprocessor directives.
    FileNotFoundError, TypeError
        On unusable input arguments.
    """
    text, src = _load_input(text_or_path)
    text = _preprocess(text)
    parser = _get_parser(xmlspec)
    try:
        return parser.parse(io.StringIO(text))
    except Exception as exc:
        raise InpParseError(_format_parse_error(src, exc)) from exc


def _stabilize_tree(node: Any) -> Any:
    """Pre-reverse repeated-section lists to compensate generator LIFO.

    cit's generator emits sections from a stack, so lists of sections come
    out reversed; reversing every list whose elements are all dicts (the
    repeated-section list form; keyword value lists never contain dicts)
    restores the original textual order after re-parsing.  Dicts and tuples
    are rebuilt unchanged (no reversal).
    """
    if isinstance(node, dict):
        return {key: _stabilize_tree(value) for key, value in node.items()}
    if isinstance(node, list):
        items = [_stabilize_tree(value) for value in node]
        if items and all(isinstance(item, dict) for item in items):
            items.reverse()
        return items
    if isinstance(node, tuple):
        return tuple(_stabilize_tree(value) for value in node)
    return node


def _canonical_json(tree: Any) -> str:
    """Order-insensitive canonical form: dict keys sorted, tuple == list,
    int != float (``1`` and ``1.0`` serialise differently in JSON)."""
    return json.dumps(tree, sort_keys=True)


def _first_difference(a: Any, b: Any, path: str = "tree") -> str | None:
    """Return a human-readable path of the first tree difference, or None."""
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b), key=str):
            if key not in a:
                return f"{path}.{key}: missing in source tree"
            if key not in b:
                return f"{path}.{key}: missing in re-parsed tree"
            diff = _first_difference(a[key], b[key], f"{path}.{key}")
            if diff is not None:
                return diff
        return None
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return f"{path}: length {len(a)} != {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            diff = _first_difference(x, y, f"{path}[{i}]")
            if diff is not None:
                return diff
        return None
    if _canonical_json(a) != _canonical_json(b):
        return f"{path}: {a!r} != {b!r}"
    return None


def generate_inp(tree: dict, *, xmlspec: str | os.PathLike | None = None) -> str:
    """Serialise a cit-style option tree back to CP2K input text.

    The tree is stabilised (repeated-section list order, see
    ``_stabilize_tree``), generated, then **re-parsed and compared against
    the caller's tree** as a fail-fast self-check: dict key order is
    ignored, list order and int/float distinctions are enforced.  Any
    mismatch raises :class:`InpGenerateError` (a generate failure is a
    programming error and is never translated into workflow errors).

    Repeated-section containers must be lists; tuples are rejected by the
    generator (wrapped into ``InpGenerateError``).
    """
    generator = _get_generator(xmlspec)
    try:
        lines = list(generator.line_iter(_stabilize_tree(tree)))
    except Exception as exc:
        raise InpGenerateError(f"cannot generate CP2K input from tree: {exc}") from exc
    text = "\n".join(lines) + "\n"
    try:
        reparsed = parse_inp(text, xmlspec=xmlspec)
    except InpParseError as exc:
        raise InpGenerateError(
            f"generated CP2K input does not re-parse (bug in tree or xmlspec): {exc}"
        ) from exc
    diff = _first_difference(tree, reparsed)
    if diff is not None:
        raise InpGenerateError(
            "generated CP2K input does not round-trip to the source tree "
            f"(first difference at {diff})"
        )
    return text


# ---------------------------------------------------------------------------
# Tree path primitives (format-level, no workflow semantics)
# ---------------------------------------------------------------------------

_PATH_PART_RE = re.compile(r"\s*([A-Za-z_][\w\-]*)\s*(?:\[\s*(\d+)\s*\])?\s*")


def _parse_path(path: str) -> list[tuple[str, int | None]]:
    """Split ``'force_eval.subsys.kind[0]'`` into ``(name, index)`` tokens."""
    tokens: list[tuple[str, int | None]] = []
    for part in path.split("."):
        match = _PATH_PART_RE.fullmatch(part)
        if match is None or not match.group(1):
            raise ValueError(f"invalid tree path: {path!r} (bad segment {part!r})")
        tokens.append(
            (match.group(1), int(match.group(2)) if match.group(2) is not None else None)
        )
    if not tokens:
        raise ValueError(f"invalid tree path: {path!r} (empty)")
    return tokens


def _walk_tokens(tree: Any, tokens: list[tuple[str, int | None]], path: str) -> Any:
    """Walk *tree* along *tokens*; descend list sections by their ``_`` param.

    Tolerates both repeated-section forms: dict-of-dicts (keyed by section
    parameter, e.g. ``kind.Pt``) and list-of-dicts (matched via the ``_``
    parameter key).  Raises ``KeyError`` with the reached path on failure.
    """
    node = tree
    consumed: list[str] = []
    for name, index in tokens:
        consumed.append(name)
        if isinstance(node, list):
            matches = [
                item
                for item in node
                if isinstance(item, dict) and item.get("_") == name
            ]
            if not matches:
                raise KeyError(f"{path!r}: no section {name!r} at {'.'.join(consumed)!r}")
            if len(matches) > 1:
                raise KeyError(
                    f"{path!r}: ambiguous section {name!r} at {'.'.join(consumed)!r}"
                )
            node = matches[0]
        elif isinstance(node, dict):
            if name not in node:
                raise KeyError(f"{path!r}: no key {name!r} at {'.'.join(consumed)!r}")
            node = node[name]
        else:
            raise KeyError(
                f"{path!r}: cannot descend into {name!r} at {'.'.join(consumed)!r}"
            )
        if index is not None:
            if not isinstance(node, (list, tuple)) or index >= len(node):
                raise KeyError(f"{path!r}: bad index [{index}] at {'.'.join(consumed)!r}")
            node = node[index]
    return node


def get_section(tree: dict, path: str) -> Any:
    """Return the section container (dict or list) at dotted *path*.

    Path syntax: ``force_eval.subsys.cell``, with optional ``[n]`` indices
    for list-form sections and name-based descent into list sections via
    their ``_`` parameter (``force_eval.subsys.kind.Pt`` works for both the
    dict and list form).  Raises ``KeyError`` if any segment is missing.
    """
    return _walk_tokens(tree, _parse_path(path), path)


def get_keyword(tree: dict, path: str) -> Any:
    """Return the keyword value at dotted *path* (same walk as get_section)."""
    return _walk_tokens(tree, _parse_path(path), path)


def set_keyword(tree: dict, path: str, value: Any) -> None:
    """Set (or create) the keyword *value* at dotted *path*, in place.

    The last path component is the keyword name; all parent sections must
    already exist (creating sections is workflow logic, not a format
    primitive).  ``path`` must not end in an ``[n]`` index.
    """
    tokens = _parse_path(path)
    name, index = tokens[-1]
    if len(tokens) < 2:
        raise ValueError(f"path must name a parent section and a keyword: {path!r}")
    if index is not None:
        raise ValueError(f"keyword path must not end in an index: {path!r}")
    parent = _walk_tokens(tree, tokens[:-1], path)
    if not isinstance(parent, dict):
        raise KeyError(f"{path!r}: parent of keyword {name!r} is not a section dict")
    parent[name] = value


def iter_collectives(tree: dict) -> Iterator[tuple[int, dict]]:
    """Yield ``(colvar_id, collective_block)`` for every &COLLECTIVE block.

    Works for both tree forms (single dict or list of dicts).  Raises
    ``KeyError`` when ``motion.constraint`` is missing.
    """
    blocks = get_section(tree, "motion.constraint.collective")
    if isinstance(blocks, dict):
        blocks = [blocks]
    for block in blocks:
        if "colvar" not in block:
            raise KeyError("motion.constraint.collective block without COLVAR id")
        yield block["colvar"], block


def set_cell_abc(tree: dict, a: float, b: float, c: float) -> None:
    """Write ``force_eval.subsys.cell.abc`` = ``(a, b, c)`` in place.

    Templates using the A/B/C vector form (or missing the ``abc`` keyword /
    ``cell`` section) are **not** modified: a warning is logged and the
    call is a no-op (same semantics as the legacy text-based helper).
    """
    try:
        cell = get_section(tree, "force_eval.subsys.cell")
        if not isinstance(cell, dict) or "abc" not in cell:
            raise KeyError("abc")
    except KeyError:
        logger.warning("No &CELL ABC line found in template; cell not updated")
        return
    cell["abc"] = (float(a), float(b), float(c))
