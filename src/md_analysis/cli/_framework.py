"""Core framework: MenuNode, MenuGroup, MenuCommand, lazy_import."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

from ..exceptions import MDAnalysisError
from ._prompt import _read, prompt_bool

if TYPE_CHECKING:
    from ._params import ParamCollector

logger = logging.getLogger(__name__)


def lazy_import(module_path: str, name: str):
    """Lazily import a callable from a dotted module path."""
    from importlib import import_module

    return getattr(import_module(module_path), name)


class MenuNode(ABC):
    """Abstract base for all menu tree nodes."""

    output_name: str = ""
    """Path segment this node contributes to the output directory."""

    def __init__(self, code: str, label: str) -> None:
        self.code = code
        self.label = label
        self.parent: MenuNode | None = None

    @abstractmethod
    def run(self) -> None: ...


class MenuGroup(MenuNode):
    """Non-leaf node: displays a menu and dispatches to children."""

    def __init__(self, code: str, label: str, *,
                 output_name: str = "") -> None:
        super().__init__(code, label)
        self.output_name = output_name
        self.children: list[MenuNode | str] = []
        self._flat_index: dict[str, MenuCommand] | None = None

    def add(self, *nodes: MenuNode | str) -> None:
        for node in nodes:
            if isinstance(node, MenuNode):
                node.parent = self
            self.children.append(node)

    def run(self) -> None:
        while True:
            print(self._render_menu())
            choice = _read(" Input: ").strip()

            if choice == "0":
                return

            # 1) Match direct child by code
            for child in self.children:
                if isinstance(child, MenuNode) and child.code == choice:
                    child.run()
                    break
            else:
                # 2) Fallback: flat index for leaf shortcut (root only)
                if self._flat_index and choice in self._flat_index:
                    self._flat_index[choice].run()
                else:
                    print(f"\n Invalid choice: {choice!r}")

    def build_flat_index(self) -> None:
        """Recursively index all descendant MenuCommand nodes by code.
        Call once on root after tree is fully assembled.

        The flat index is shared with all descendant ``MenuGroup`` nodes
        so that leaf-code shortcuts work from any level in the menu tree.
        """
        self._flat_index = {}
        self._collect_commands(self._flat_index)
        self._propagate_flat_index(self._flat_index)

    def _propagate_flat_index(self, index: dict[str, MenuCommand]) -> None:
        """Share *index* with all descendant MenuGroup nodes."""
        for child in self.children:
            if isinstance(child, MenuGroup):
                child._flat_index = index
                child._propagate_flat_index(index)

    def _collect_commands(self, index: dict[str, MenuCommand]) -> None:
        for child in self.children:
            if isinstance(child, MenuCommand):
                index[child.code] = child
            elif isinstance(child, MenuGroup):
                child._collect_commands(index)

    def _render_menu(self) -> str:
        lines = [f"\n ---------- {self.label} ----------\n"]
        for child in self.children:
            if isinstance(child, str):
                if child:
                    lines.append(f"\n --- {child} ---")
                else:
                    lines.append("")
            else:
                lines.append(f" {child.code}) {child.label}")
        lines.append("\n   0) Back / Exit\n")
        return "\n".join(lines)


class MenuCommand(MenuNode):
    """Leaf node: collects parameters and executes an analysis command."""

    params: tuple[ParamCollector, ...] = ()
    advanced_params: tuple[ParamCollector, ...] = ()

    @property
    def output_subdir(self) -> str:
        """Derive output sub-directory by walking up the menu tree.

        Each node's ``output_name`` contributes one path segment.
        """
        parts: list[str] = []
        node: MenuNode | None = self
        while node is not None:
            if node.output_name:
                parts.append(node.output_name)
            node = node.parent
        parts.reverse()
        return "/".join(parts) if parts else ""

    def run(self) -> None:
        try:
            ctx = self._collect_all_params()
            subdir = self.output_subdir
            if subdir:
                outdir = Path(ctx.get("outdir", "analysis")) / subdir
                outdir.mkdir(parents=True, exist_ok=True)
                ctx["outdir_resolved"] = outdir
            self.execute(ctx)
        except (MDAnalysisError, FileNotFoundError, ValueError, RuntimeError) as exc:
            print(f"\n  Error: {exc}")
        except Exception as exc:
            logger.error("Unexpected error in %s: %s", self.label, exc, exc_info=True)
            print(f"\n  Unexpected error ({type(exc).__name__}): {exc}")

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        for p in self.params:
            p.collect(ctx)
        if self.advanced_params:
            if prompt_bool("Modify advanced parameters?", default=False):
                for p in self.advanced_params:
                    p.collect(ctx)
            else:
                for p in self.advanced_params:
                    p.apply_default(ctx)
        return ctx

    @abstractmethod
    def execute(self, ctx: dict) -> None:
        """Subclasses implement pure business logic here."""
        ...
