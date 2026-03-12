"""Enhanced Sampling command classes (501-502)."""

from __future__ import annotations

import re
from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import K
from ._prompt import prompt_int, prompt_str


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------

def _discover_restart_file(workdir: Path) -> str | None:
    """Glob for ``*.restart``, excluding ``_<digits>.restart`` checkpoints."""
    checkpoint_re = re.compile(r"_\d+\.restart$")
    candidates = [
        p for p in sorted(workdir.glob("*.restart"))
        if not checkpoint_re.search(p.name)
    ]
    if len(candidates) == 1:
        return str(candidates[0])
    return None


def _discover_log_file(workdir: Path) -> str | None:
    """Glob for ``*.LagrangeMultLog``."""
    candidates = sorted(workdir.glob("*.LagrangeMultLog"))
    if len(candidates) == 1:
        return str(candidates[0])
    return None


def _print_sg_info(restart_path: str, log_path: str) -> None:
    """Parse and display trajectory metadata."""
    ColvarMDInfo = lazy_import(
        "md_analysis.utils.RestartParser.ColvarParser", "ColvarMDInfo",
    )
    try:
        info = ColvarMDInfo.from_paths(restart_path, log_path)
    except Exception as exc:
        print(f"  (Could not parse metadata: {exc})")
        return

    cv = info.restart.colvars.primary
    print(f"\n  Trajectory info:")
    print(f"    Steps:        {info.n_steps}")
    print(f"    Timestep:     {info.restart.timestep_fs} fs")
    print(f"    CV target:    {cv.target_au:.6f} a.u.")
    print(f"    CV growth:    {cv.target_growth_au:.6e} a.u./step")
    target_start = cv.target_au + (0 - info.restart.step_start) * cv.target_growth_au
    target_end = cv.target_au + (info.n_steps - 1 - info.restart.step_start) * cv.target_growth_au
    print(f"    CV range:     [{target_start:.6f}, {target_end:.6f}] a.u.")
    print(f"    Valid index:   0 .. {info.n_steps - 1}")
    print()


# ---------------------------------------------------------------------------
# Command classes
# ---------------------------------------------------------------------------

class _SlowgrowthPlotCmd(MenuCommand):
    """Base for slow-growth plot commands."""

    output_subdir = ""
    _plot_style: str = "both"

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}

        workdir = Path(".").resolve()

        # Discover or prompt for restart file
        default_restart = _discover_restart_file(workdir)
        if default_restart:
            print(f"  Found restart: {Path(default_restart).name}")
        ctx[K.RESTART_PATH] = prompt_str(
            "COLVAR restart file", default=default_restart,
        )
        if not ctx[K.RESTART_PATH]:
            raise FileNotFoundError("No restart file specified and none discovered.")

        # Discover or prompt for log file
        default_log = _discover_log_file(workdir)
        if default_log:
            print(f"  Found log:     {Path(default_log).name}")
        ctx[K.LOG_PATH] = prompt_str(
            "LagrangeMultLog file", default=default_log,
        )
        if not ctx[K.LOG_PATH]:
            raise FileNotFoundError("No LagrangeMultLog file specified and none discovered.")

        # Show trajectory info
        _print_sg_info(ctx[K.RESTART_PATH], ctx[K.LOG_PATH])

        # Step range
        ctx[K.INITIAL_STEP] = prompt_int("Initial step (0-based)", default=0) or 0
        ctx[K.FINAL_STEP] = prompt_int(
            "Final step (exclusive, empty=all)", default=None,
        )
        ctx[K.COLVAR_ID] = prompt_int("Colvar ID (empty=primary)", default=None)
        ctx[K.OUTDIR] = prompt_str("Output directory", default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        slowgrowth_analysis = lazy_import(
            "md_analysis.enhanced_sampling.slowgrowth.SlowGrowthPlot",
            "slowgrowth_analysis",
        )
        outdir = Path(ctx[K.OUTDIR]).resolve()
        slowgrowth_analysis(
            ctx[K.RESTART_PATH],
            ctx[K.LOG_PATH],
            initial_step=ctx[K.INITIAL_STEP],
            final_step=ctx[K.FINAL_STEP],
            output_dir=outdir,
            plot_style=self._plot_style,
            colvar_id=ctx[K.COLVAR_ID],
        )


class SGQuickPlotCmd(_SlowgrowthPlotCmd):
    _plot_style = "quick"


class SGPublicationPlotCmd(_SlowgrowthPlotCmd):
    _plot_style = "publication"
