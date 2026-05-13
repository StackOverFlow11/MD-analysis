"""Enhanced Sampling command classes (301-302)."""

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
    import numpy as np

    from ..utils.constants import AU_TIME_TO_FS

    ColvarMDInfo = lazy_import(
        "md_analysis.utils.formats.cp2k.colvar", "ColvarMDInfo",
    )
    try:
        info = ColvarMDInfo.from_paths(restart_path, log_path)
    except Exception as exc:
        print(f"  (Could not parse metadata: {exc})")
        return

    cv = info.restart.colvars.primary
    dt_au = info.restart.timestep_fs / AU_TIME_TO_FS
    growth_per_step = cv.target_growth_au * dt_au
    print(f"\n  Trajectory info:")
    print(f"    Steps:        {info.n_steps}")
    print(f"    Timestep:     {info.restart.timestep_fs} fs")
    print(f"    CV target:    {cv.target_au:.6f} a.u.")
    print(f"    CV growth:    {growth_per_step:.6e} a.u./step")
    xi = info.target_series_au()
    print(f"    CV range:     [{xi[0]:.6f}, {xi[-1]:.6f}] a.u.")
    print(f"    Valid index:   0 .. {info.n_steps - 1}")

    # Warn about overflow (nan) steps
    nan_mask = np.isnan(info.lagrange.collective_shake)
    n_nan = int(np.sum(nan_mask))
    if n_nan > 0:
        nan_indices = np.where(nan_mask)[0]
        preview = ", ".join(str(i) for i in nan_indices[:10])
        if n_nan > 10:
            preview += f", ... ({n_nan} total)"
        print(f"    WARNING: {n_nan} overflow (***) steps: [{preview}]")
        print(f"    Use initial/final step to avoid these indices.")
    print()


# ---------------------------------------------------------------------------
# Command classes
# ---------------------------------------------------------------------------

class _SlowgrowthPlotCmd(MenuCommand):
    """Base for slow-growth plot commands."""

    # Workflow name to look up on md_analysis.workflows.enhanced_sampling.
    # Subclasses override; the base value matches the historical
    # ``plot_style="both"`` semantics by combining the two plot styles
    # but is never instantiated directly.
    _workflow_name: str = ""

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
        if not self._workflow_name:
            raise RuntimeError(
                f"{type(self).__name__} must set _workflow_name"
            )
        run_sg = lazy_import(
            "md_analysis.workflows.enhanced_sampling",
            self._workflow_name,
        )
        run_sg(
            restart_path=ctx[K.RESTART_PATH],
            log_path=ctx[K.LOG_PATH],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            initial_step=ctx[K.INITIAL_STEP],
            final_step=ctx[K.FINAL_STEP],
            colvar_id=ctx[K.COLVAR_ID],
        )


class SGQuickPlotCmd(_SlowgrowthPlotCmd):
    _workflow_name = "run_slowgrowth_quick_plot"


class SGPublicationPlotCmd(_SlowgrowthPlotCmd):
    _workflow_name = "run_slowgrowth_publication_plot"
