"""Slow-growth free-energy plotting and CSV export.

Public API
----------
- ``slowgrowth_analysis``          — parse, slice, plot + CSV in one call
- ``plot_slowgrowth_quick``        — quick dual-axis plot with MD-step top axis
- ``plot_slowgrowth_publication``   — publication-quality dual-axis plot
- ``write_slowgrowth_csv``         — tabular CSV export
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ...utils._io_helpers import _write_csv
from ...utils.config import HA_TO_EV
from .SlowGrowth import Slowgrowth, SlowgrowthFull
from .config import (
    DEFAULT_SG_CSV_NAME,
    DEFAULT_SG_PUBLICATION_PNG_NAME,
    DEFAULT_SG_QUICK_PNG_NAME,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Secondary-axis factories (CV <-> step linear mapping)
# ---------------------------------------------------------------------------

def _cv_to_step_factory(
    cv: np.ndarray, steps: np.ndarray,
):
    """Return a callable that maps CV values to step numbers (linear)."""
    cv0, cv1 = float(cv[0]), float(cv[-1])
    s0, s1 = float(steps[0]), float(steps[-1])
    if abs(cv1 - cv0) < 1e-30:
        return lambda x: np.full_like(np.asarray(x, dtype=float), s0)
    slope = (s1 - s0) / (cv1 - cv0)

    def forward(x):
        return s0 + (np.asarray(x, dtype=float) - cv0) * slope
    return forward


def _step_to_cv_factory(
    cv: np.ndarray, steps: np.ndarray,
):
    """Return a callable that maps step numbers to CV values (linear)."""
    cv0, cv1 = float(cv[0]), float(cv[-1])
    s0, s1 = float(steps[0]), float(steps[-1])
    if abs(s1 - s0) < 1e-30:
        return lambda x: np.full_like(np.asarray(x, dtype=float), cv0)
    slope = (cv1 - cv0) / (s1 - s0)

    def inverse(x):
        return cv0 + (np.asarray(x, dtype=float) - s0) * slope
    return inverse


# ---------------------------------------------------------------------------
# Quick plot
# ---------------------------------------------------------------------------

def plot_slowgrowth_quick(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    png_name: str = DEFAULT_SG_QUICK_PNG_NAME,
    absolute_steps: np.ndarray | None = None,
    ma_window: int = 50,
) -> Path:
    """Dual-axis quick plot with CV bottom axis and MD-step top axis.

    Left y-axis: free energy (eV), right y-axis: Lagrange multiplier (a.u.).

    Parameters
    ----------
    absolute_steps : np.ndarray, optional
        Original (pre-reversal) step indices for annotating the barrier
        peak with its absolute MD step number.  Falls back to ``sg.steps``.
    ma_window : int, optional
        Window size for the Lagrange multiplier moving average (default 50).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    fe_ev = sg.free_energy_au * HA_TO_EV
    abs_steps = absolute_steps if absolute_steps is not None else sg.steps

    fig, ax_left = plt.subplots(figsize=(8, 5), dpi=180)

    # Left axis — free energy
    ax_left.plot(
        sg.target_au, fe_ev,
        color="tab:red", lw=2.0,
    )
    ax_left.set_xlabel("Collective Variable (a.u.)")
    ax_left.set_ylabel("Free Energy (eV)", color="tab:red")
    ax_left.tick_params(axis="y", labelcolor="tab:red")

    # Right axis — Lagrange multiplier
    ax_right = ax_left.twinx()
    ax_right.plot(
        sg.target_au, sg.lagrange_shake,
        color="tab:blue", lw=0.8, alpha=0.7,
    )
    # Moving average of Lagrange multiplier
    if ma_window > 1 and len(sg.lagrange_shake) >= ma_window:
        kernel = np.ones(ma_window) / ma_window
        lm_ma = np.convolve(sg.lagrange_shake, kernel, mode="valid")
        ma_x = sg.target_au[ma_window - 1:]
        ax_right.plot(ma_x, lm_ma, color="darkblue", lw=1.5, label=f"MA({ma_window})")
        # Annotate last MA value
        last_ma = float(lm_ma[-1])
        ax_right.annotate(
            f"MA = {last_ma:.4f}",
            xy=(ma_x[-1], last_ma),
            xytext=(ma_x[-1], last_ma + 0.1 * (ax_right.get_ylim()[1] - ax_right.get_ylim()[0])),
            arrowprops=dict(arrowstyle="->", color="darkblue", lw=1.5),
            fontsize=9, color="darkblue", ha="center",
        )
    ax_right.set_ylabel("Lagrange Multiplier (a.u.)", color="tab:blue")
    ax_right.tick_params(axis="y", labelcolor="tab:blue")

    # Top axis — MD step (secondary x)
    ax_top = ax_left.secondary_xaxis(
        "top",
        functions=(
            _cv_to_step_factory(sg.target_au, sg.steps),
            _step_to_cv_factory(sg.target_au, sg.steps),
        ),
    )
    ax_top.set_xlabel("MD Step")

    # Ensure time flows left-to-right: invert x-axis when CV is decreasing
    if len(sg.target_au) > 1 and sg.target_au[-1] < sg.target_au[0]:
        ax_left.invert_xaxis()

    # Barrier annotation — arrow pointing to peak with absolute step
    peak_idx = int(np.nanargmax(fe_ev))
    delta_f_barrier = float(fe_ev[peak_idx] - fe_ev[0])
    peak_abs_step = int(abs_steps[peak_idx])
    ax_left.annotate(
        f"$\\Delta F^\\dagger$ = {delta_f_barrier:.4f} eV\n(step {peak_abs_step})",
        xy=(sg.target_au[peak_idx], fe_ev[peak_idx]),
        xytext=(sg.target_au[peak_idx], fe_ev[peak_idx] + 0.15 * abs(delta_f_barrier or 0.01)),
        arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.5),
        fontsize=9, color="tab:red", ha="center",
    )

    # Total free-energy annotation — arrow pointing to endpoint
    delta_f_total = float(fe_ev[-1] - fe_ev[0])
    ax_left.annotate(
        f"$\\Delta F$ = {delta_f_total:.4f} eV",
        xy=(sg.target_au[-1], fe_ev[-1]),
        xytext=(sg.target_au[-1], fe_ev[-1] + 0.15 * abs(delta_f_total or 0.01)),
        arrowprops=dict(arrowstyle="->", color="darkred", lw=1.5),
        fontsize=9, color="darkred", ha="center",
    )

    ax_left.grid(True, alpha=0.25)
    fig.tight_layout()
    out_path = outdir / png_name
    fig.savefig(out_path)
    plt.close(fig)

    logger.info("Quick plot saved: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Publication plot
# ---------------------------------------------------------------------------

def plot_slowgrowth_publication(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    png_name: str = DEFAULT_SG_PUBLICATION_PNG_NAME,
) -> Path:
    """Publication-quality dual-axis plot with energy values in legend.

    Left y-axis: free energy (eV), right y-axis: Lagrange multiplier (a.u.).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import rc_context

    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    fe_ev = sg.free_energy_au * HA_TO_EV

    rc_overrides = {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
    }

    with rc_context(rc_overrides):
        fig, ax_left = plt.subplots(figsize=(8, 5), dpi=180)

        # Left axis — free energy
        ax_left.plot(
            sg.target_au, fe_ev,
            color="tab:red", lw=2.0,
        )
        ax_left.set_xlabel("Collective Variable (a.u.)")
        ax_left.set_ylabel("Free Energy (eV)", color="tab:red")
        ax_left.tick_params(axis="y", labelcolor="tab:red")

        # Right axis — Lagrange multiplier
        ax_right = ax_left.twinx()
        ax_right.plot(
            sg.target_au, sg.lagrange_shake,
            color="tab:blue", lw=0.8, alpha=0.7,
        )
        ax_right.set_ylabel("Lagrange Multiplier (a.u.)", color="tab:blue")
        ax_right.tick_params(axis="y", labelcolor="tab:blue")

        # Ensure time flows left-to-right: invert x-axis when CV is decreasing
        if len(sg.target_au) > 1 and sg.target_au[-1] < sg.target_au[0]:
            ax_left.invert_xaxis()

        # Legend with energy values (invisible handles)
        peak_idx = int(np.nanargmax(fe_ev))
        delta_f_barrier = float(fe_ev[peak_idx] - fe_ev[0])
        delta_f_total = float(fe_ev[-1] - fe_ev[0])
        legend_handles = [
            Line2D([], [], color="none",
                   label=f"$\\Delta F^\\dagger$ = {delta_f_barrier:.4f} eV"),
            Line2D([], [], color="none",
                   label=f"$\\Delta F$ = {delta_f_total:.4f} eV"),
        ]
        ax_left.legend(handles=legend_handles, loc="best",
                       frameon=True, handlelength=0)

        ax_left.grid(True, alpha=0.25)
        fig.tight_layout()
        out_path = outdir / png_name
        fig.savefig(out_path)
        plt.close(fig)

    logger.info("Publication plot saved: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def write_slowgrowth_csv(
    sg: Slowgrowth,
    *,
    output_dir: Path | None = None,
    csv_name: str = DEFAULT_SG_CSV_NAME,
) -> Path:
    """Write slow-growth data to a CSV file."""
    outdir = (output_dir or Path(".")).resolve()
    out_path = outdir / csv_name

    fe_ev = sg.free_energy_au * HA_TO_EV
    fieldnames = [
        "step", "time_fs", "target_au",
        "lagrange_au", "free_energy_au", "free_energy_ev",
    ]
    rows = [
        {
            "step": int(sg.steps[i]),
            "time_fs": float(sg.times_fs[i]),
            "target_au": float(sg.target_au[i]),
            "lagrange_au": float(sg.lagrange_shake[i]),
            "free_energy_au": float(sg.free_energy_au[i]),
            "free_energy_ev": float(fe_ev[i]),
        }
        for i in range(sg.n_steps)
    ]
    _write_csv(out_path, rows, fieldnames)
    logger.info("CSV saved: %s (%d rows)", out_path, len(rows))
    return out_path


# ---------------------------------------------------------------------------
# Unified entry point
# ---------------------------------------------------------------------------

def slowgrowth_analysis(
    restart_path: str,
    log_path: str,
    *,
    initial_step: int = 0,
    final_step: int | None = None,
    output_dir: Path | None = None,
    plot_style: str = "both",
    colvar_id: int | None = None,
) -> dict[str, Path]:
    """Parse, slice, optionally reverse, then plot + export CSV.

    Parameters
    ----------
    restart_path, log_path : str
        Paths to the CP2K COLVAR restart and LagrangeMultLog files.
    initial_step, final_step : int
        Array indices (0-based) for the segment ``[initial, final)``.
        If *initial_step* > *final_step*, the segment is reversed.
    output_dir : Path, optional
        Output directory (default: current directory).
    plot_style : str
        ``"quick"``, ``"publication"``, or ``"both"``.
    colvar_id : int, optional
        Which collective variable to use (default: primary).

    Returns
    -------
    dict[str, Path]
        Keys may include ``"quick_png"``, ``"publication_png"``, ``"csv"``.
    """
    full = SlowgrowthFull.from_paths(restart_path, log_path, colvar_id=colvar_id)
    logger.info(
        "Parsed slow-growth: %d steps, timestep=%.2f fs, "
        "CV range [%.6f, %.6f] a.u.",
        full.n_steps, full.timestep_fs,
        float(full.target_au[0]), float(full.target_au[-1]),
    )

    if final_step is None:
        final_step = full.n_steps

    if initial_step > final_step:
        pre_rev = full.segment(final_step, initial_step)
        absolute_steps = pre_rev.steps[::-1].copy()  # preserve original indices
        seg = pre_rev.reversed()
        logger.info(
            "Reversed segment [%d, %d) → %d steps", final_step, initial_step, seg.n_steps,
        )
    else:
        seg = full.segment(initial_step, final_step)
        absolute_steps = seg.steps.copy()
        logger.info("Segment [%d, %d) → %d steps", initial_step, final_step, seg.n_steps)

    outdir = (output_dir or Path(".")).resolve()
    results: dict[str, Path] = {}

    results["csv"] = write_slowgrowth_csv(seg, output_dir=outdir)

    if plot_style in ("quick", "both"):
        results["quick_png"] = plot_slowgrowth_quick(
            seg, output_dir=outdir, absolute_steps=absolute_steps,
        )
    if plot_style in ("publication", "both"):
        results["publication_png"] = plot_slowgrowth_publication(seg, output_dir=outdir)

    fe_ev = seg.free_energy_au * HA_TO_EV
    print(f"\n  Steps: {seg.n_steps}")
    print(f"  CV range: [{float(seg.target_au[0]):.6f}, {float(seg.target_au[-1]):.6f}] a.u.")
    print(f"  ΔF† = {float(np.max(fe_ev) - fe_ev[0]):.4f} eV")
    print(f"  ΔF  = {float(fe_ev[-1] - fe_ev[0]):.4f} eV")

    for key, path in results.items():
        print(f"  {key}: {path}")

    return results
