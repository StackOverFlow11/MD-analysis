"""Pure plotting module for constrained-TI diagnostics.

Never imports workflow.py — dependency is one-way (workflow -> plot).
Follows SlowGrowthPlot.py conventions: lazy matplotlib, Agg backend,
dpi=180, returns Path.
"""

from __future__ import annotations

from pathlib import Path

from .config import DEFAULT_DIAGNOSTICS_PNG_PREFIX, DEFAULT_FE_PROFILE_PNG_NAME
from .models import ConstraintPointReport, TIReport
from ...utils.config import HA_TO_EV


def plot_point_diagnostics(
    report: ConstraintPointReport,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Generate 2x2 diagnostic plot for one constraint point.

    Parameters
    ----------
    report : ConstraintPointReport
    output_dir : Path or None

    Returns
    -------
    Path to saved PNG.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    filename = f"{DEFAULT_DIAGNOSTICS_PNG_PREFIX}_xi{report.xi:.4f}.png"
    out_path = out / filename

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=180)

    # --- Top-left: Running average ---
    ax = axes[0, 0]
    rm = report.running_avg.running_mean
    n = len(rm)
    ax.plot(np.arange(1, n + 1), rm, linewidth=0.5, color="C0")
    final_mean = report.lambda_mean
    ax.axhline(final_mean, color="grey", linestyle="--", linewidth=0.8)
    ax.axhline(
        final_mean + report.sem_final,
        color="C1",
        linestyle=":",
        linewidth=0.7,
        label=f"mean ± SEM ({report.sem_final:.2e})",
    )
    ax.axhline(final_mean - report.sem_final, color="C1", linestyle=":", linewidth=0.7)
    ax.set_xlabel("Frame count")
    ax.set_ylabel("Cumulative mean λ")
    ax.set_title(f"Running Average (D={report.running_avg.drift_D:.2e})")
    ax.legend(fontsize=7)

    # --- Top-right: ACF ---
    ax = axes[0, 1]
    acf = report.autocorr.acf
    max_lag = min(len(acf), int(10 * report.autocorr.tau_corr) + 1, len(acf))
    ax.plot(np.arange(max_lag), acf[:max_lag], linewidth=0.8, color="C0")
    ax.axhline(0, color="grey", linestyle="-", linewidth=0.5)
    # Mark truncation point at alpha * tau
    m_cut = min(int(5 * report.autocorr.tau_corr), max_lag - 1)
    ax.axvline(m_cut, color="C3", linestyle="--", linewidth=0.8, label=f"M=5τ ({m_cut})")
    ax.set_xlabel("Lag j")
    ax.set_ylabel("C(j)")
    ax.set_title(f"ACF (τ_corr={report.autocorr.tau_corr:.1f})")
    ax.legend(fontsize=7)

    # --- Bottom-left: Block average SEM(B) — Flyvbjerg-Petersen ---
    ax = axes[1, 0]
    ba = report.block_avg
    bs = ba.block_sizes
    valid = ~np.isnan(ba.sem_curve)

    # Split points into plateau (green squares) and non-plateau (blue circles)
    if ba.plateau_reached and ba.plateau_index is not None:
        pre = np.zeros(len(bs), dtype=bool)
        pre[:ba.plateau_index] = True
        pre &= valid
        plat = np.zeros(len(bs), dtype=bool)
        plat[ba.plateau_index:] = True
        plat &= valid
        # Non-plateau points
        ax.errorbar(
            bs[pre], ba.sem_curve[pre], yerr=ba.delta_sem[pre],
            fmt="o", markersize=4, color="C0", capsize=3,
            ecolor="C0", elinewidth=0.8,
        )
        # Plateau points (highlighted)
        ax.errorbar(
            bs[plat], ba.sem_curve[plat], yerr=ba.delta_sem[plat],
            fmt="s", markersize=5, color="C2", capsize=3,
            ecolor="C2", elinewidth=0.8, label="plateau",
        )
        # Plateau line and band
        ax.axhline(ba.plateau_sem, color="C2", linestyle="--", linewidth=0.8,
                    label=f"SEM_block ({ba.plateau_sem:.2e})")
        ax.axhspan(
            ba.plateau_sem - ba.plateau_delta,
            ba.plateau_sem + ba.plateau_delta,
            alpha=0.12, color="C2",
        )
    else:
        ax.errorbar(
            bs[valid], ba.sem_curve[valid], yerr=ba.delta_sem[valid],
            fmt="o-", markersize=4, color="C0", capsize=3,
            ecolor="C0", elinewidth=0.8, label="SEM(B) ± δSEM",
        )

    ax.set_xscale("log", base=2)
    ax.set_xlabel("Block size B")
    ax.set_ylabel("SEM(B)")
    if report.sem_max is not None:
        ax.axhline(
            report.sem_max, color="C3", linestyle="--", linewidth=0.8,
            label=f"SEM_max ({report.sem_max:.2e})",
        )
    # Annotate n_b on each point
    for i in range(len(bs)):
        if not valid[i]:
            continue
        nb = ba.n_total // bs[i]
        ax.annotate(
            f"n={nb}", (bs[i], ba.sem_curve[i]),
            textcoords="offset points", xytext=(0, 8),
            fontsize=5, ha="center", color="0.5",
        )
    ax.set_title(
        f"Block Average (plateau={'yes' if ba.plateau_reached else 'no'})"
    )
    ax.legend(fontsize=7)

    # --- Bottom-right: Summary table (4 sections) ---
    ax = axes[1, 1]
    ax.axis("off")
    r = report
    ba = r.block_avg
    n = ba.n_total
    sigma = r.sigma_lambda

    # Block-average implied tau and N_eff:
    # SEM_block = sigma * sqrt(2*tau_implied/N)  =>  tau_implied = N*(SEM/sigma)^2/2
    if sigma > 0 and ba.plateau_sem > 0:
        tau_block = n * (ba.plateau_sem / sigma) ** 2 / 2.0
        neff_block = n / (2.0 * tau_block) if tau_block > 0 else float("inf")
    else:
        tau_block = float("nan")
        neff_block = float("nan")

    # λ = dA/dξ has units of Hartree/(ξ unit), NOT pure energy.
    # Keep λ-related quantities in a.u. (consistent with plot axes).

    # Section 1: Overview
    sec1 = [
        "── Overview ──",
        f"ξ = {r.xi:.6f}    N = {r.n_analyzed}",
        f"t = {r.time_start_fs:.1f} – {r.time_end_fs:.1f} fs",
        f"⟨λ⟩ = {r.lambda_mean:.6f} a.u.",
        f"SEM  = {r.sem_final:.6f} a.u.",
    ]
    if r.passed is not None:
        sec1.append(f"OVERALL: {'PASS' if r.passed else 'FAIL'}")
    else:
        sec1.append("OVERALL: N/A (no target)")

    # Section 2: ACF
    sec2 = [
        "── ACF ──",
        f"τ_corr = {r.autocorr.tau_corr:.1f} frames",
        f"N_eff  = {r.autocorr.n_eff:.1f}  ({'PASS' if r.autocorr.passed_neff else 'FAIL'})",
        f"SEM    = {r.autocorr.sem_auto:.6f} a.u.",
    ]

    # Section 3: Block Average (F&P)
    plat_info = (
        f"yes, B={ba.plateau_block_size}" if ba.plateau_reached else "no"
    )
    sec3 = [
        "── Block Average ──",
        f"plateau: {plat_info}",
        f"SEM    = {ba.plateau_sem:.6f} ± {ba.plateau_delta:.6f} a.u.",
        f"τ_impl = {tau_block:.1f} frames",
        f"N_impl = {neff_block:.1f}",
    ]

    # Section 4: Geweke
    sec4 = [
        "── Geweke ──",
        f"|z| = {abs(r.geweke.z):.3f}  ({'PASS' if r.geweke.passed else 'FAIL'})"
        + ("" if r.geweke.reliable else " [unreliable]"),
        f"mean_A = {r.geweke.mean_a:.6f} a.u.",
        f"mean_B = {r.geweke.mean_b:.6f} a.u.",
    ]

    lines = sec1 + [""] + sec2 + [""] + sec3 + [""] + sec4

    ax.text(
        0.05,
        0.98,
        "\n".join(lines),
        transform=ax.transAxes,
        fontsize=7,
        verticalalignment="top",
        fontfamily="monospace",
    )

    fig.suptitle(f"Constrained TI Diagnostics — ξ = {report.xi:.4f}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_free_energy_profile(
    ti_report: TIReport,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Plot free-energy profile with error bars.

    Parameters
    ----------
    ti_report : TIReport
    output_dir : Path or None

    Returns
    -------
    Path to saved PNG.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    out_path = out / DEFAULT_FE_PROFILE_PNG_NAME

    fig, ax1 = plt.subplots(figsize=(8, 5), dpi=180)

    xi = ti_report.xi_values
    forces = ti_report.forces  # a.u. (Hartree / ξ_unit)
    errors = ti_report.force_errors

    # Left axis: dA/dxi with error bars (a.u.)
    colors = [
        "C2" if r.passed else "C3"
        for r in ti_report.point_reports
    ]
    ax1.errorbar(xi, forces, yerr=errors, fmt="o", markersize=4, color="C0", ecolor="C0", capsize=3)
    for i, (x, y) in enumerate(zip(xi, forces)):
        ax1.plot(x, y, "o", markersize=5, color=colors[i], zorder=5)
    ax1.set_xlabel("ξ (a.u.)")
    if len(xi) >= 2 and xi[0] > xi[-1]:
        ax1.invert_xaxis()
    ax1.set_ylabel("dA/dξ (a.u.)", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")

    # Right axis: integrated A(xi) in eV
    ax2 = ax1.twinx()
    cumul_A = np.cumsum(ti_report.weights * forces) * HA_TO_EV
    cumul_sigma = np.sqrt(np.cumsum(ti_report.weights**2 * errors**2)) * HA_TO_EV
    ax2.plot(xi, cumul_A, "-s", color="C1", markersize=3, linewidth=1.2, label="A(ξ)")
    ax2.fill_between(
        xi,
        cumul_A - cumul_sigma,
        cumul_A + cumul_sigma,
        alpha=0.2,
        color="C1",
    )
    ax2.set_ylabel("A(ξ) (eV)", color="C1")
    ax2.tick_params(axis="y", labelcolor="C1")

    delta_A_eV = ti_report.delta_A * HA_TO_EV
    sigma_A_eV = ti_report.sigma_A * HA_TO_EV
    title = (
        f"ΔA = {delta_A_eV:.6f} ± {sigma_A_eV:.6f} eV"
        f"  ({'ALL PASS' if ti_report.all_passed else f'{len(ti_report.failing_indices)} FAILED'})"
    )
    ax1.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
