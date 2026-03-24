"""Constrained TI analysis command classes (311-312)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._framework import MenuCommand, lazy_import
from ._params import K
from ._prompt import prompt_float, prompt_int, prompt_str
from ._enhanced_sampling import _discover_restart_file, _discover_log_file


_VALID_PATTERNS = ("ti_target", "xi", "auto")


# ---------------------------------------------------------------------------
# 311 — Single-Point Diagnostics
# ---------------------------------------------------------------------------

class TISingleDiagCmd(MenuCommand):
    """Single constraint-point convergence diagnostics."""

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

        ctx[K.EQUILIBRATION] = prompt_int("Equilibration frames to discard", default=0) or 0

        # SEM target: nullable float via prompt_str
        sem_raw = prompt_str("SEM target (a.u., empty=none)", default="")
        ctx[K.SEM_TARGET] = float(sem_raw) if sem_raw else None

        ctx[K.COLVAR_ID] = prompt_int("Colvar ID (empty=primary)", default=None)
        ctx[K.OUTDIR] = prompt_str("Output directory", default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        standalone_diagnostics = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "standalone_diagnostics",
        )
        outdir = ctx[K.OUTDIR_RESOLVED]
        result = standalone_diagnostics(
            restart_path=ctx[K.RESTART_PATH],
            log_path=ctx[K.LOG_PATH],
            equilibration=ctx[K.EQUILIBRATION],
            sem_target=ctx[K.SEM_TARGET],
            colvar_id=ctx[K.COLVAR_ID],
            output_dir=outdir,
        )

        # Console summary
        r = result["report"]
        status = "N/A" if r.passed is None else ("PASS" if r.passed else "FAIL")
        method = "plateau" if r.block_avg.plateau_reached else "acf"
        print(f"\n  ξ = {r.xi:.6f} a.u.")
        print(f"  ⟨λ⟩ = {r.lambda_mean:.6f} a.u.")
        print(f"  τ_corr = {r.autocorr.tau_corr:.1f} frames, N_eff = {r.autocorr.n_eff:.1f}")
        print(f"  SEM_final = {r.sem_final:.6f} a.u. (method: {method})")
        print(f"  Geweke |z| = {abs(r.geweke.z):.3f}  ({'PASS' if r.geweke.passed else 'FAIL'})")
        print(f"  Overall: {status}")

        if r.failure_reasons:
            for reason in r.failure_reasons:
                print(f"    - {reason}")

        for key in ("diagnostics_png", "csv"):
            if key in result:
                print(f"  {key}: {result[key]}")


# ---------------------------------------------------------------------------
# 312 — Full TI Analysis
# ---------------------------------------------------------------------------

class TIFullAnalysisCmd(MenuCommand):
    """Multi-point constrained TI convergence analysis + free-energy integration."""

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}

        ctx[K.TI_ROOT_DIR] = prompt_str("TI root directory", default=".") or "."

        # Pattern with validation
        while True:
            pattern = prompt_str(
                f"Directory pattern ({'/'.join(_VALID_PATTERNS)})",
                default="auto",
            ) or "auto"
            if pattern in _VALID_PATTERNS:
                break
            print(f"  Invalid pattern '{pattern}'. Must be one of: {', '.join(_VALID_PATTERNS)}")
        ctx[K.TI_DIR_PATTERN] = pattern

        ctx[K.EQUILIBRATION] = prompt_int(
            "Equilibration frames to discard (all points)", default=0,
        ) or 0
        ctx[K.EPSILON_TOL_EV] = prompt_float(
            "Free-energy tolerance ε (eV)", default=0.05,
        )
        ctx[K.OUTDIR] = prompt_str("Output directory", default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        discover_ti_points = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.io",
            "discover_ti_points",
        )
        load_ti_series = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.io",
            "load_ti_series",
        )
        analyze_ti = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "analyze_ti",
        )
        write_convergence_csv = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "write_convergence_csv",
        )
        write_free_energy_csv = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.workflow",
            "write_free_energy_csv",
        )
        plot_point_diagnostics = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.plot",
            "plot_point_diagnostics",
        )
        plot_free_energy_profile = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.plot",
            "plot_free_energy_profile",
        )

        outdir = ctx[K.OUTDIR_RESOLVED]
        root_dir = Path(ctx[K.TI_ROOT_DIR])

        # 1. Discover constraint points
        point_defs = discover_ti_points(root_dir, pattern=ctx[K.TI_DIR_PATTERN])
        print(f"\n  Found {len(point_defs)} constraint points:")
        for p in point_defs:
            print(f"    ξ = {p.xi:.6f}")

        # 2. Load series
        series_data = load_ti_series(point_defs)
        xi_values = np.array([x for x, _, _ in series_data])
        lambda_list = [s for _, s, _ in series_data]

        # 3. dt consistency check
        dts = [d for _, _, d in series_data]
        unique_dts = set(f"{d:.10f}" for d in dts)
        if len(unique_dts) > 1:
            print(f"\n  WARNING: inconsistent timesteps detected: {set(dts)}")
            print(f"  Using dt = {dts[0]:.6f} fs from first point")
        dt = dts[0]
        print(f"  Timestep: {dt:.6f} fs")

        # 4. Analyze
        ti_report = analyze_ti(
            xi_values, lambda_list, dt,
            epsilon_tol_ev=ctx[K.EPSILON_TOL_EV],
            equilibration=ctx[K.EQUILIBRATION],
        )

        # 5. Console summary table
        from ...utils.config import HA_TO_EV
        print(f"\n  {'Point':<6} {'ξ':<12} {'⟨λ⟩':<14} {'SEM':<14} {'Status'}")
        print(f"  {'─' * 58}")
        for i, r in enumerate(ti_report.point_reports):
            status = "PASS" if r.passed else "FAIL"
            print(f"  {i:<6} {r.xi:<12.6f} {r.lambda_mean:<14.6f} {r.sem_final:<14.6f} {status}")

        delta_A_ev = ti_report.delta_A * HA_TO_EV
        sigma_A_ev = ti_report.sigma_A * HA_TO_EV
        print(f"\n  ΔA = {delta_A_ev:.6f} ± {sigma_A_ev:.6f} eV")

        if ti_report.all_passed:
            print("  Status: ALL PASS")
            print("  All points converged. Use 311 for per-point diagnostics.")
        else:
            failing = ti_report.failing_indices
            print(f"  Status: {len(failing)} FAILED (indices: {', '.join(str(i) for i in failing)})")
            if ti_report.suggested_time_ratios is not None:
                print("  Suggested time allocation (relative):")
                parts = [f"ξ={xi_values[i]:.4f}: {ti_report.suggested_time_ratios[i]:.2f}"
                         for i in range(len(xi_values))]
                print(f"    {', '.join(parts)}")

        # 6. Write output files
        write_convergence_csv(ti_report, output_dir=outdir)
        write_free_energy_csv(ti_report, output_dir=outdir)
        plot_free_energy_profile(ti_report, output_dir=outdir)

        # 7. Diagnostics plots for failing points only
        for i in ti_report.failing_indices:
            plot_point_diagnostics(ti_report.point_reports[i], output_dir=outdir)

        # 8. List output files
        print(f"\n  Output files:")
        print(f"    {outdir / 'ti_free_energy.png'}")
        print(f"    {outdir / 'ti_free_energy.csv'}")
        print(f"    {outdir / 'ti_convergence_report.csv'}")
        for i in ti_report.failing_indices:
            xi = ti_report.point_reports[i].xi
            print(f"    {outdir / f'ti_diag_xi{xi:.4f}.png'}")
