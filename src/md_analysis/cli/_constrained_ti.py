"""Constrained TI analysis command classes (311-313)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._framework import MenuCommand, lazy_import
from ._params import K
from ._prompt import prompt_bool, prompt_choice, prompt_float, prompt_int, prompt_str
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
        n_total = result["n_total"]
        n_analyzed = result["n_analyzed"]
        equil = result["equilibration"]
        t_start = result["time_start_fs"]
        t_end = result["time_end_fs"]
        print(f"  Frames: {n_total} total, {equil} discarded, {n_analyzed} analyzed")
        print(f"  Time range: {t_start:.1f} – {t_end:.1f} fs")
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
            "Default equilibration frames to discard", default=0,
        ) or 0
        ctx[K.EPSILON_TOL_EV] = prompt_float(
            "Free-energy tolerance ε (eV)", default=0.05,
        )
        ctx[K.TI_REVERSE] = prompt_bool(
            "Reverse integration direction (initial state = max ξ)?",
            default=False,
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
        point_defs = discover_ti_points(
            root_dir,
            pattern=ctx[K.TI_DIR_PATTERN],
            reverse=ctx[K.TI_REVERSE],
        )
        print(f"\n  Found {len(point_defs)} constraint points:")
        for p in point_defs:
            print(f"    ξ = {p.xi:.6f}")

        # Per-point equilibration
        default_equil = ctx[K.EQUILIBRATION]
        if len(point_defs) > 1 and prompt_bool(
            "Set per-point equilibration frames?", default=False,
        ):
            equil_list = []
            for p in point_defs:
                val = prompt_int(
                    f"  ξ={p.xi:.6f} equilibration frames",
                    default=default_equil,
                )
                equil_list.append(val if val is not None else default_equil)
            equilibration = equil_list
        else:
            equilibration = default_equil

        # 2. Load series + parse time_start for each point
        parse_colvar_restart = lazy_import(
            "md_analysis.utils.RestartParser.ColvarParser",
            "parse_colvar_restart",
        )
        series_data = load_ti_series(point_defs)
        xi_values = np.array([x for x, _, _ in series_data])
        lambda_list = [s for _, s, _ in series_data]
        time_starts = [
            parse_colvar_restart(str(p.restart_path)).time_start_fs
            for p in point_defs
        ]

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
            equilibration=equilibration,
            time_starts=time_starts,
        )

        # 5. Console summary table
        from ..utils.config import HA_TO_EV
        print(f"\n  {'Point':<6} {'ξ':<12} {'⟨λ⟩':<14} {'SEM':<14} {'N':<8} {'Time range (fs)':<24} {'Status'}")
        print(f"  {'─' * 82}")
        for i, r in enumerate(ti_report.point_reports):
            status = "PASS" if r.passed else "FAIL"
            print(
                f"  {i:<6} {r.xi:<12.6f} {r.lambda_mean:<14.6f} {r.sem_final:<14.6f} "
                f"{r.n_analyzed:<8} {r.time_start_fs:.1f} – {r.time_end_fs:.1f}"
                f"{'':>4}{status}"
            )

        delta_A_ev = ti_report.delta_A * HA_TO_EV
        sigma_A_ev = ti_report.sigma_A * HA_TO_EV
        print(f"\n  ΔA = {delta_A_ev:.6f} ± {sigma_A_ev:.6f} eV")

        if ti_report.all_passed:
            print("  Status: ALL PASS")
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

        # 7. Diagnostics plots for all points
        for r in ti_report.point_reports:
            plot_point_diagnostics(r, output_dir=outdir)

        # 8. List output files
        print(f"\n  Output files:")
        print(f"    {outdir / 'ti_free_energy.png'}")
        print(f"    {outdir / 'ti_free_energy.csv'}")
        print(f"    {outdir / 'ti_convergence_report.csv'}")
        for r in ti_report.point_reports:
            print(f"    {outdir / f'ti_diag_xi{r.xi:.4f}.png'}")


# ---------------------------------------------------------------------------
# 313 — Constant-Potential Free Energy Correction
# ---------------------------------------------------------------------------

_VALID_SIDES = ("aligned", "opposed")
_VALID_METHODS = ("counterion", "layer")


class TIConstPotCorrectionCmd(MenuCommand):
    """Constant-potential free energy correction (Norskov) on top of TI analysis."""

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}

        # --- Same as 312 ---
        ctx[K.TI_ROOT_DIR] = prompt_str("TI root directory", default=".") or "."

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
            "Default equilibration frames to discard", default=0,
        ) or 0
        ctx[K.EPSILON_TOL_EV] = prompt_float(
            "Free-energy tolerance ε (eV)", default=0.05,
        )
        ctx[K.TI_REVERSE] = prompt_bool(
            "Reverse integration direction (initial state = max ξ)?",
            default=False,
        )

        # --- Correction-specific ---
        ctx[K.TARGET_SIDE] = prompt_choice(
            "Target electrode surface", list(_VALID_SIDES), default="aligned",
        )
        cal_raw = prompt_str("Calibration JSON path (empty=default)", default="")
        ctx[K.CALIBRATION_JSON] = cal_raw if cal_raw else None

        # Advanced
        if prompt_bool("Modify advanced charge parameters?", default=False):
            ctx[K.NORMAL] = prompt_choice(
                "Surface normal axis", ["a", "b", "c"], default="c",
            )
            ctx[K.METHOD] = prompt_choice(
                "Charge method", list(_VALID_METHODS), default="counterion",
            )
        else:
            ctx[K.NORMAL] = "c"
            ctx[K.METHOD] = "counterion"

        ctx[K.OUTDIR] = prompt_str("Output directory", default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        # --- Lazy imports ---
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
        compute_constant_potential_correction = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.correction",
            "compute_constant_potential_correction",
        )
        write_corrected_free_energy_csv = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.correction",
            "write_corrected_free_energy_csv",
        )
        plot_corrected_free_energy_profile = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.correction",
            "plot_corrected_free_energy_profile",
        )

        outdir = ctx[K.OUTDIR_RESOLVED]
        root_dir = Path(ctx[K.TI_ROOT_DIR])

        # ===== Phase 1: Full 312 TI analysis =====

        # 1. Discover constraint points
        point_defs = discover_ti_points(
            root_dir,
            pattern=ctx[K.TI_DIR_PATTERN],
            reverse=ctx[K.TI_REVERSE],
        )
        print(f"\n  Found {len(point_defs)} constraint points:")
        for p in point_defs:
            print(f"    ξ = {p.xi:.6f}")

        # Per-point equilibration
        default_equil = ctx[K.EQUILIBRATION]
        if len(point_defs) > 1 and prompt_bool(
            "Set per-point equilibration frames?", default=False,
        ):
            equil_list = []
            for p in point_defs:
                val = prompt_int(
                    f"  ξ={p.xi:.6f} equilibration frames",
                    default=default_equil,
                )
                equil_list.append(val if val is not None else default_equil)
            equilibration = equil_list
        else:
            equilibration = default_equil

        # 2. Load series
        parse_colvar_restart = lazy_import(
            "md_analysis.utils.RestartParser.ColvarParser",
            "parse_colvar_restart",
        )
        series_data = load_ti_series(point_defs)
        xi_values = np.array([x for x, _, _ in series_data])
        lambda_list = [s for _, s, _ in series_data]
        time_starts = [
            parse_colvar_restart(str(p.restart_path)).time_start_fs
            for p in point_defs
        ]

        # 3. dt consistency check
        dts = [d for _, _, d in series_data]
        unique_dts = set(f"{d:.10f}" for d in dts)
        if len(unique_dts) > 1:
            print(f"\n  WARNING: inconsistent timesteps detected: {set(dts)}")
            print(f"  Using dt = {dts[0]:.6f} fs from first point")
        dt = dts[0]
        print(f"  Timestep: {dt:.6f} fs")

        # 4. Analyze TI
        ti_report = analyze_ti(
            xi_values, lambda_list, dt,
            epsilon_tol_ev=ctx[K.EPSILON_TOL_EV],
            equilibration=equilibration,
            time_starts=time_starts,
        )

        # 5. Console summary (same as 312)
        from ..utils.config import HA_TO_EV
        print(f"\n  {'Point':<6} {'ξ':<12} {'⟨λ⟩':<14} {'SEM':<14} {'N':<8} {'Time range (fs)':<24} {'Status'}")
        print(f"  {'─' * 82}")
        for i, r in enumerate(ti_report.point_reports):
            status = "PASS" if r.passed else "FAIL"
            print(
                f"  {i:<6} {r.xi:<12.6f} {r.lambda_mean:<14.6f} {r.sem_final:<14.6f} "
                f"{r.n_analyzed:<8} {r.time_start_fs:.1f} – {r.time_end_fs:.1f}"
                f"{'':>4}{status}"
            )

        delta_A_ev = ti_report.delta_A * HA_TO_EV
        sigma_A_ev = ti_report.sigma_A * HA_TO_EV
        print(f"\n  ΔA (const-q) = {delta_A_ev:.6f} ± {sigma_A_ev:.6f} eV")

        if ti_report.all_passed:
            print("  Status: ALL PASS")
        else:
            failing = ti_report.failing_indices
            print(f"  Status: {len(failing)} FAILED (indices: {', '.join(str(i) for i in failing)})")
            if ti_report.suggested_time_ratios is not None:
                print("  Suggested time allocation (relative):")
                parts = [f"ξ={xi_values[i]:.4f}: {ti_report.suggested_time_ratios[i]:.2f}"
                         for i in range(len(xi_values))]
                print(f"    {', '.join(parts)}")

        # 6. Write 312 output files
        write_convergence_csv(ti_report, output_dir=outdir)
        write_free_energy_csv(ti_report, output_dir=outdir)
        plot_free_energy_profile(ti_report, output_dir=outdir)
        for r in ti_report.point_reports:
            plot_point_diagnostics(r, output_dir=outdir)

        # ===== Phase 2: Constant-potential correction =====

        # 7. Load calibration
        load_calibration_json = lazy_import(
            "md_analysis.electrochemical.calibration._data",
            "load_calibration_json",
        )
        mapper_from_dict = lazy_import(
            "md_analysis.electrochemical.calibration._mapper",
            "mapper_from_dict",
        )
        cal_path = ctx.get(K.CALIBRATION_JSON)
        _cal_data, fit_params = load_calibration_json(cal_path)
        mapper = mapper_from_dict(fit_params)

        # 8. Compute correction
        result = compute_constant_potential_correction(
            ti_report,
            point_defs,
            mapper,
            target_side=ctx[K.TARGET_SIDE],
            method=ctx[K.METHOD],
            normal=ctx[K.NORMAL],
        )

        # 9. Correction console summary
        corr = result.correction
        print(f"\n  --- Constant-Potential Correction (Norskov) ---")
        print(f"  Electrode surface: {ctx[K.TARGET_SIDE]}, area = {corr.area_A2:.2f} Å²")
        print(f"\n  {'Point':<6} {'ξ':<12} {'σ(μC/cm²)':<14} {'Φ(V/SHE)':<14} {'correction(eV)'}")
        print(f"  {'─' * 60}")
        for i in range(len(xi_values)):
            print(
                f"  {i:<6} {xi_values[i]:<12.6f} "
                f"{corr.sigma_uC_cm2[i]:<14.4f} {corr.phi_V_SHE[i]:<14.4f} "
                f"{corr.correction_eV[i]:<14.6f}"
            )
        total_corr = result.delta_A_const_phi_eV - float(result.A_const_q_eV[-1])
        print(f"\n  ΔA (const-Φ) = {result.delta_A_const_phi_eV:.6f} eV"
              f"  (correction = {total_corr:+.6f} eV)")

        # 10. Write corrected output files
        write_corrected_free_energy_csv(result, output_dir=outdir)
        plot_corrected_free_energy_profile(result, output_dir=outdir)

        # 11. List all output files
        print(f"\n  Output files:")
        print(f"    {outdir / 'ti_free_energy.png'}")
        print(f"    {outdir / 'ti_free_energy.csv'}")
        print(f"    {outdir / 'ti_corrected_free_energy.png'}")
        print(f"    {outdir / 'ti_corrected_free_energy.csv'}")
        print(f"    {outdir / 'ti_convergence_report.csv'}")
        for r in ti_report.point_reports:
            print(f"    {outdir / f'ti_diag_xi{r.xi:.4f}.png'}")
