"""Constrained TI analysis command classes (311-313).

Phase 6.5: routed through ``workflows.enhanced_sampling.run_ti_*``
facades. The interactive prompts that don't fit a single end-to-end
call (point-slice selection, per-point equilibration override) stay
here: each command pre-discovers points with ``strict=False`` (the
historical lenient behavior on broken ``ti_target_*`` subdirs),
prompts, then dispatches to the workflow with ``strict=False`` so the
analysis pass sees the same point set the user saw.
"""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import K
from ._prompt import prompt_bool, prompt_choice, prompt_float, prompt_int, prompt_str
from ._enhanced_sampling import _discover_restart_file, _discover_log_file


_VALID_SIDES = ("aligned", "opposed")
_VALID_METHODS = ("counterion", "layer")


# ---------------------------------------------------------------------------
# Shared helpers for 312 and 313
# ---------------------------------------------------------------------------

def _collect_ti_base_params(ctx: dict) -> None:
    """Collect TI parameters shared by 312 and 313.

    Populates: TI_ROOT_DIR, EQUILIBRATION, EPSILON_TOL_EV, TI_REVERSE.
    Does NOT collect OUTDIR.
    """
    ctx[K.TI_ROOT_DIR] = prompt_str("TI root directory", default=".") or "."

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
    ctx[K.AUTO_EQUILIBRATION] = prompt_bool(
        "Auto-equilibration (iteratively discard first half until converged)?",
        default=False,
    )


def _prompt_ti_slice_and_equil(ctx: dict, point_defs):
    """Prompt for an optional Python-slice and per-point equilibration.

    ``point_slice`` is validated eagerly with the SAME parser the
    workflow uses (``_parse_point_slice``) so an ambiguous / illegal
    spec (e.g. ``"2"`` with no colon) raises ``ValueError`` BEFORE any
    per-point equilibration is collected. The string is forwarded
    verbatim to the workflow; the local slice is only used to size the
    per-point prompt and the selection echo.

    Returns
    -------
    (point_slice, equilibration)
        ``point_slice``: ``str | None`` (forwarded as-is to the workflow)
        ``equilibration``: ``int | list[int]`` (per-point list sized by
        the *sliced* point count so the workflow's re-discovered +
        re-sliced point set aligns 1:1)
    """
    _parse_point_slice = lazy_import(
        "md_analysis.enhanced_sampling.constrained_ti.workflow",
        "_parse_point_slice",
    )

    slice_str = prompt_str(
        "Select points (Python slice, e.g. 3:8, :8, 3::2, empty=all)",
        default="",
    )
    point_slice: str | None
    if slice_str:
        # Validate with the workflow's parser. Raises ValueError on
        # ambiguous / illegal specs; the MenuCommand.run() try-except
        # catches it and aborts before per-point prompts run.
        sl = _parse_point_slice(slice_str)
        sliced = point_defs[sl]
        point_slice = slice_str
        print(f"  Selected {len(sliced)} points:")
        for i, p in enumerate(sliced):
            print(f"    [{i}] ξ = {p.xi:.6f}")
    else:
        sliced = point_defs
        point_slice = None

    default_equil = ctx[K.EQUILIBRATION]
    if len(sliced) > 1 and prompt_bool(
        "Set per-point equilibration frames?", default=False,
    ):
        equil_list = []
        for p in sliced:
            val = prompt_int(
                f"  ξ={p.xi:.6f} equilibration frames",
                default=default_equil,
            )
            equil_list.append(val if val is not None else default_equil)
        equilibration: int | list[int] = equil_list
    else:
        equilibration = default_equil

    return point_slice, equilibration


def _print_ti_summary_table(ti_report) -> None:
    """Console convergence summary for a :class:`TIReport`.

    Shared by 312 and 313; reads ``ti_report.point_reports`` plus the
    pass/fail roll-up and the optional suggested time-allocation hint.
    """
    print(f"\n  {'Point':<6} {'ξ':<12} {'⟨λ⟩':<14} {'SEM':<14} "
          f"{'N':<8} {'Time range (fs)':<24} {'Status'}")
    print(f"  {'─' * 82}")
    for i, r in enumerate(ti_report.point_reports):
        status = "PASS" if r.passed else "FAIL"
        print(
            f"  {i:<6} {r.xi:<12.6f} {r.lambda_mean:<14.6f} "
            f"{r.sem_final:<14.6f} "
            f"{r.n_analyzed:<8} {r.time_start_fs:.1f} – {r.time_end_fs:.1f}"
            f"{'':>4}{status}"
        )

    if ti_report.all_passed:
        print("  Status: ALL PASS")
    else:
        failing = ti_report.failing_indices
        print(f"  Status: {len(failing)} FAILED "
              f"(indices: {', '.join(str(i) for i in failing)})")
        if ti_report.suggested_time_ratios is not None:
            xi_values = [r.xi for r in ti_report.point_reports]
            print("  Suggested time allocation (relative):")
            parts = [
                f"ξ={xi_values[i]:.4f}: "
                f"{ti_report.suggested_time_ratios[i]:.2f}"
                for i in range(len(xi_values))
            ]
            print(f"    {', '.join(parts)}")


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
        run_single = lazy_import(
            "md_analysis.workflows.enhanced_sampling",
            "run_ti_single_diagnostics",
        )
        outdir = ctx[K.OUTDIR_RESOLVED]
        result = run_single(
            restart_path=ctx[K.RESTART_PATH],
            log_path=ctx[K.LOG_PATH],
            output_dir=outdir,
            equilibration=ctx[K.EQUILIBRATION],
            sem_target=ctx[K.SEM_TARGET],
            colvar_id=ctx[K.COLVAR_ID],
        )

        # Console summary from result.extra (ConstraintPointReport) +
        # result.metadata (frame counts / time range).
        r = result.extra
        m = result.metadata
        status = "N/A" if r.passed is None else ("PASS" if r.passed else "FAIL")
        method = "plateau" if r.block_avg.plateau_reached else "acf"
        print(f"\n  ξ = {r.xi:.6f} a.u.")
        print(f"  Frames: {m['n_total']} total, {m['equilibration']} discarded, "
              f"{m['n_analyzed']} analyzed")
        print(f"  Time range: {m['time_start_fs']:.1f} – {m['time_end_fs']:.1f} fs")
        print(f"  ⟨λ⟩ = {r.lambda_mean:.6f} a.u.")
        print(f"  τ_corr = {r.autocorr.tau_corr:.1f} frames, N_eff = {r.autocorr.n_eff:.1f}")
        print(f"  SEM_final = {r.sem_final:.6f} a.u. (method: {method})")
        print(f"  Geweke |z| = {abs(r.geweke.z):.3f}  ({'PASS' if r.geweke.passed else 'FAIL'})")
        print(f"  Overall: {status}")

        if r.failure_reasons:
            for reason in r.failure_reasons:
                print(f"    - {reason}")

        for key in ("diagnostics_png", "csv"):
            if key in result.artifacts:
                print(f"  {key}: {result.artifacts[key]}")


# ---------------------------------------------------------------------------
# 312 — Full TI Analysis
# ---------------------------------------------------------------------------

class TIFullAnalysisCmd(MenuCommand):
    """Multi-point constrained TI convergence analysis + free-energy integration."""

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        _collect_ti_base_params(ctx)
        ctx[K.OUTDIR] = prompt_str("Output directory",
                                   default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        discover_ti_points = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.io",
            "discover_ti_points",
        )
        run_full = lazy_import(
            "md_analysis.workflows.enhanced_sampling",
            "run_ti_full_analysis",
        )
        outdir = ctx[K.OUTDIR_RESOLVED]
        root = Path(ctx[K.TI_ROOT_DIR])

        # Pre-discover for display + per-point equil prompt sizing.
        # strict=False is explicit (not the io default) so the wiring
        # test can pin the D4 invariant: CLI sees the same N points the
        # workflow analyses (workflow also passes strict=False).
        point_defs = discover_ti_points(
            root, reverse=ctx[K.TI_REVERSE], strict=False,
        )
        print(f"\n  Found {len(point_defs)} constraint points:")
        for i, p in enumerate(point_defs):
            print(f"    [{i}] ξ = {p.xi:.6f}")

        point_slice, equilibration = _prompt_ti_slice_and_equil(ctx, point_defs)

        result = run_full(
            root_dir=root,
            output_dir=outdir,
            reverse=ctx[K.TI_REVERSE],
            equilibration=equilibration,
            epsilon_tol_ev=ctx[K.EPSILON_TOL_EV],
            auto_equilibration=ctx.get(K.AUTO_EQUILIBRATION, False),
            point_slice=point_slice,
            strict=False,
        )

        _print_ti_summary_table(result.extra.ti_report)

        delta_A = result.metadata["delta_A_eV"]
        sigma_A = result.metadata["sigma_A_eV"]
        print(f"\n  ΔA = {delta_A:.6f} ± {sigma_A:.6f} eV")

        print(f"\n  Output files:")
        print(f"    {result.artifacts['free_energy_png']}")
        print(f"    {result.artifacts['free_energy_csv']}")
        print(f"    {result.artifacts['convergence_csv']}")
        for key, path in result.artifacts.items():
            if key.startswith("diagnostics_png_"):
                print(f"    {path}")


# ---------------------------------------------------------------------------
# 313 — Constant-Potential Free Energy Correction
# ---------------------------------------------------------------------------

class TIConstPotCorrectionCmd(MenuCommand):
    """Constant-potential free energy correction (Norskov) on top of TI analysis."""

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        _collect_ti_base_params(ctx)

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

        ctx[K.OUTDIR] = prompt_str("Output directory",
                                   default="analysis") or "analysis"
        return ctx

    def execute(self, ctx: dict) -> None:
        import numpy as np

        discover_ti_points = lazy_import(
            "md_analysis.enhanced_sampling.constrained_ti.io",
            "discover_ti_points",
        )
        run_corr = lazy_import(
            "md_analysis.workflows.enhanced_sampling",
            "run_ti_constant_potential_correction",
        )
        default_json = lazy_import(
            "md_analysis.electrochemical.calibration.config",
            "DEFAULT_CALIBRATION_FILE",
        )
        outdir = ctx[K.OUTDIR_RESOLVED]
        root = Path(ctx[K.TI_ROOT_DIR])

        # Pre-discover for display + per-point equil prompt. strict=False
        # explicit (D4 invariant — same as 312).
        point_defs = discover_ti_points(
            root, reverse=ctx[K.TI_REVERSE], strict=False,
        )
        print(f"\n  Found {len(point_defs)} constraint points:")
        for i, p in enumerate(point_defs):
            print(f"    [{i}] ξ = {p.xi:.6f}")

        point_slice, equilibration = _prompt_ti_slice_and_equil(ctx, point_defs)

        cal_path = ctx.get(K.CALIBRATION_JSON) or default_json

        result = run_corr(
            root_dir=root,
            output_dir=outdir,
            calibration_json_path=Path(cal_path),
            target_side=ctx[K.TARGET_SIDE],
            method=ctx[K.METHOD],
            normal=ctx[K.NORMAL],
            reverse=ctx[K.TI_REVERSE],
            equilibration=equilibration,
            epsilon_tol_ev=ctx[K.EPSILON_TOL_EV],
            auto_equilibration=ctx.get(K.AUTO_EQUILIBRATION, False),
            point_slice=point_slice,
            strict=False,
        )

        from ..utils.constants import HA_TO_EV

        correction_result = result.extra  # ConstantPotentialResult
        ti_report = correction_result.ti_report
        _print_ti_summary_table(ti_report)

        # Preserve the pre-migration numeric formula exactly: ΔA (const-q)
        # is the raw integrated TI free energy converted Hartree → eV
        # (NOT metadata["delta_A_const_q_eV"], whose definition is the
        # cumulative const-q value at the last point — same physical
        # quantity but a different code path; keep byte-equal output).
        delta_A_ev = ti_report.delta_A * HA_TO_EV
        sigma_A_ev = ti_report.sigma_A * HA_TO_EV
        print(f"\n  ΔA (const-q) = {delta_A_ev:.6f} ± {sigma_A_ev:.6f} eV")

        corr = correction_result.correction
        print(f"\n  --- Constant-Potential Correction (Norskov) ---")
        print(f"  Electrode surface: {ctx[K.TARGET_SIDE]}, "
              f"area = {corr.area_A2:.2f} Å²")
        print(f"\n  {'Point':<6} {'ξ':<12} {'σ(μC/cm²)':<14} "
              f"{'Φ(V/SHE)':<14} {'correction(eV)'}")
        print(f"  {'─' * 60}")
        xi_values = np.array([r.xi for r in ti_report.point_reports])
        for i in range(len(xi_values)):
            print(
                f"  {i:<6} {xi_values[i]:<12.6f} "
                f"{corr.sigma_uC_cm2[i]:<14.4f} "
                f"{corr.phi_V_SHE[i]:<14.4f} "
                f"{corr.correction_eV[i]:<14.6f}"
            )
        print(f"\n  ΔA (const-Φ) = {result.metadata['delta_A_const_phi_eV']:.6f} eV"
              f"  (correction = {result.metadata['total_correction_eV']:+.6f} eV)")

        # Output files
        print(f"\n  Output files:")
        for key in ("free_energy_png", "free_energy_csv",
                    "corrected_free_energy_png", "corrected_free_energy_csv",
                    "convergence_csv"):
            if key in result.artifacts:
                print(f"    {result.artifacts[key]}")
        for k, p in result.artifacts.items():
            if k.startswith("diagnostics_png_"):
                print(f"    {p}")
