"""Programmatic entry points for running analysis workflows.

Public API
----------
- ``run_water_analysis``  — water density/orientation/adsorbed-layer analysis
- ``run_water_analysis_with_report``  — agent-facing water wrapper
- ``run_potential_analysis`` — Hartree potential / Fermi / electrode potential analysis
- ``run_potential_analysis_with_report`` — agent-facing potential wrapper
- ``run_charge_analysis`` — Bader surface charge density time series
- ``run_charge_analysis_with_report``          — agent-facing surface charge wrapper
- ``run_tracked_charge_analysis_with_report``  — agent-facing tracked-atom charge wrapper
- ``run_counterion_charge_analysis_with_report`` — agent-facing counterion charge wrapper
- ``run_all`` — both water + potential
- ``run_all_with_report`` — agent-facing composite wrapper
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .electrochemical.potential.config import DEFAULT_THICKNESS_ANG
from .utils.constants import CHARGE_METHOD_COUNTERION, DEFAULT_LAYER_TOL_A

logger = logging.getLogger(__name__)


def run_water_analysis(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run water analysis (three-panel plot + CSVs).

    All outputs are written directly into *output_dir*. Callers that want
    the canonical ``<root>/water/`` layout should pass
    ``output_dir=<root> / "water"`` (``run_all`` does this automatically).

    Returns a dict mapping output names to file paths.
    """
    logger.info("Starting water analysis: xyz=%s, output_dir=%s", xyz_path, output_dir)

    from .water import plot_water_three_panel_analysis
    from .water.config import (
        DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
        DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME,
    )

    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    png_path = plot_water_three_panel_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        cell_abc=cell_abc,
        output_dir=water_dir,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **kwargs,
    )

    return {
        "density_csv": water_dir / DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        "orientation_csv": water_dir / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        "adsorbed_profile_csv": water_dir / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        "adsorbed_range_txt": water_dir / DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        "adsorbed_theta_csv": water_dir / DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
        "plot_png": png_path,
    }


def run_potential_analysis(
    *,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    compute_u: bool = True,
    compute_phi_z: bool = True,
    max_curves: int = 0,
    thickness_end: float = 15.0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> dict[str, Path]:
    """Run all potential analysis workflows.

    Per-sub-analysis outputs are written into ``output_dir/<sub>/`` where
    ``<sub>`` is ``center``/``fermi``/``electrode``/``phi_z``/``thickness_sensitivity``.
    Callers that want the canonical ``<root>/electrochemical/potential/``
    layout should pass ``output_dir=<root> / "electrochemical" / "potential"``
    (``run_all`` does this automatically).

    Returns a dict mapping output names to file paths.
    """
    logger.info("Starting potential analysis: output_dir=%s", output_dir)

    from .electrochemical.potential import (
        center_slab_potential_analysis,
        fermi_energy_analysis,
        electrode_potential_analysis,
        phi_z_planeavg_analysis,
        thickness_sensitivity_analysis,
    )

    pot_dir = Path(output_dir)
    results: dict[str, Path] = {}

    # Shared distributed-mode kwargs
    _dist = {
        "input_mode": input_mode,
        "sp_root_dir": sp_root_dir,
        "sp_dir_pattern": sp_dir_pattern,
        "sp_cube_filename": sp_cube_filename,
        "sp_out_filename": sp_out_filename,
    }
    is_distributed = input_mode == "distributed"

    # In distributed mode, Fermi data comes from sp.out in each subdir
    has_fermi = is_distributed or md_out_path is not None

    if compute_u and has_fermi:
        # Full electrode potential analysis (includes center + fermi)
        electrode_dir = pot_dir / "electrode"
        electrode_dir.mkdir(parents=True, exist_ok=True)
        u_csv = electrode_potential_analysis(
            cube_pattern,
            md_out_path,
            output_dir=electrode_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        results["electrode_csv"] = u_csv
    else:
        # Run individually
        center_dir = pot_dir / "center"
        center_dir.mkdir(parents=True, exist_ok=True)
        center_csv = center_slab_potential_analysis(
            cube_pattern,
            output_dir=center_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        results["center_csv"] = center_csv

        if md_out_path is not None:
            fermi_dir = pot_dir / "fermi"
            fermi_dir.mkdir(parents=True, exist_ok=True)
            fermi_csv = fermi_energy_analysis(
                md_out_path,
                output_dir=fermi_dir,
                fermi_unit=fermi_unit,
                frame_start=frame_start,
                frame_end=frame_end,
                frame_step=frame_step,
                **_dist,
            )
            results["fermi_csv"] = fermi_csv

    if compute_phi_z:
        phi_z_dir = pot_dir / "phi_z"
        phi_z_dir.mkdir(parents=True, exist_ok=True)
        phi_z_png = phi_z_planeavg_analysis(
            cube_pattern,
            output_dir=phi_z_dir,
            max_curves=max_curves,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        results["phi_z_png"] = phi_z_png

    # Thickness sensitivity sweep
    if has_fermi:
        ts_dir = pot_dir / "thickness_sensitivity"
        ts_dir.mkdir(parents=True, exist_ok=True)
        ts_csv = thickness_sensitivity_analysis(
            cube_pattern,
            md_out_path,
            output_dir=ts_dir,
            thickness_end=thickness_end,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        results["thickness_sensitivity_csv"] = ts_csv

    return results


def run_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = CHARGE_METHOD_COUNTERION,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = 1,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Run surface charge density analysis (CSV + PNG).

    Outputs are written into ``output_dir/<method>/`` (the method sub-dir
    is mandatory so aligned/counterion/layer runs don't collide). Callers
    that want the canonical ``<root>/electrochemical/charge/<method>/``
    layout should pass ``output_dir=<root> / "electrochemical" / "charge"``
    (``run_all`` does this automatically).

    Returns a dict mapping output names to file paths.
    """
    logger.info("Starting charge analysis: method=%s, output_dir=%s", method, output_dir)

    from .electrochemical.charge import surface_charge_analysis
    from .electrochemical.charge.config import (
        DEFAULT_SURFACE_CHARGE_CSV_NAME,
        DEFAULT_SURFACE_CHARGE_PNG_NAME,
    )

    charge_dir = Path(output_dir) / method
    charge_dir.mkdir(parents=True, exist_ok=True)

    result = surface_charge_analysis(
        root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        method=method,
        layer_tol_A=layer_tol_A,
        n_surface_layers=n_surface_layers,
        dir_pattern=dir_pattern,
        output_dir=charge_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return {
        "charge_csv": result.csv_path,
        "charge_png": result.csv_path.parent / DEFAULT_SURFACE_CHARGE_PNG_NAME,
    }


def run_tracked_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    atom_indices_xyz: Iterable[int],
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Track Bader net charges for specified XYZ atoms (CSV + PNG).

    Outputs are written into ``output_dir/tracked/``. Callers that want
    the canonical ``<root>/electrochemical/charge/tracked/`` layout should
    pass ``output_dir=<root> / "electrochemical" / "charge"`` (``run_all``
    does this automatically).

    Returns a dict mapping output names to file paths.
    """
    logger.info("Starting tracked charge analysis: output_dir=%s", output_dir)

    from .electrochemical.charge import tracked_atom_charge_analysis
    from .electrochemical.charge.Bader.AtomCharges import (
        DEFAULT_TRACKED_CHARGE_CSV,
        DEFAULT_TRACKED_CHARGE_PNG,
    )

    tracked_dir = Path(output_dir) / "tracked"
    tracked_dir.mkdir(parents=True, exist_ok=True)

    csv_path = tracked_atom_charge_analysis(
        root_dir,
        atom_indices_xyz=atom_indices_xyz,
        dir_pattern=dir_pattern,
        output_dir=tracked_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return {
        "tracked_charge_csv": csv_path,
        "tracked_charge_png": csv_path.parent / DEFAULT_TRACKED_CHARGE_PNG,
    }


def run_counterion_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Detect counterions per-frame and track their Bader charges (CSV + PNG).

    Outputs are written into ``output_dir/counterion_tracking/``. Callers
    that want the canonical ``<root>/electrochemical/charge/counterion_tracking/``
    layout should pass ``output_dir=<root> / "electrochemical" / "charge"``
    (``run_all`` does this automatically).

    Returns a dict mapping output names to file paths.
    """
    logger.info("Starting counterion charge analysis: output_dir=%s", output_dir)

    from .electrochemical.charge import counterion_charge_analysis
    from .electrochemical.charge.Bader.AtomCharges import (
        DEFAULT_COUNTERION_CHARGE_CSV,
        DEFAULT_COUNTERION_CHARGE_PNG,
        DEFAULT_COUNTERION_SUMMARY_CSV,
    )

    ci_dir = Path(output_dir) / "counterion_tracking"
    ci_dir.mkdir(parents=True, exist_ok=True)

    csv_path = counterion_charge_analysis(
        root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        layer_tol_A=layer_tol_A,
        dir_pattern=dir_pattern,
        output_dir=ci_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return {
        "counterion_charge_csv": csv_path,
        "counterion_summary_csv": csv_path.parent / DEFAULT_COUNTERION_SUMMARY_CSV,
        "counterion_charge_png": csv_path.parent / DEFAULT_COUNTERION_CHARGE_PNG,
    }


# ---------------------------------------------------------------------------
# Agent-facing charge wrappers: structured reports
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ChargeSurfaceReport:
    """Return value of :func:`run_charge_analysis_with_report`."""

    charge_csv: Path
    charge_png: Path
    n_frames: int
    method: str
    sigma_aligned_mean: float
    sigma_aligned_std: float
    sigma_opposed_mean: float
    sigma_opposed_std: float
    phi_cumavg_last: float | None
    phi_reference: str | None
    output_dir: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "charge_csv": str(self.charge_csv),
            "charge_png": str(self.charge_png),
            "n_frames": int(self.n_frames),
            "method": str(self.method),
            "sigma_aligned_mean": float(self.sigma_aligned_mean),
            "sigma_aligned_std": float(self.sigma_aligned_std),
            "sigma_opposed_mean": float(self.sigma_opposed_mean),
            "sigma_opposed_std": float(self.sigma_opposed_std),
            "phi_cumavg_last": (
                float(self.phi_cumavg_last)
                if self.phi_cumavg_last is not None else None
            ),
            "phi_reference": (
                str(self.phi_reference)
                if self.phi_reference is not None else None
            ),
            "output_dir": str(self.output_dir),
        }


@dataclass(frozen=True)
class TrackedChargeReport:
    """Return value of :func:`run_tracked_charge_analysis_with_report`."""

    tracked_charge_csv: Path
    tracked_charge_png: Path
    n_frames: int
    atom_indices_xyz: tuple[int, ...]
    n_atoms_tracked: int
    output_dir: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "tracked_charge_csv": str(self.tracked_charge_csv),
            "tracked_charge_png": str(self.tracked_charge_png),
            "n_frames": int(self.n_frames),
            "atom_indices_xyz": [int(i) for i in self.atom_indices_xyz],
            "n_atoms_tracked": int(self.n_atoms_tracked),
            "output_dir": str(self.output_dir),
        }


@dataclass(frozen=True)
class CounterionChargeReport:
    """Return value of :func:`run_counterion_charge_analysis_with_report`.

    ``counterion_charge_png`` is ``None`` when the underlying detection
    produced zero counterions (:func:`plot_counterion_charges`
    short-circuits in that case).  The agent layer must treat it as
    an optional artifact and only surface it when non-``None``.
    """

    counterion_charge_csv: Path
    counterion_summary_csv: Path
    counterion_charge_png: Path | None
    n_frames: int
    n_unique_counterions: int
    output_dir: Path

    def to_dict(self) -> dict[str, object]:
        return {
            "counterion_charge_csv": str(self.counterion_charge_csv),
            "counterion_summary_csv": str(self.counterion_summary_csv),
            "counterion_charge_png": (
                str(self.counterion_charge_png)
                if self.counterion_charge_png is not None else None
            ),
            "n_frames": int(self.n_frames),
            "n_unique_counterions": int(self.n_unique_counterions),
            "output_dir": str(self.output_dir),
        }


def run_charge_analysis_with_report(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = CHARGE_METHOD_COUNTERION,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = 1,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> ChargeSurfaceReport:
    """Agent-safe wrapper around :func:`run_charge_analysis`.

    Preserves the legacy ``output_dir / method`` layout and reads
    metrics straight from the :class:`SurfaceChargeResult` returned by
    :func:`surface_charge_analysis`, so the agent never re-parses the
    generated CSV.  This wrapper **only analyses existing Bader frame
    directories** — it does not generate Bader dirs, submit VASP jobs,
    or decide whether constant-potential correction is needed.  The
    optional calibration-potential append is existing reporting
    behaviour (it writes a ``phi_*`` column when a calibration JSON
    exists at the default location); missing calibration is logged and
    ignored.
    """
    from .electrochemical.charge import surface_charge_analysis
    from .electrochemical.charge.config import DEFAULT_SURFACE_CHARGE_PNG_NAME

    charge_dir = Path(output_dir) / method
    charge_dir.mkdir(parents=True, exist_ok=True)

    result = surface_charge_analysis(
        root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        method=method,
        layer_tol_A=layer_tol_A,
        n_surface_layers=n_surface_layers,
        dir_pattern=dir_pattern,
        output_dir=charge_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    png_path = result.csv_path.parent / DEFAULT_SURFACE_CHARGE_PNG_NAME

    return ChargeSurfaceReport(
        charge_csv=result.csv_path,
        charge_png=png_path,
        n_frames=int(result.n_frames),
        method=str(method),
        sigma_aligned_mean=float(result.sigma_aligned_mean),
        sigma_aligned_std=float(result.sigma_aligned_std),
        sigma_opposed_mean=float(result.sigma_opposed_mean),
        sigma_opposed_std=float(result.sigma_opposed_std),
        phi_cumavg_last=(
            float(result.phi_cumavg_last)
            if result.phi_cumavg_last is not None else None
        ),
        phi_reference=(
            str(result.phi_reference)
            if result.phi_reference is not None else None
        ),
        output_dir=charge_dir,
    )


def run_tracked_charge_analysis_with_report(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    atom_indices_xyz: Iterable[int],
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> TrackedChargeReport:
    """Agent-safe wrapper around :func:`run_tracked_charge_analysis`.

    Pulls ``n_frames`` / ``atom_indices_xyz`` / ``n_atoms_tracked`` from
    in-memory metadata via
    :func:`tracked_atom_charge_analysis_with_report` so the agent layer
    never re-parses the generated CSV.  **Analysis only** — does not
    generate Bader dirs or submit jobs.
    """
    from .electrochemical.charge.Bader.AtomCharges import (
        tracked_atom_charge_analysis_with_report,
    )

    tracked_dir = Path(output_dir) / "tracked"
    tracked_dir.mkdir(parents=True, exist_ok=True)

    result = tracked_atom_charge_analysis_with_report(
        root_dir,
        atom_indices_xyz=atom_indices_xyz,
        dir_pattern=dir_pattern,
        output_dir=tracked_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return TrackedChargeReport(
        tracked_charge_csv=result.csv_path,
        tracked_charge_png=result.png_path,
        n_frames=int(result.n_frames),
        atom_indices_xyz=tuple(int(i) for i in result.atom_indices_xyz),
        n_atoms_tracked=int(result.n_atoms_tracked),
        output_dir=tracked_dir,
    )


def run_counterion_charge_analysis_with_report(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> CounterionChargeReport:
    """Agent-safe wrapper around :func:`run_counterion_charge_analysis`.

    Uses :func:`counterion_charge_analysis_with_report` so ``n_frames``
    and ``n_unique_counterions`` come from in-memory ``frame_records``
    rather than reparsing the generated summary CSV.  **Analysis only**
    — does not generate Bader dirs, submit jobs, or decide whether
    detected ions are chemically relevant.
    """
    from .electrochemical.charge.Bader.AtomCharges import (
        counterion_charge_analysis_with_report,
    )

    ci_dir = Path(output_dir) / "counterion_tracking"
    ci_dir.mkdir(parents=True, exist_ok=True)

    result = counterion_charge_analysis_with_report(
        root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        layer_tol_A=layer_tol_A,
        dir_pattern=dir_pattern,
        output_dir=ci_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return CounterionChargeReport(
        counterion_charge_csv=result.csv_path,
        counterion_summary_csv=result.summary_path,
        counterion_charge_png=result.png_path,  # may be None when no detection
        n_frames=int(result.n_frames),
        n_unique_counterions=int(result.n_unique_counterions),
        output_dir=ci_dir,
    )


# ---------------------------------------------------------------------------
# Agent-facing water / potential wrappers: structured reports
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WaterThreePanelReport:
    """Return value of :func:`run_water_analysis_with_report`.

    ``artifacts`` carries one :class:`Path` per expected water output
    key (``density_csv``, ``orientation_csv``, ``adsorbed_profile_csv``,
    ``adsorbed_range_txt``, ``adsorbed_theta_csv``, ``plot_png``).
    The wrapper verifies every path exists on disk before returning,
    so downstream agents never see a dangling artifact.
    """

    artifacts: dict[str, Path]
    output_dir: Path
    n_artifacts: int
    frame_start: int | None
    frame_end: int | None
    frame_step: int | None

    def to_dict(self) -> dict[str, object]:
        return {
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "output_dir": str(self.output_dir),
            "n_artifacts": int(self.n_artifacts),
            "frame_start": (
                int(self.frame_start)
                if self.frame_start is not None else None
            ),
            "frame_end": (
                int(self.frame_end)
                if self.frame_end is not None else None
            ),
            "frame_step": (
                int(self.frame_step)
                if self.frame_step is not None else None
            ),
        }


@dataclass(frozen=True)
class PotentialFullReport:
    """Return value of :func:`run_potential_analysis_with_report`.

    Sub-analyses are conditional on inputs (Fermi availability,
    ``compute_u`` / ``compute_phi_z`` flags, ``input_mode``).  The
    boolean ``ran_*`` fields reflect which artifacts were **actually**
    produced and verified on disk, not merely which flags the caller
    requested.  ``artifacts`` therefore only contains existing files.
    """

    artifacts: dict[str, Path]
    output_dir: Path
    input_mode: str
    has_fermi: bool
    ran_electrode: bool
    ran_center: bool
    ran_fermi: bool
    ran_phi_z: bool
    ran_thickness_sensitivity: bool
    n_artifacts: int

    def to_dict(self) -> dict[str, object]:
        return {
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "output_dir": str(self.output_dir),
            "input_mode": str(self.input_mode),
            "has_fermi": bool(self.has_fermi),
            "ran_electrode": bool(self.ran_electrode),
            "ran_center": bool(self.ran_center),
            "ran_fermi": bool(self.ran_fermi),
            "ran_phi_z": bool(self.ran_phi_z),
            "ran_thickness_sensitivity": bool(self.ran_thickness_sensitivity),
            "n_artifacts": int(self.n_artifacts),
        }


_WATER_EXPECTED_KEYS = (
    "density_csv",
    "orientation_csv",
    "adsorbed_profile_csv",
    "adsorbed_range_txt",
    "adsorbed_theta_csv",
    "plot_png",
)


def run_water_analysis_with_report(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> WaterThreePanelReport:
    """Agent-safe wrapper around :func:`run_water_analysis`.

    Preserves every artifact filename and layout; the only difference
    from the legacy function is the structured report + a strict
    post-condition that **every** expected artifact exists on disk
    before the wrapper returns.  A missing artifact after a successful
    lower-level call is surfaced as :class:`RuntimeError` (mapped to
    ``analysis``) rather than as a success with a dangling output
    path.

    This wrapper does not parse generated CSVs for metrics.
    """
    result = run_water_analysis(
        xyz_path,
        md_inp_path,
        cell_abc=cell_abc,
        output_dir=output_dir,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **kwargs,
    )

    artifacts: dict[str, Path] = {}
    missing: list[str] = []
    for key in _WATER_EXPECTED_KEYS:
        path = result.get(key)
        if path is None:
            missing.append(key)
            continue
        p = Path(path)
        if not p.is_file():
            missing.append(f"{key}={p}")
            continue
        artifacts[key] = p

    if missing:
        raise RuntimeError(
            "Water three-panel analysis did not produce expected artifacts: "
            + ", ".join(missing)
        )

    return WaterThreePanelReport(
        artifacts=artifacts,
        output_dir=Path(output_dir),
        n_artifacts=len(artifacts),
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )


_POTENTIAL_EXPECTED_KEYS = (
    "electrode_csv",
    "center_csv",
    "fermi_csv",
    "phi_z_png",
    "thickness_sensitivity_csv",
)


def run_potential_analysis_with_report(
    *,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    compute_u: bool = True,
    compute_phi_z: bool = True,
    max_curves: int = 0,
    thickness_end: float = 15.0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> PotentialFullReport:
    """Agent-safe wrapper around :func:`run_potential_analysis`.

    Sub-analyses run conditionally based on inputs and flags.  The
    wrapper filters the returned artifact dict to only files that
    actually exist on disk and derives the ``ran_*`` booleans from
    which artifacts were produced, not from which flags were set —
    this is important because silent no-ops (e.g. Fermi energy
    unavailable) otherwise leak as "ran" signals to the agent.

    ``has_fermi`` reflects whether any Fermi-dependent artifact
    (electrode / fermi / thickness sensitivity) was actually produced.

    **Continuous mode cwd semantics preserved**: ``cube_pattern`` is
    resolved against the process current working directory, matching
    the legacy :func:`run_potential_analysis` behaviour.  The wrapper
    does **not** introduce a ``workdir`` argument.

    This wrapper does not parse generated CSVs for metrics.
    """
    # Boundary validation — legacy run_potential_analysis silently
    # treats unknown input_mode as non-distributed (i.e. continuous),
    # which would surface as FileNotFoundError later.  For agent calls
    # we want the cleaner ``validation`` classification.
    if input_mode not in ("continuous", "distributed"):
        raise ValueError(
            f"input_mode must be 'continuous' or 'distributed', got "
            f"{input_mode!r}"
        )
    if center_mode not in ("interface", "cell"):
        raise ValueError(
            f"center_mode must be 'interface' or 'cell', got "
            f"{center_mode!r}"
        )
    if fermi_unit not in ("au", "ev"):
        raise ValueError(
            f"fermi_unit must be 'au' or 'ev', got {fermi_unit!r}"
        )

    result = run_potential_analysis(
        output_dir=output_dir,
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        thickness_ang=thickness_ang,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        fermi_unit=fermi_unit,
        compute_u=compute_u,
        compute_phi_z=compute_phi_z,
        max_curves=max_curves,
        thickness_end=thickness_end,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        input_mode=input_mode,
        sp_root_dir=sp_root_dir,
        sp_dir_pattern=sp_dir_pattern,
        sp_cube_filename=sp_cube_filename,
        sp_out_filename=sp_out_filename,
    )

    artifacts: dict[str, Path] = {}
    for key in _POTENTIAL_EXPECTED_KEYS:
        path = result.get(key)
        if path is None:
            continue
        p = Path(path)
        if p.is_file():
            artifacts[key] = p

    ran_electrode = "electrode_csv" in artifacts
    ran_center = "center_csv" in artifacts
    ran_fermi = "fermi_csv" in artifacts
    ran_phi_z = "phi_z_png" in artifacts
    ran_thickness_sensitivity = "thickness_sensitivity_csv" in artifacts
    has_fermi = ran_electrode or ran_fermi or ran_thickness_sensitivity

    return PotentialFullReport(
        artifacts=artifacts,
        output_dir=Path(output_dir),
        input_mode=str(input_mode),
        has_fermi=has_fermi,
        ran_electrode=ran_electrode,
        ran_center=ran_center,
        ran_fermi=ran_fermi,
        ran_phi_z=ran_phi_z,
        ran_thickness_sensitivity=ran_thickness_sensitivity,
        n_artifacts=len(artifacts),
    )


def run_all(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run all analysis workflows (water + potential).

    Returns a merged dict of all output file paths.
    """
    logger.info("Starting full analysis: output_dir=%s", output_dir)

    root = Path(output_dir)
    results: dict[str, Path] = {}

    results.update(run_water_analysis(
        xyz_path, md_inp_path, cell_abc=cell_abc,
        output_dir=root / "water",
        frame_start=frame_start, frame_end=frame_end, frame_step=frame_step,
        verbose=verbose,
    ))

    pot_kwargs = {
        k: v for k, v in kwargs.items()
        if k in {
            "thickness_ang", "center_mode", "metal_elements",
            "layer_tol_ang", "fermi_unit", "compute_u",
            "compute_phi_z", "max_curves", "thickness_end",
        }
    }
    results.update(run_potential_analysis(
        output_dir=root / "electrochemical" / "potential",
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **pot_kwargs,
    ))

    return results


# ---------------------------------------------------------------------------
# Agent-facing composite wrapper — water + potential only.
# ---------------------------------------------------------------------------


# Whitelist of potential kwargs `run_all()` forwards to
# ``run_potential_analysis``; mirrored here so the composite wrapper
# keeps the exact same forwarding behaviour (no new surface area).
_RUN_ALL_POTENTIAL_KWARGS = frozenset({
    "thickness_ang", "center_mode", "metal_elements",
    "layer_tol_ang", "fermi_unit", "compute_u",
    "compute_phi_z", "max_curves", "thickness_end",
})


@dataclass(frozen=True)
class RunAllReport:
    """Return value of :func:`run_all_with_report`.

    ``run_all`` is a **composite convenience wrapper** — it just runs
    ``water_three_panel`` + ``potential_full`` together against the
    standard output layout.  It does NOT run charge / Bader /
    calibration / TI / constant-potential correction, and makes NO
    scientific judgements about sampling, convergence, or quality.

    Summary booleans / counts are all derived from real artifacts that
    actually exist on disk, not from the caller's requested flags.
    """

    artifacts: dict[str, Path]
    output_dir: Path
    water_output_dir: Path
    potential_output_dir: Path
    n_artifacts: int
    water_n_artifacts: int
    potential_n_artifacts: int
    ran_water: bool
    ran_potential: bool
    potential_has_fermi: bool
    potential_ran_electrode: bool
    potential_ran_center: bool
    potential_ran_fermi: bool
    potential_ran_phi_z: bool
    potential_ran_thickness_sensitivity: bool
    frame_start: int | None
    frame_end: int | None
    frame_step: int | None

    def to_dict(self) -> dict[str, object]:
        return {
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "output_dir": str(self.output_dir),
            "water_output_dir": str(self.water_output_dir),
            "potential_output_dir": str(self.potential_output_dir),
            "n_artifacts": int(self.n_artifacts),
            "water_n_artifacts": int(self.water_n_artifacts),
            "potential_n_artifacts": int(self.potential_n_artifacts),
            "ran_water": bool(self.ran_water),
            "ran_potential": bool(self.ran_potential),
            "potential_has_fermi": bool(self.potential_has_fermi),
            "potential_ran_electrode": bool(self.potential_ran_electrode),
            "potential_ran_center": bool(self.potential_ran_center),
            "potential_ran_fermi": bool(self.potential_ran_fermi),
            "potential_ran_phi_z": bool(self.potential_ran_phi_z),
            "potential_ran_thickness_sensitivity": bool(
                self.potential_ran_thickness_sensitivity
            ),
            "frame_start": (
                int(self.frame_start) if self.frame_start is not None else None
            ),
            "frame_end": (
                int(self.frame_end) if self.frame_end is not None else None
            ),
            "frame_step": (
                int(self.frame_step) if self.frame_step is not None else None
            ),
        }


def run_all_with_report(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> RunAllReport:
    """Agent-safe composite wrapper: water + potential.

    Reuses the Batch 4 leaf wrappers so all artifact verification and
    metric extraction lives in one place per leaf:

    - :func:`run_water_analysis_with_report`
    - :func:`run_potential_analysis_with_report`

    Layout is preserved:

    - water  → ``output_dir / "water"``
    - potential → ``output_dir / "electrochemical" / "potential"``

    Potential kwargs follow the same whitelist as legacy
    :func:`run_all` (``_RUN_ALL_POTENTIAL_KWARGS``); anything else in
    ``**kwargs`` is silently dropped to match historical behaviour.

    Raises
    ------
    RuntimeError
        If water and potential produce an artifact under the same key
        (no collision exists today, but the composite refuses to
        silently overwrite instead of guessing).
    """
    root = Path(output_dir)

    water_report = run_water_analysis_with_report(
        xyz_path,
        md_inp_path,
        cell_abc=cell_abc,
        output_dir=root / "water",
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    pot_kwargs = {
        k: v for k, v in kwargs.items() if k in _RUN_ALL_POTENTIAL_KWARGS
    }
    potential_report = run_potential_analysis_with_report(
        output_dir=root / "electrochemical" / "potential",
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **pot_kwargs,
    )

    # Merge artifacts; refuse to silently clobber on a key collision.
    collisions = set(water_report.artifacts) & set(potential_report.artifacts)
    if collisions:
        raise RuntimeError(
            "run_all composite refuses to merge: water + potential "
            "reports share artifact keys "
            f"{sorted(collisions)!r}.  Fix upstream filename convention "
            "instead of silently overwriting one leaf's file."
        )

    merged: dict[str, Path] = {}
    merged.update(water_report.artifacts)
    merged.update(potential_report.artifacts)

    return RunAllReport(
        artifacts=merged,
        output_dir=root,
        water_output_dir=water_report.output_dir,
        potential_output_dir=potential_report.output_dir,
        n_artifacts=len(merged),
        water_n_artifacts=int(water_report.n_artifacts),
        potential_n_artifacts=int(potential_report.n_artifacts),
        ran_water=water_report.n_artifacts > 0,
        ran_potential=potential_report.n_artifacts > 0,
        potential_has_fermi=bool(potential_report.has_fermi),
        potential_ran_electrode=bool(potential_report.ran_electrode),
        potential_ran_center=bool(potential_report.ran_center),
        potential_ran_fermi=bool(potential_report.ran_fermi),
        potential_ran_phi_z=bool(potential_report.ran_phi_z),
        potential_ran_thickness_sensitivity=bool(
            potential_report.ran_thickness_sensitivity
        ),
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )
