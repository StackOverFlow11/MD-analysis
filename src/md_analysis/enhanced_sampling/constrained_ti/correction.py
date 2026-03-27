"""Constant-potential free energy correction for constrained TI.

Implements the Norskov parabolic correction that converts constant-charge
TI free energy to constant-potential free energy:

    dF_phi(xi) = dF_q(xi) + [sigma(xi) - sigma(IS)] * [Phi(xi) - Phi(IS)] * A / 2

where sigma is the ensemble-averaged surface charge density from Bader analysis,
Phi is the electrode potential predicted via sigma->phi calibration, and
A is the electrode surface area.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import DEFAULT_CORRECTED_FE_CSV_NAME, DEFAULT_CORRECTED_FE_PNG_NAME
from .models import TIPointDefinition, TIReport
from ...utils.config import AREA_VECTOR_INDICES, HA_TO_EV

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstantPotentialCorrection:
    """Per-point constant-potential correction data."""

    sigma_uC_cm2: np.ndarray  # (K,) ensemble-avg surface charge density
    phi_V_SHE: np.ndarray  # (K,) predicted electrode potential
    correction_eV: np.ndarray  # (K,) cumulative correction at each point
    area_A2: float  # electrode surface area (Angstrom^2)


@dataclass(frozen=True)
class ConstantPotentialResult:
    """Combined constant-charge TI + constant-potential correction."""

    ti_report: TIReport
    correction: ConstantPotentialCorrection
    A_const_q_eV: np.ndarray  # (K,) raw cumulative free energy (eV)
    A_const_phi_eV: np.ndarray  # (K,) corrected cumulative free energy (eV)
    delta_A_const_phi_eV: float  # final corrected free energy difference (eV)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_electrode_area(bader_dir: Path, normal: str) -> float:
    """Compute electrode surface area from first Bader frame's POSCAR."""
    from ...utils.BaderParser import load_bader_atoms
    from ...electrochemical.charge.config import (
        DEFAULT_ACF_FILENAME,
        DEFAULT_POTCAR_FILENAME,
        DEFAULT_STRUCTURE_FILENAME,
    )
    from ...electrochemical.charge.Bader._frame_utils import _sorted_frame_dirs

    frame_dirs = _sorted_frame_dirs(bader_dir)
    if not frame_dirs:
        raise FileNotFoundError(f"No bader frame directories in {bader_dir}")

    first_dir = frame_dirs[0]
    atoms = load_bader_atoms(
        first_dir / DEFAULT_STRUCTURE_FILENAME,
        first_dir / DEFAULT_ACF_FILENAME,
        first_dir / DEFAULT_POTCAR_FILENAME,
    )
    i0, i1 = AREA_VECTOR_INDICES[normal]
    cell = atoms.cell.array
    return float(np.linalg.norm(np.cross(cell[i0], cell[i1])))


_SIDE_INDEX = {"aligned": 0, "opposed": 1}


def _collect_bader_sigma(
    point_defs: list[TIPointDefinition],
    *,
    target_side: str,
    method: str = "counterion",
    normal: str = "c",
    **charge_kwargs: object,
) -> tuple[np.ndarray, float, list[bool]]:
    """Compute ensemble-averaged surface charge at each TI constraint point.

    Returns
    -------
    sigma_mean : np.ndarray, shape (K,)
        Ensemble-averaged sigma in uC/cm^2.
    area_A2 : float
        Electrode surface area in Angstrom^2.
    bader_available : list[bool]
        True for each point that has a bader/ directory.
    """
    from ...electrochemical.charge.Bader.SurfaceCharge import (
        trajectory_surface_charge,
    )

    if target_side not in _SIDE_INDEX:
        raise ValueError(
            f"target_side must be 'aligned' or 'opposed', got {target_side!r}"
        )
    side_idx = _SIDE_INDEX[target_side]

    k = len(point_defs)
    sigma_mean = np.zeros(k)
    bader_available: list[bool] = []
    area_A2 = 0.0

    for i, pdef in enumerate(point_defs):
        ti_dir = pdef.restart_path.parent
        bader_dir = ti_dir / "bader"

        if not bader_dir.is_dir():
            logger.warning(
                "bader/ not found for xi=%.6f (%s)", pdef.xi, ti_dir.name
            )
            bader_available.append(False)
            continue

        bader_available.append(True)

        # Compute area from first point that has bader data
        if area_A2 == 0.0:
            area_A2 = _get_electrode_area(bader_dir, normal)

        sigma_traj = trajectory_surface_charge(
            bader_dir,
            method=method,
            normal=normal,
            **charge_kwargs,
        )  # (t, 2) uC/cm^2
        sigma_mean[i] = float(np.mean(sigma_traj[:, side_idx]))

    return sigma_mean, area_A2, bader_available


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def compute_constant_potential_correction(
    ti_report: TIReport,
    point_defs: list[TIPointDefinition],
    mapper: object,
    *,
    target_side: str = "aligned",
    method: str = "counterion",
    normal: str = "c",
    **charge_kwargs: object,
) -> ConstantPotentialResult:
    """Compute constant-potential corrected free energy profile.

    Parameters
    ----------
    ti_report : TIReport
        Result of ``analyze_ti`` (constant-charge analysis).
    point_defs : list[TIPointDefinition]
        Discovered TI point definitions (same order as ti_report).
    mapper : ChargePotentialMapper
        Fitted sigma->phi mapper from calibration.
    target_side : str
        Electrode surface: ``"aligned"`` or ``"opposed"``.
    method : str
        Bader charge method (``"counterion"`` or ``"layer"``).
    normal : str
        Surface normal axis (``"a"``, ``"b"``, or ``"c"``).

    Returns
    -------
    ConstantPotentialResult
    """
    from ...electrochemical.charge.config import E_PER_A2_TO_UC_PER_CM2

    sigma_mean, area_A2, bader_available = _collect_bader_sigma(
        point_defs,
        target_side=target_side,
        method=method,
        normal=normal,
        **charge_kwargs,
    )

    # Cumulative raw free energy (eV)
    forces = ti_report.forces
    weights = ti_report.weights
    errors = ti_report.force_errors
    cumul_A_q = np.cumsum(weights * forces) * HA_TO_EV

    # If any point lacks bader data, return uncorrected
    if not all(bader_available):
        missing = [
            point_defs[i].xi
            for i, avail in enumerate(bader_available)
            if not avail
        ]
        logger.warning(
            "Bader data missing for %d point(s): %s. "
            "Skipping constant-potential correction.",
            len(missing),
            ", ".join(f"xi={x:.6f}" for x in missing),
        )
        zero_corr = ConstantPotentialCorrection(
            sigma_uC_cm2=sigma_mean,
            phi_V_SHE=np.zeros_like(sigma_mean),
            correction_eV=np.zeros_like(sigma_mean),
            area_A2=area_A2,
        )
        return ConstantPotentialResult(
            ti_report=ti_report,
            correction=zero_corr,
            A_const_q_eV=cumul_A_q,
            A_const_phi_eV=cumul_A_q.copy(),
            delta_A_const_phi_eV=float(cumul_A_q[-1]),
        )

    # Predict potential from sigma
    phi = mapper.predict(sigma_mean)  # (K,) V vs SHE

    # Norskov correction: [sigma(xi) - sigma(IS)] * [phi(xi) - phi(IS)] * A / 2
    sigma_IS = sigma_mean[0]
    phi_IS = phi[0]
    delta_sigma_e_A2 = (sigma_mean - sigma_IS) / E_PER_A2_TO_UC_PER_CM2
    delta_phi = phi - phi_IS
    correction_eV = delta_sigma_e_A2 * delta_phi * area_A2 / 2.0

    cumul_A_phi = cumul_A_q + correction_eV

    corr = ConstantPotentialCorrection(
        sigma_uC_cm2=sigma_mean,
        phi_V_SHE=phi,
        correction_eV=correction_eV,
        area_A2=area_A2,
    )
    return ConstantPotentialResult(
        ti_report=ti_report,
        correction=corr,
        A_const_q_eV=cumul_A_q,
        A_const_phi_eV=cumul_A_phi,
        delta_A_const_phi_eV=float(cumul_A_phi[-1]),
    )


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_corrected_free_energy_csv(
    result: ConstantPotentialResult,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Write corrected free-energy CSV."""
    from ...utils._io_helpers import _write_csv

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    out_path = out / DEFAULT_CORRECTED_FE_CSV_NAME

    tr = result.ti_report
    corr = result.correction
    cumul_sigma = np.sqrt(
        np.cumsum(tr.weights**2 * tr.force_errors**2)
    ) * HA_TO_EV

    fieldnames = [
        "xi", "weight", "dA_dxi", "sem",
        "A_const_q_eV", "sigma_A_cumulative_eV",
        "sigma_uC_cm2", "phi_V_SHE", "correction_eV", "A_const_phi_eV",
    ]
    rows = []
    for i in range(len(tr.xi_values)):
        rows.append({
            "xi": f"{tr.xi_values[i]:.8f}",
            "weight": f"{tr.weights[i]:.8f}",
            "dA_dxi": f"{tr.forces[i]:.8f}",
            "sem": f"{tr.force_errors[i]:.8f}",
            "A_const_q_eV": f"{result.A_const_q_eV[i]:.8f}",
            "sigma_A_cumulative_eV": f"{cumul_sigma[i]:.8f}",
            "sigma_uC_cm2": f"{corr.sigma_uC_cm2[i]:.6f}",
            "phi_V_SHE": f"{corr.phi_V_SHE[i]:.6f}",
            "correction_eV": f"{corr.correction_eV[i]:.8f}",
            "A_const_phi_eV": f"{result.A_const_phi_eV[i]:.8f}",
        })

    _write_csv(out_path, rows, fieldnames)
    return out_path


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------


def plot_corrected_free_energy_profile(
    result: ConstantPotentialResult,
    *,
    output_dir: Path | None = None,
) -> Path:
    """Plot free-energy profile with constant-potential correction overlay."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(output_dir) if output_dir else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    out_path = out / DEFAULT_CORRECTED_FE_PNG_NAME

    tr = result.ti_report
    corr = result.correction
    xi = tr.xi_values
    forces = tr.forces
    errors = tr.force_errors

    fig, ax1 = plt.subplots(figsize=(8, 5), dpi=180)

    # Left axis: dA/dxi with error bars (same as 312)
    colors = ["C2" if r.passed else "C3" for r in tr.point_reports]
    ax1.errorbar(
        xi, forces, yerr=errors, fmt="o", markersize=4,
        color="C0", ecolor="C0", capsize=3,
    )
    for i, (x, y) in enumerate(zip(xi, forces)):
        ax1.plot(x, y, "o", markersize=5, color=colors[i], zorder=5)
    ax1.set_xlabel("ξ (a.u.)")
    if len(xi) >= 2 and xi[0] > xi[-1]:
        ax1.invert_xaxis()
    ax1.set_ylabel("dA/dξ (a.u.)", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")

    # Right axis: cumulative free energy (eV)
    ax2 = ax1.twinx()
    cumul_sigma = np.sqrt(
        np.cumsum(tr.weights**2 * errors**2)
    ) * HA_TO_EV

    # Const-q (dashed orange with error band)
    ax2.plot(
        xi, result.A_const_q_eV, "--s", color="C1",
        markersize=3, linewidth=1.0, label="A(ξ) const-q", alpha=0.7,
    )
    ax2.fill_between(
        xi,
        result.A_const_q_eV - cumul_sigma,
        result.A_const_q_eV + cumul_sigma,
        alpha=0.15, color="C1",
    )

    # Const-phi (solid red)
    ax2.plot(
        xi, result.A_const_phi_eV, "-o", color="C3",
        markersize=3, linewidth=1.2, label="A(ξ) const-Φ",
    )

    ax2.set_ylabel("A(ξ) (eV)")
    ax2.legend(loc="best", fontsize=8)

    delta_q = float(result.A_const_q_eV[-1])
    sigma_q = cumul_sigma[-1]
    delta_phi = result.delta_A_const_phi_eV
    title = (
        f"ΔA(const-q) = {delta_q:.6f} ± {sigma_q:.6f} eV"
        f"  |  ΔA(const-Φ) = {delta_phi:.6f} eV"
    )
    ax1.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
