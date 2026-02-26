"""Ensemble-average c profiles + auto-window theta-PDF plotting utility."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import iread

REPO_ROOT = Path(__file__).resolve().parents[4]

from src.structure.utils.WaterParser import _compute_water_mass_density_z_distribution as compute_water_mass_density_z_distribution
from src.structure.utils.WaterParser import _compute_water_orientation_theta_pdf_in_c_fraction_window as compute_water_orientation_theta_pdf_in_c_fraction_window
from src.structure.utils.WaterParser import _compute_water_orientation_weighted_density_z_distribution as compute_water_orientation_weighted_density_z_distribution
from src.structure.utils.WaterParser import detect_water_molecule_indices
from src.structure.utils.WaterParser import get_water_oxygen_indices_array
from src.structure.utils.config import DEFAULT_THETA_BIN_DEG
from src.structure.utils.config import DEFAULT_Z_BIN_WIDTH_A

ANGSTROM3_TO_CM3 = 1.0e-24


def _parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def _z_edges_and_centers_A(lz_A: float, dz_A: float) -> tuple[np.ndarray, np.ndarray]:
    if dz_A <= 0.0:
        raise ValueError("dz_A must be > 0")
    nbins = max(1, int(np.ceil(lz_A / dz_A)))
    if nbins == 1:
        edges = np.array([0.0, lz_A], dtype=float)
    else:
        interior = np.arange(1, nbins, dtype=float) * dz_A
        edges = np.concatenate(([0.0], interior, [lz_A]))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers


def _fraction_window_mask(values: np.ndarray, *, start: float, end: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    start_mod = float(np.mod(start, 1.0))
    end_mod = float(np.mod(end, 1.0))
    if np.isclose(start_mod, end_mod):
        return np.ones(values.shape, dtype=bool)
    if start_mod < end_mod:
        return (values >= start_mod) & (values < end_mod)
    return (values >= start_mod) | (values < end_mod)


def _find_first_two_near_zero_c_fractions(
    rho_profile: np.ndarray,
    c_fraction_centers: np.ndarray,
    *,
    near_zero_ratio: float = 0.05,
    min_separation_bins: int = 2,
) -> tuple[float, float]:
    """
    Find first two near-zero locations scanning from c-fraction 0 -> 1.

    near-zero is defined by threshold = near_zero_ratio * max(rho_profile),
    and candidates are local minima below that threshold.
    """
    rho = np.asarray(rho_profile, dtype=float).reshape(-1)
    cfrac = np.asarray(c_fraction_centers, dtype=float).reshape(-1)
    if rho.size != cfrac.size:
        raise ValueError("rho_profile and c_fraction_centers must have same length")
    if rho.size < 2:
        raise ValueError("rho_profile needs at least 2 bins")

    rho_smooth = rho.copy()
    if rho.size >= 5:
        kernel = np.array([1.0, 2.0, 3.0, 2.0, 1.0], dtype=float)
        kernel /= np.sum(kernel)
        rho_smooth = np.convolve(rho, kernel, mode="same")

    threshold = float(max(near_zero_ratio * float(np.max(rho_smooth)), 0.0))

    candidates: list[int] = []
    for i in range(rho_smooth.size):
        left = rho_smooth[i - 1] if i > 0 else rho_smooth[i]
        right = rho_smooth[i + 1] if i < rho_smooth.size - 1 else rho_smooth[i]
        if rho_smooth[i] <= left and rho_smooth[i] <= right and rho_smooth[i] <= threshold:
            candidates.append(i)

    selected: list[int] = []
    for idx in candidates:
        if not selected or (idx - selected[-1]) >= min_separation_bins:
            selected.append(idx)
        if len(selected) >= 2:
            break

    if len(selected) < 2:
        # Fallback: pick first two smallest values (ordered by c coordinate).
        smallest = np.argsort(rho_smooth)[: max(2, min(10, rho_smooth.size))]
        selected = sorted(int(i) for i in smallest)[:2]
        if len(selected) < 2:
            selected = [0, rho_smooth.size - 1]

    return float(cfrac[selected[0]]), float(cfrac[selected[1]])


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute ensemble-averaged rho(c), orientation(c), and theta-PDF in a c-fraction window. "
            "Use --c-start/--c-end for manual near-surface window; otherwise use auto near-zero detection."
        )
    )
    parser.add_argument("--c-start", type=float, default=None, help="manual c-fraction window start")
    parser.add_argument("--c-end", type=float, default=None, help="manual c-fraction window end")
    return parser.parse_args()


def _theta_bin_count_from_ndeg(ndeg: float) -> int:
    ndeg = float(ndeg)
    if ndeg <= 0.0:
        raise ValueError("ndeg must be > 0")
    n_bins_float = 180.0 / ndeg
    n_bins = int(round(n_bins_float))
    if not np.isclose(n_bins_float, float(n_bins), rtol=0.0, atol=1.0e-12):
        raise ValueError("theta bin width must divide 180 exactly")
    return n_bins


def main() -> None:
    args = _parse_args()
    if (args.c_start is None) != (args.c_end is None):
        raise ValueError("Please provide both --c-start and --c-end, or provide neither")

    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)
    z_edges_A, z_centers_A = _z_edges_and_centers_A(c_A, DEFAULT_Z_BIN_WIDTH_A)
    c_fraction_centers = z_centers_A / c_A

    bin_widths_A = np.diff(z_edges_A)
    area_xy_A2 = a_A * b_A
    bin_volumes_A3 = area_xy_A2 * bin_widths_A
    bin_volumes_cm3 = bin_volumes_A3 * ANGSTROM3_TO_CM3

    # Step 1: ensemble-averaged rho(c) and orientation-weighted profile along c.
    mass_sum_g = np.zeros_like(z_centers_A, dtype=float)
    cos_sum = np.zeros_like(z_centers_A, dtype=float)
    volume_time_sum_cm3 = np.zeros_like(z_centers_A, dtype=float)

    n_frames = 0
    n_waters_total = 0

    for atoms in iread(str(xyz_path), index=":"):
        atoms.set_cell([a_A, b_A, c_A])
        atoms.set_pbc([True, True, True])

        water_idx = detect_water_molecule_indices(atoms)
        oxygen_n1 = get_water_oxygen_indices_array(water_idx)

        rho = compute_water_mass_density_z_distribution(atoms, oxygen_n1, dz_A=DEFAULT_Z_BIN_WIDTH_A).reshape(-1)
        orient = compute_water_orientation_weighted_density_z_distribution(
            atoms,
            oxygen_n1,
            water_molecule_indices=water_idx,
            dz_A=DEFAULT_Z_BIN_WIDTH_A,
        ).reshape(-1)
        if rho.size != z_centers_A.size or orient.size != z_centers_A.size:
            raise ValueError("Profile bin size mismatch across frames")

        mass_sum_g += rho * bin_volumes_cm3
        cos_sum += orient * bin_volumes_A3
        volume_time_sum_cm3 += bin_volumes_cm3

        n_frames += 1
        n_waters_total += int(water_idx.shape[0])

    if n_frames == 0:
        raise RuntimeError("No frames were read from trajectory")

    rho_ensemble = mass_sum_g / volume_time_sum_cm3
    orient_ensemble = cos_sum / (volume_time_sum_cm3 / ANGSTROM3_TO_CM3)

    # Step 2: choose c-fraction window.
    if args.c_start is not None and args.c_end is not None:
        cfrac_start = float(args.c_start)
        cfrac_end = float(args.c_end)
        window_mode = "manual"
    else:
        cfrac_start, cfrac_end = _find_first_two_near_zero_c_fractions(rho_ensemble, c_fraction_centers)
        window_mode = "auto_near_zero"

    # Step 3: compute ensemble-averaged theta-PDF within auto window.
    n_theta_bins = _theta_bin_count_from_ndeg(DEFAULT_THETA_BIN_DEG)
    theta_pdf_weighted_sum = np.zeros(n_theta_bins, dtype=float)
    theta_weight_sum = 0.0

    for atoms in iread(str(xyz_path), index=":"):
        atoms.set_cell([a_A, b_A, c_A])
        atoms.set_pbc([True, True, True])

        water_idx = detect_water_molecule_indices(atoms)
        oxygen_n1 = get_water_oxygen_indices_array(water_idx)
        oxygen_flat = oxygen_n1.reshape(-1)

        scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
        oxygen_frac_c = scaled[oxygen_flat, 2]
        selected_mask = _fraction_window_mask(oxygen_frac_c, start=cfrac_start, end=cfrac_end)
        n_selected = int(np.sum(selected_mask))

        theta_pdf = compute_water_orientation_theta_pdf_in_c_fraction_window(
            atoms,
            oxygen_n1,
            [cfrac_start, cfrac_end],
            water_molecule_indices=water_idx,
            ndeg=DEFAULT_THETA_BIN_DEG,
        )

        if theta_pdf.size != n_theta_bins:
            raise ValueError("Theta PDF bin size mismatch")
        if n_selected > 0:
            theta_pdf_weighted_sum += theta_pdf * n_selected
            theta_weight_sum += float(n_selected)

    if theta_weight_sum > 0.0:
        theta_pdf_ensemble = theta_pdf_weighted_sum / theta_weight_sum
    else:
        theta_pdf_ensemble = np.zeros(n_theta_bins, dtype=float)
    avg_selected_waters_per_frame = theta_weight_sum / n_frames if n_frames > 0 else 0.0

    theta_edges = np.linspace(0.0, 180.0, n_theta_bins + 1, dtype=float)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])

    # Plot: 3 panels (rho(c), orient(c), theta-PDF in auto c-window).
    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=False)

    axes[0].plot(c_fraction_centers, rho_ensemble, color="tab:blue", lw=1.5)
    axes[0].axvline(cfrac_start, color="gray", ls="--", lw=1.0)
    axes[0].axvline(cfrac_end, color="gray", ls="--", lw=1.0)
    axes[0].set_xlim(0.0, 1.0)
    axes[0].set_ylabel("rho_water (g/cm^3)")
    axes[0].set_title("Ensemble-Averaged Water Density Along c")
    axes[0].grid(True, ls="--", alpha=0.3)

    axes[1].plot(c_fraction_centers, orient_ensemble, color="tab:orange", lw=1.5)
    axes[1].axvline(cfrac_start, color="gray", ls="--", lw=1.0)
    axes[1].axvline(cfrac_end, color="gray", ls="--", lw=1.0)
    axes[1].set_xlim(0.0, 1.0)
    axes[1].set_ylabel("orientation weighted (1/Angstrom^3)")
    axes[1].set_xlabel("c fractional coordinate")
    axes[1].set_title("Ensemble-Averaged Orientation-Weighted Profile Along c")
    axes[1].grid(True, ls="--", alpha=0.3)

    axes[2].plot(theta_centers, theta_pdf_ensemble, color="tab:green", lw=1.5)
    axes[2].set_xlim(0.0, 180.0)
    axes[2].set_xlabel("theta (degree)")
    axes[2].set_ylabel("PDF (degree^-1)")
    axes[2].set_title(
        f"Ensemble-Averaged Theta PDF in {window_mode} c Window [{cfrac_start:.4f}, {cfrac_end:.4f}]"
    )
    axes[2].grid(True, ls="--", alpha=0.3)

    fig.tight_layout()
    prefix = "all_frames_manual_c_window_theta_pdf" if window_mode == "manual" else "all_frames_auto_c_window_theta_pdf"
    png_path = out_dir / f"{prefix}.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    profile_csv_path = out_dir / "all_frames_ensemble_c_profiles.csv"
    profile_data = np.column_stack([c_fraction_centers, rho_ensemble, orient_ensemble])
    np.savetxt(
        profile_csv_path,
        profile_data,
        delimiter=",",
        header="c_fraction_center,rho_g_cm3,orientation_weighted_1_A3",
        comments="",
    )

    theta_csv_path = out_dir / f"{prefix}.csv"
    theta_data = np.column_stack([theta_centers, theta_pdf_ensemble])
    np.savetxt(theta_csv_path, theta_data, delimiter=",", header="theta_degree,pdf_degree_inv", comments="")

    meta_path = out_dir / f"{prefix}_meta.txt"
    meta_path.write_text(
        "\n".join(
            [
                f"window_mode={window_mode}",
                f"frames_processed={n_frames}",
                f"avg_water_molecules_per_frame={n_waters_total / n_frames:.6f}",
                f"window_c_fraction_start={cfrac_start:.10f}",
                f"window_c_fraction_end={cfrac_end:.10f}",
                f"avg_selected_waters_per_frame={avg_selected_waters_per_frame:.6f}",
                f"z_bin_width_A={DEFAULT_Z_BIN_WIDTH_A}",
                f"theta_bin_width_deg={DEFAULT_THETA_BIN_DEG}",
            ]
        ),
        encoding="utf-8",
    )

    print(f"frames processed: {n_frames}")
    print(f"average water molecules per frame: {n_waters_total / n_frames:.3f}")
    print(f"window mode: {window_mode}")
    print(f"c-fraction window: [{cfrac_start:.6f}, {cfrac_end:.6f}]")
    print(f"average selected waters per frame in window: {avg_selected_waters_per_frame:.3f}")
    print(f"saved figure: {png_path}")
    print(f"saved c-profile data: {profile_csv_path}")
    print(f"saved theta-pdf data: {theta_csv_path}")
    print(f"saved metadata: {meta_path}")


if __name__ == "__main__":
    main()
