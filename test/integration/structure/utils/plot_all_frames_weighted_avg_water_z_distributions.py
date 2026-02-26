"""Plot all-frames weighted-average water z-distribution profiles."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import iread

REPO_ROOT = Path(__file__).resolve().parents[4]

from src.structure.utils.WaterParser import _compute_water_mass_density_z_distribution as compute_water_mass_density_z_distribution
from src.structure.utils.WaterParser import _compute_water_orientation_weighted_density_z_distribution as compute_water_orientation_weighted_density_z_distribution
from src.structure.utils.WaterParser import detect_water_molecule_indices
from src.structure.utils.WaterParser import get_water_oxygen_indices_array
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


def main() -> None:
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    a_A, b_A, c_A = _parse_abc_from_md_inp(md_inp_path)
    z_edges_A, z_centers_A = _z_edges_and_centers_A(c_A, DEFAULT_Z_BIN_WIDTH_A)
    bin_widths_A = np.diff(z_edges_A)
    area_xy_A2 = a_A * b_A
    bin_volumes_A3 = area_xy_A2 * bin_widths_A
    bin_volumes_cm3 = bin_volumes_A3 * ANGSTROM3_TO_CM3

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

        # Volume-time weighted accumulation:
        # rho = mass / volume, orient = cos_sum / volume
        mass_sum_g += rho * bin_volumes_cm3
        cos_sum += orient * bin_volumes_A3
        volume_time_sum_cm3 += bin_volumes_cm3

        n_frames += 1
        n_waters_total += int(water_idx.shape[0])

    if n_frames == 0:
        raise RuntimeError("No frames were read from trajectory")

    rho_weighted_avg = mass_sum_g / volume_time_sum_cm3
    orient_weighted_avg = cos_sum / (volume_time_sum_cm3 / ANGSTROM3_TO_CM3)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].plot(z_centers_A, rho_weighted_avg, color="tab:blue", lw=1.5)
    axes[0].set_ylabel("rho_water weighted avg (g/cm^3)")
    axes[0].set_title("All-Frames Weighted-Average Water Mass Density z Distribution")
    axes[0].grid(True, ls="--", alpha=0.3)

    axes[1].plot(z_centers_A, orient_weighted_avg, color="tab:orange", lw=1.5)
    axes[1].set_xlabel("z (Angstrom)")
    axes[1].set_ylabel("orientation weighted avg (1/Angstrom^3)")
    axes[1].set_title("All-Frames Weighted-Average Water Orientation z Distribution")
    axes[1].grid(True, ls="--", alpha=0.3)

    fig.tight_layout()
    png_path = out_dir / "all_frames_weighted_avg_water_z_distributions.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    csv_path = out_dir / "all_frames_weighted_avg_water_z_distributions.csv"
    data = np.column_stack([z_centers_A, rho_weighted_avg, orient_weighted_avg])
    np.savetxt(csv_path, data, delimiter=",", header="z_A,rho_weighted_avg_g_cm3,orientation_weighted_avg_1_A3", comments="")

    print(f"frames processed: {n_frames}")
    print(f"average water molecules per frame: {n_waters_total / n_frames:.3f}")
    print(f"saved figure: {png_path}")
    print(f"saved data: {csv_path}")


if __name__ == "__main__":
    main()
