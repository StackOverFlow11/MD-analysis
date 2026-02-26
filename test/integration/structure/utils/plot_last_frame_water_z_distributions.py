"""Plot water z-distribution profiles for the last frame of md-pos-1.xyz."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read

REPO_ROOT = Path(__file__).resolve().parents[4]

from src.structure.utils.WaterParser import _compute_water_mass_density_z_distribution as compute_water_mass_density_z_distribution
from src.structure.utils.WaterParser import _compute_water_orientation_weighted_density_z_distribution as compute_water_orientation_weighted_density_z_distribution
from src.structure.utils.WaterParser import detect_water_molecule_indices
from src.structure.utils.WaterParser import get_water_oxygen_indices_array
from src.structure.utils.config import DEFAULT_Z_BIN_WIDTH_A


def _parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def _z_bin_centers_A(lz_A: float, dz_A: float, nbins: int) -> np.ndarray:
    if nbins == 1:
        edges = np.array([0.0, lz_A], dtype=float)
    else:
        interior = np.arange(1, nbins, dtype=float) * dz_A
        edges = np.concatenate(([0.0], interior, [lz_A]))
    return 0.5 * (edges[:-1] + edges[1:])


def main() -> None:
    repo_root = REPO_ROOT
    xyz_path = repo_root / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = repo_root / "data_example" / "potential" / "md.inp"

    atoms = read(str(xyz_path), index=-1)
    a, b, c = _parse_abc_from_md_inp(md_inp_path)
    atoms.set_cell([a, b, c])
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

    z_centers_A = _z_bin_centers_A(lz_A=c, dz_A=DEFAULT_Z_BIN_WIDTH_A, nbins=rho.size)

    out_dir = repo_root / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].plot(z_centers_A, rho, color="tab:blue", lw=1.5)
    axes[0].set_ylabel("rho_water (g/cm^3)")
    axes[0].set_title("Last Frame Water Mass Density z Distribution")
    axes[0].grid(True, ls="--", alpha=0.3)

    axes[1].plot(z_centers_A, orient, color="tab:orange", lw=1.5)
    axes[1].set_xlabel("z (Angstrom)")
    axes[1].set_ylabel("orientation weighted density (1/Angstrom^3)")
    axes[1].set_title("Last Frame Water Orientation Weighted z Distribution")
    axes[1].grid(True, ls="--", alpha=0.3)

    fig.tight_layout()
    png_path = out_dir / "last_frame_water_z_distributions.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    csv_path = out_dir / "last_frame_water_z_distributions.csv"
    data = np.column_stack([z_centers_A, rho, orient])
    np.savetxt(csv_path, data, delimiter=",", header="z_A,rho_g_cm3,orientation_weighted_density_1_A3", comments="")

    print(f"water molecules detected: {water_idx.shape[0]}")
    print(f"saved figure: {png_path}")
    print(f"saved data: {csv_path}")


if __name__ == "__main__":
    main()
