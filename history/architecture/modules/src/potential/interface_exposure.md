# src.potential — Interface Exposure

## Public API

| Function                           | Module               | Description                                      |
|------------------------------------|----------------------|--------------------------------------------------|
| `center_slab_potential_analysis`   | CenterPotential.py   | Multi-frame center slab potential → CSV + PNG     |
| `fermi_energy_analysis`            | CenterPotential.py   | Fermi energy time series → CSV + PNG              |
| `electrode_potential_analysis`     | CenterPotential.py   | Full U vs SHE pipeline → CSV + PNG                |
| `thickness_sensitivity_analysis`   | CenterPotential.py   | Thickness sweep → U vs SHE mean + spatial std φ(z) |
| `phi_z_planeavg_analysis`          | PhiZProfile.py       | Full-frame φ(z) profiles → CSV + PNG              |

## Config Constants

All default output filenames are defined in `config.py`:

| Constant                                  | Value                                 |
|-------------------------------------------|---------------------------------------|
| `DEFAULT_CENTER_POTENTIAL_CSV_NAME`       | `center_potential.csv`                |
| `DEFAULT_CENTER_POTENTIAL_PNG_NAME`       | `center_potential.png`                |
| `DEFAULT_FERMI_ENERGY_CSV_NAME`           | `fermi_energy.csv`                    |
| `DEFAULT_FERMI_ENERGY_PNG_NAME`           | `fermi_energy.png`                    |
| `DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME`    | `electrode_potential_U_vs_SHE.csv`    |
| `DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME`    | `electrode_potential_U_vs_SHE.png`    |
| `DEFAULT_SLAB_CENTER_CSV_NAME`            | `slab_center_and_interfaces.csv`      |
| `DEFAULT_PHI_Z_STATS_CSV_NAME`            | `phi_z_planeavg_stats.csv`            |
| `DEFAULT_PHI_Z_PNG_NAME`                  | `phi_z_planeavg_all_frames.png`       |
| `DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME`  | `thickness_sensitivity.csv`           |
| `DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME`  | `thickness_sensitivity.png`           |
