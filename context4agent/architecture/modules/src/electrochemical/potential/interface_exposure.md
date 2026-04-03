# md_analysis.potential — Interface Exposure

## Public API

| Function                           | Module               | Description                                      |
|------------------------------------|----------------------|--------------------------------------------------|
| `center_slab_potential_analysis`   | CenterPotential.py   | Multi-frame center slab potential → CSV + PNG     |
| `fermi_energy_analysis`            | CenterPotential.py   | Fermi energy time series → CSV + PNG              |
| `electrode_potential_analysis`     | CenterPotential.py   | Full U vs SHE pipeline → CSV + PNG                |
| `thickness_sensitivity_analysis`   | CenterPotential.py   | Thickness sweep → U vs SHE mean + spatial std φ(z) |
| `phi_z_planeavg_analysis`          | PhiZProfile.py       | Full-frame φ(z) profiles → CSV + PNG              |

## Frame Discovery API (new)

| Symbol                           | Module               | Description                                      |
|------------------------------------|----------------------|--------------------------------------------------|
| `PotentialFrame`                 | _frame_source.py     | Frozen dataclass: step, time_fs, cube_path, header, values, fermi_raw, atoms |
| `discover_continuous_frames`     | _frame_source.py     | Discover frames from single directory (mode A)   |
| `discover_distributed_frames`    | _frame_source.py     | Discover frames from `potential_t*_i*/` subdirs (mode B) |

## Input Mode Parameters (all 5 analysis functions)

All public analysis functions now accept these optional keyword arguments (default = continuous mode behavior):

| Parameter          | Type           | Default                              |
|--------------------|----------------|--------------------------------------|
| `input_mode`       | `str`          | `"continuous"`                       |
| `sp_root_dir`      | `Path \| str`  | `None`                               |
| `sp_dir_pattern`   | `str`          | `"potential_t*_i*"`                  |
| `sp_cube_filename` | `str`          | `"sp_potential-v_hartree-1_0.cube"`  |
| `sp_out_filename`  | `str`          | `"sp.out"`                           |

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
| `DEFAULT_THICKNESS_ANG`                  | `7.5`                                 |
| `DEFAULT_SP_DIR_PATTERN`                 | `"potential_t*_i*"`                   |
| `DEFAULT_SP_CUBE_FILENAME`               | `"sp_potential-v_hartree-1_0.cube"`   |
| `DEFAULT_SP_OUT_FILENAME`                | `"sp.out"`                            |
