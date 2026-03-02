# src.potential — Interface Exposure

## Public API

| Function                           | Module               | Description                                     |
|------------------------------------|----------------------|-------------------------------------------------|
| `center_slab_potential_analysis`   | CenterPotential.py   | Multi-frame center slab potential → CSV + PNG    |
| `fermi_energy_analysis`            | CenterPotential.py   | Fermi energy time series → CSV + PNG             |
| `electrode_potential_analysis`     | CenterPotential.py   | Full U vs SHE pipeline → CSV + PNG               |
| `phi_z_planeavg_analysis`          | PhiZProfile.py       | Full-frame φ(z) profiles → CSV + PNG             |

## Config Constants

All default output filenames are defined in `config.py`.
