# md_analysis.charge — Interface Exposure

## Public API

| Symbol                           | Module              | Description                                           |
|----------------------------------|---------------------|-------------------------------------------------------|
| `AtomSelector`                   | ChargeAnalysis.py   | ABC for atom selection (`select(atoms) -> ndarray`)   |
| `ElementSelector`                | ChargeAnalysis.py   | Select atoms by element symbol(s)                     |
| `IndexSelector`                  | ChargeAnalysis.py   | Select atoms by fixed 0-based indices                 |
| `TrajectoryChargeResult`         | ChargeAnalysis.py   | Frozen dataclass for multi-frame results              |
| `compute_frame_surface_charge`   | ChargeAnalysis.py   | Single-frame surface charge density → `atoms.info`    |
| `trajectory_charge_analysis`     | ChargeAnalysis.py   | Multi-frame trajectory analysis → result + CSV        |
| `E_PER_A2_TO_UC_PER_CM2`        | config.py           | Unit conversion: 1 e/Å² = 1.602176634×10³ μC/cm²     |

## Config Constants

All default filenames and constants are defined in `config.py`:

| Constant                                  | Value                            |
|-------------------------------------------|----------------------------------|
| `E_PER_A2_TO_UC_PER_CM2`                 | `1.602176634e3`                  |
| `DEFAULT_DIR_PATTERN`                     | `calc_t*_i*`                     |
| `DEFAULT_STRUCTURE_FILENAME`              | `POSCAR`                         |
| `DEFAULT_ACF_FILENAME`                    | `ACF.dat`                        |
| `DEFAULT_POTCAR_FILENAME`                 | `POTCAR`                         |
| `DEFAULT_SURFACE_CHARGE_CSV_NAME`         | `surface_charge_density.csv`     |
| `DEFAULT_SELECTED_ATOM_CHARGES_CSV_NAME`  | `selected_atom_charges.csv`      |

## TrajectoryChargeResult Fields

| Field                                 | Shape / Type          | Description                            |
|---------------------------------------|-----------------------|----------------------------------------|
| `frame_labels`                        | `tuple[str, ...]`     | Directory names (sorted)               |
| `surface_charge_density_uC_cm2`      | `(n_frames, 2)`      | [σ_bottom, σ_top] per frame            |
| `mean_surface_charge_density_uC_cm2` | `(2,)`               | Ensemble mean                          |
| `std_surface_charge_density_uC_cm2`  | `(2,)`               | Ensemble std                           |
| `selected_atom_net_charges`           | `(n_frames, n_sel)`  | Per-atom net charge (or empty)         |
| `mean_selected_atom_net_charges`      | `(n_sel,)`           | Mean per-atom net charge (or empty)    |
| `selected_atom_indices`               | `tuple[int, ...]`    | 0-based indices (or empty tuple)       |
