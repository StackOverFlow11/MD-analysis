# charge.Bader — Interface Exposure

## Public API

### SurfaceCharge.py

| Function | Description |
|----------|-------------|
| `compute_frame_surface_charge(atoms, *, method, normal, ...)` | Single-frame surface charge -> atoms.info |
| `trajectory_surface_charge(root_dir, *, method, normal, ...)` | Multi-frame -> (t, 2) uC/cm2 |
| `surface_charge_analysis(root_dir, *, method, ...)` | End-to-end CSV+PNG (with optional calibration phi) |

### AtomCharges.py

| Function | Description |
|----------|-------------|
| `frame_indexed_atom_charges(atoms, atom_indices)` | Single-frame -> (N, 2) |
| `trajectory_indexed_atom_charges(root_dir, atom_index_matrix, ...)` | Multi-frame -> (t, N, 2) |
| `tracked_atom_charge_analysis(root_dir, *, atom_indices_xyz, ...)` | XYZ-indexed tracking -> CSV+PNG |
| `counterion_charge_analysis(root_dir, *, normal, ...)` | Per-frame counterion detection -> CSV+PNG |

### BaderData.py

| Symbol | Description |
|--------|-------------|
| `BaderTrajectoryData` | Frozen dataclass: steps, times, atom_indices_xyz, net_charges (XYZ order) |
| `load_bader_trajectory(root_dir, ...)` | Load all frames -> BaderTrajectoryData |

## Stability

Stable — used by CLI 221-226 and constrained_ti.correction.
