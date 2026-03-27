# charge.Bader — Implementation Guidelines

## Role

Core Bader charge analysis sub-package. Handles surface charge density (two methods), per-atom charge tracking, and counterion charge detection across MD trajectories.

## File Organization

- `BaderData.py` — `BaderTrajectoryData` frozen dataclass + `load_bader_trajectory()` (remap via IndexMap to XYZ order)
- `SurfaceCharge.py` — `compute_frame_surface_charge()` (counterion/layer methods), `trajectory_surface_charge()`, `surface_charge_analysis()` (end-to-end CSV+PNG with optional calibration)
- `AtomCharges.py` — `frame_indexed_atom_charges()`, `trajectory_indexed_atom_charges()`, `tracked_atom_charge_analysis()`, `counterion_charge_analysis()`
- `_frame_utils.py` — `_sorted_frame_dirs()` (numeric sort by `_t(\d+)` regex), step/time extraction

## Two Charge Methods

- **counterion**: non-water, non-metal species only; sigma = -sum(q_counterion) / area (charge neutrality)
- **layer**: sum net charges of N outermost metal layers per interface / area

## Key Design

- Output ordering: `[sigma_aligned, sigma_opposed]` — by stable `interface_label`, not `center_frac`
- MIC-based directional assignment for counterion method
- `_sorted_frame_dirs()` sorts by numeric `_t(\d+)` value (not lexicographic)
- Calibration integration: `surface_charge_analysis()` auto-loads `calibration.json` for sigma->phi extrapolation

## Dependencies

- `utils.BaderParser` (load_bader_atoms)
- `utils.StructureParser` (LayerParser, WaterParser)
- `utils.config` (AREA_VECTOR_INDICES, AXIS_MAP)
- `scripts.utils.IndexMapper` (remap_array, read_index_map_from_poscar)
- `electrochemical.calibration` (optional, for phi extrapolation)
