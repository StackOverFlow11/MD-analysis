# CLAUDE.md

## Project Overview

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations. Standard src-layout, Python 3.10+, distributed as `md-analysis` (pip), importable as `md_analysis`.

## Quick Start

```bash
pip install numpy matplotlib ase pytest tqdm
pip install .          # install package (required before running tests)
pytest test/           # run all tests
```

Entry point: `md-analysis` console script → `md_analysis.cli:main` (VASPKIT-style interactive menu).

## Project Map

### `src/md_analysis/` — Package Root

| File | Purpose |
|---|---|
| `__init__.py` | Re-exports sub-packages (`utils`, `water`, `potential`, `charge`), `MDAnalysisError`; `__version__ = "0.1.0"` |
| `exceptions.py` | `MDAnalysisError(Exception)` — base class for all domain-specific errors; all 9 custom exceptions inherit from it |
| `config.py` | Persistent user configuration (`~/.config/md_analysis/config.json`): `load_config`, `save_config`, `get_config`, `set_config`, `delete_config`, `ConfigError`, `KEY_VASP_SCRIPT_PATH`, `KEY_LAYER_TOL_A`, `KEY_Z_BIN_WIDTH_A`, `KEY_THETA_BIN_DEG`, `KEY_WATER_OH_CUTOFF_A`, `CONFIGURABLE_DEFAULTS` |
| `main.py` | Programmatic entry points: `run_water_analysis()`, `run_potential_analysis()`, `run_charge_analysis()`, `run_all()` |

### `src/md_analysis/cli/` — Interactive CLI (no argparse)

| File | Purpose | Key Functions |
|---|---|---|
| `__init__.py` | Top menu (1=Water, 2=Potential, 3=Charge, 4=Scripts, 9=Settings, 0=Exit) | `main()` |
| `_prompt.py` | Reusable input helpers + error decorator | `_get_effective_default`, `_prompt_str`, `_prompt_int`, `_prompt_float`, `_prompt_choice`, `_prompt_bool`, `_parse_metal_elements`, `_prompt_global_params`, `_handle_cmd_error` (decorator: catches exceptions in `_cmd_*` handlers, prints user-friendly message, returns 1) |
| `_water.py` | Water sub-menu (101-105) | `water_menu()`, dispatches to `water` module or `main.run_water_analysis` |
| `_potential.py` | Potential sub-menu (201-206) | `potential_menu()`, dispatches to `potential` module or `main.run_potential_analysis` |
| `_charge.py` | Charge sub-menu (301-303) | `charge_menu()`, `_collect_params()`, `_run_charge()`, `_print_ensemble_summary()` |
| `_scripts.py` | Scripts/Tools sub-menu (401-402) | `scripts_menu()`, 401: single-frame Bader workdir, 402: batch Bader workdirs |
| `_settings.py` | Settings sub-menu (901-907) | `settings_menu()`, 901: set VASP script path, 902: show config, 903-906: set analysis defaults, 907: reset all defaults |

**Menu codes:**
- 101: mass density, 102: orientation-weighted density, 103: adsorbed orientation, 104: theta distribution, 105: full three-panel
- 201: center slab potential, 202: Fermi energy, 203: electrode potential (U vs SHE), 204: phi(z) profile, 205: thickness sensitivity, 206: full potential
- 301: surface charge (counterion), 302: surface charge (layer), 303: full charge (method prompted)
- 401: generate Bader work directory (single frame), 402: batch generate Bader work directories
- 901: set VASP submission script path, 902: show current configuration
- 903: set layer clustering tolerance, 904: set z-axis bin width, 905: set theta bin width, 906: set water O-H cutoff, 907: reset all analysis defaults

### `src/md_analysis/utils/` — Shared Low-Level Utilities

| File | Key Exports | Purpose |
|---|---|---|
| `config.py` | `TRANSITION_METAL_SYMBOLS`, `DEFAULT_METAL_SYMBOLS`, `DEFAULT_Z_BIN_WIDTH_A`, `DEFAULT_THETA_BIN_DEG`, `DEFAULT_WATER_OH_CUTOFF_A`, `WATER_MOLAR_MASS_G_PER_MOL`, `HA_TO_EV`, `BOHR_TO_ANG`, `DP_A_H3O_W_EV`, `MU_HPLUS_G0_EV`, `DELTA_E_ZP_EV`, `AXIS_MAP`, `AREA_VECTOR_INDICES`, `INTERFACE_NORMAL_ALIGNED`, `INTERFACE_NORMAL_OPPOSED`, `CHARGE_METHOD_COUNTERION`, `CHARGE_METHOD_LAYER` | Physical constants + defaults |
| `_io_helpers.py` | `_cumulative_average(values)`, `_write_csv(path, rows, fieldnames)` | Private shared I/O helpers (used by potential + charge workflows) |
| `CubeParser.py` | `CubeHeader` (dataclass), `read_cube_header_and_values(path)`, `slab_average_potential_ev(header, values, thickness, *, z_center)`, `plane_avg_phi_z_ev()`, `z_coords_ang()`, `discover_cube_files(cube_pattern, *, workdir, frame_start, frame_end, frame_step)`, `extract_step_from_cube_filename()` | Gaussian cube file I/O + potential utilities |
| `BaderParser.py` | `load_bader_atoms(structure, acf, potcar)` → Atoms with `bader_charge`/`bader_net_charge` arrays, `BaderParseError`, `_read_acf()`, `_read_potcar_zval()` | VASP Bader charge parsing |
| `StructureParser/ClusterUtils.py` | `cluster_1d_periodic(values, period, tol)`, `find_largest_gap_periodic(centers, period)`, `gap_midpoint_periodic(low, high, period)` | 1D periodic clustering and gap detection |
| `StructureParser/LayerParser.py` | `Layer` (dataclass: `atom_indices`, `center_frac`, `is_interface`, `interface_label`, `normal_unit`), `SurfaceDetectionResult` (`.interface_normal_aligned()`, `.interface_normal_opposed()`), `SurfaceGeometryError`, `detect_interface_layers(atoms, *, metal_symbols, normal, layer_tol_A)`, `format_detection_summary()`, `_circular_mean_fractional()`, `_mic_delta_fractional()` | Metal layer detection + interface labeling |
| `StructureParser/WaterParser.py` | `detect_water_molecule_indices(atoms)` → `(n_water, 3)`, `get_water_oxygen_indices_array()`, `WaterTopologyError`, `_compute_bisector_cos_theta_vec()`, `_oxygen_to_hydrogen_map()`, `_compute_water_mass_density_z_distribution()`, `_compute_water_orientation_weighted_density_z_distribution()`, `_compute_water_orientation_theta_pdf_in_c_fraction_window()` | Water topology + single-frame z-profiles |
| `RestartParser/CellParser.py` | `parse_abc_from_restart(restart_path)` → `(a, b, c)`, `parse_abc_from_md_inp(md_inp_path)` → `(a, b, c)`, `CellParseError` | CP2K cell parameter parsing (`.restart` + `md.inp`) |
| `RestartParser/ColvarParser.py` | `ColvarParseError`, `ConstraintInfo`, `ColvarInfo`, `ColvarRestart`, `LagrangeMultLog`, `parse_colvar_restart(restart_path)` → `ColvarRestart`, `parse_lagrange_mult_log(log_path)` → `LagrangeMultLog`, `compute_target_series(restart, n_steps)` → `np.ndarray` | CP2K COLVAR restart + LagrangeMultLog parsing |

### `src/md_analysis/water/` — Water Analysis Workflows

| File | Key Exports | Purpose |
|---|---|---|
| `config.py` | `DEFAULT_START_INTERFACE = "normal_aligned"`, output CSV/PNG name constants | Water analysis defaults |
| `__init__.py` | Re-exports all public water API + utils types | Package interface |
| `Water.py` | `plot_water_three_panel_analysis(xyz_path, md_inp_path, *, output_dir, ...)` → PNG path | Integrated 3-panel figure (density + orientation + theta PDF) |
| `WaterAnalysis/__init__.py` | Re-exports from sub-modules; internal: `StartInterface`, `_compute_density_orientation_ensemble` | Sub-package interface |
| `WaterAnalysis/_common.py` | `_detect_interface_fractions()`, `_iter_trajectory()`, `_single_frame_density_and_orientation()`, `_compute_density_orientation_ensemble()` (re-imports `_parse_abc_from_md_inp` from `utils.RestartParser.CellParser`) | Shared private trajectory/frame helpers |
| `WaterAnalysis/WaterDensity.py` | `water_mass_density_z_distribution_analysis(xyz, md_inp, *, output_dir, ...)` → CSV path | Ensemble-averaged mass density profile |
| `WaterAnalysis/WaterOrientation.py` | `water_orientation_weighted_density_z_distribution_analysis(xyz, md_inp, *, output_dir, ...)` → CSV path | Ensemble-averaged orientation-weighted density |
| `WaterAnalysis/AdWaterOrientation.py` | `ad_water_orientation_analysis(xyz, md_inp, *, output_dir, ...)` → (CSV, TXT), `compute_adsorbed_water_theta_distribution(...)` → (centers, pdf, CSV), `detect_adsorbed_layer_range_from_density_profile(distance, rho)` | Adsorbed-layer orientation analysis |

### `src/md_analysis/potential/` — Potential Analysis Workflows

| File | Key Exports | Purpose |
|---|---|---|
| `config.py` | `DEFAULT_THICKNESS_ANG = 7.5`, output CSV/PNG name constants | Potential analysis defaults |
| `__init__.py` | Re-exports all public potential API + config constants | Package interface |
| `CenterPotential.py` | `center_slab_potential_analysis(cube_pattern, *, output_dir, thickness_ang, center_mode, ...)` → CSV path, `fermi_energy_analysis(md_out, *, output_dir, fermi_unit, ...)` → CSV path, `electrode_potential_analysis(cube_pattern, md_out, *, ...)` → CSV path, `thickness_sensitivity_analysis(cube_pattern, md_out, *, thickness_end, ...)` → CSV path | Multi-frame potential + electrode potential (U vs SHE) |
| `PhiZProfile.py` | `phi_z_planeavg_analysis(cube_pattern, *, output_dir, max_curves, ...)` → PNG path | phi(z) plane-averaged overlay plot |

**cSHE formula:** `U = -E_Fermi + phi_center + DP_A(H3O+/w) - mu(H+,g0) - DELTA_E_ZP`

### `src/md_analysis/charge/` — Bader Charge Analysis

| File | Key Exports | Purpose |
|---|---|---|
| `config.py` | `E_PER_A2_TO_UC_PER_CM2 = 1602.176634`, `DEFAULT_DIR_PATTERN = "bader_t*_i*"`, `DEFAULT_STRUCTURE_FILENAME = "POSCAR"`, `DEFAULT_ACF_FILENAME = "ACF.dat"`, `DEFAULT_POTCAR_FILENAME = "POTCAR"`, output CSV/PNG names | Charge analysis defaults |
| `__init__.py` | Re-exports 5 functions + 3 constants | Package interface |
| `BaderAnalysis.py` | `compute_frame_surface_charge(atoms, *, metal_symbols, normal, method)` → mutated Atoms, `frame_indexed_atom_charges(atoms, indices)` → `(N,2)`, `trajectory_indexed_atom_charges(root_dir, matrix, ...)` → `(t,N,2)`, `trajectory_surface_charge(root_dir, *, method, ...)` → `(t,2)` uC/cm2, `surface_charge_analysis(root_dir, *, method, output_dir, ...)` → CSV path | Surface charge density (counterion/layer methods) |

**Charge methods:** `"counterion"` (exclude water+metal, use solute species) / `"layer"` (net charge of interface metal layers).
**atoms.info keys set:** `surface_charge_density_e_A2`, `surface_charge_density_uC_cm2`, `n_charged_atoms_per_surface`, `charge_per_surface_e` — all `[aligned, opposed]`.

### `src/md_analysis/scripts/` — Automation Scripts & Utilities

| File | Key Exports | Purpose |
|---|---|---|
| `__init__.py` | Re-exports `BaderGenError`, `generate_bader_workdir`, `batch_generate_bader_workdirs` | Package interface |
| `BaderGen.py` | `BaderGenError`, `DEFAULT_WORKDIR_NAME`, `generate_bader_workdir(atoms, output_dir, *, ...)` → `Path`, `batch_generate_bader_workdirs(xyz_path, cell_abc, output_dir, *, frame_start, frame_end, frame_step, ...)` → `list[Path]` | Generate VASP Bader work directories (single + batch) |
| `template/__init__.py` | (empty) | Package marker for `importlib.resources` |
| `template/INCAR` | — | VASP INCAR template for Bader single-point |
| `template/KPOINTS` | — | VASP KPOINTS template (Gamma-only) |
| `utils/__init__.py` | Re-exports all `IndexMapper` public symbols | Sub-package interface |
| `utils/IndexMapper.py` | `IndexMap` (dataclass), `IndexMapError`, `IndexMapParseError`, `compute_index_map(atoms_xyz, *, frame, source, element_order)` → `IndexMap`, `write_poscar_with_map(atoms_xyz, output_path, index_map, *, direct)` → `Path`, `read_index_map_from_poscar(poscar_path)` → `IndexMap`, `encode_comment_line(index_map)` → `str`, `decode_comment_line(comment)` → `IndexMap`, `remap_array(data, index_map, direction)` → `np.ndarray` | CP2K XYZ ↔ VASP POSCAR bijective index mapping |

**POSCAR comment format:** `md_analysis::v1 frame=<int> source=<urlencoded> n=<int> order=<csv> p2x=<base64>`
**Not re-exported** from top-level `md_analysis.__init__`.

### `test/` — Test Suite

| Path | Coverage |
|---|---|
| `test/unit/utils/test_cluster_utils.py` | `cluster_1d_periodic`, `find_largest_gap_periodic`, `gap_midpoint_periodic` |
| `test/unit/utils/test_layer_parser.py` | `detect_interface_layers`, `Layer`, `SurfaceDetectionResult` |
| `test/unit/utils/test_water_parser.py` | `detect_water_molecule_indices`, density/orientation helpers |
| `test/unit/utils/test_bader_parser.py` | `_read_acf`, `_read_potcar_zval`, `load_bader_atoms` |
| `test/unit/utils/test_cube_parser.py` | `read_cube_header_and_values`, `slab_average_potential_ev`, `plane_avg_phi_z_ev`, `discover_cube_files` |
| `test/unit/utils/test_cell_parser.py` | `parse_abc_from_restart` + `parse_abc_from_md_inp`: valid/missing/error cases |
| `test/unit/utils/test_slowgrowth_parser.py` | `ColvarInfo`, `parse_colvar_restart` (4 scenarios), `parse_lagrange_mult_log` (single/multi), `compute_target_series`, multi-COLLECTIVE blocks, edge cases |
| `test/unit/test_config.py` | `load_config`, `save_config`, `get_config`, `set_config`, `delete_config`, `ConfigError`, `CONFIGURABLE_DEFAULTS` |
| `test/unit/test_logging_setup.py` | NullHandler on library logger, StreamHandler setup in CLI |
| `test/unit/cli/test_handle_cmd_error.py` | `_handle_cmd_error` decorator: known/unexpected exceptions, normal return |
| `test/unit/cli/test_settings_defaults.py` | `_get_effective_default`, settings handlers 903-907, 902 enhanced display |
| `test/unit/charge/test_charge_analysis.py` | All 5 charge public functions + internal `_extract_t_value`, `_sorted_frame_dirs` |
| `test/unit/potential/test_single_frame_electrode.py` | Single-frame potential pipeline |
| `test/unit/scripts/test_bader_gen.py` | `generate_bader_workdir` + `batch_generate_bader_workdirs`: directory structure, frame slicing, metadata |
| `test/unit/scripts/utils/test_index_mapper.py` | `compute_index_map`, `encode/decode_comment_line`, `write/read_poscar_with_map`, `remap_array` |
| `test/integration/test_main.py` | Programmatic entry points: `run_water_analysis`, `run_potential_analysis`, `run_charge_analysis` |
| `test/integration/utils/test_water_layer_pipeline.py` | Water + layer detection combined pipeline |
| `test/integration/charge/test_charge_trajectory.py` | `trajectory_*` functions, `surface_charge_analysis` end-to-end |
| `test/integration/potential/test_center_potential.py` | `center_slab_potential_analysis` with real cube files |
| `test/integration/potential/test_phi_z_profile.py` | `phi_z_planeavg_analysis` with real cube files |
| `test/integration/water/test_water_*.py` | Water density/orientation/theta plots, three-panel integration |
| `test/conftest.py` | Exports `parse_abc_from_md_inp` helper |

**Test data:** `data_example/bader/bader_work_dir/` (POSCAR, ACF.dat, POTCAR — 274 atoms: 62 Cu + 2 Ag + 70 O + 140 H), `data_example/potential/` (cube files, md.out, xyz), `data_example/sg/` (COLVAR restart + LagrangeMultLog, 4 scenarios).

### `context4agent/` — Documentation Mirror (Single Source of Truth)

Mirrors `src/md_analysis/` structure. Each sub-package has `interface_exposure.md` + `implementation_guidelines.md`. Also contains `data_contract.md` (CSV specs), `glossary_units.md`, `requirements/short_term.md`.

## Logging Conventions

- **Library code** uses `logging.getLogger(__name__)` per module; `NullHandler` is set on the root `md_analysis` logger in `__init__.py` (PEP 282 best practice)
- **CLI** configures `StreamHandler(stderr)` at `INFO` level in `main()` — only when no application handler exists
- **`verbose` parameter** controls tqdm progress bars only; logging messages are always emitted (visibility depends on handler configuration)
- **Log placement**: only at workflow boundaries and before long operations; never in tight loops or low-level math functions
- **`_handle_cmd_error`**: logs unexpected exceptions with `logger.error(..., exc_info=True)` for full traceback; prints concise message to stdout

## Key Conventions

- **Interface labels:** `"normal_aligned"` (+axis facing) / `"normal_opposed"` (-axis facing)
- **Layer ordering:** `[normal_aligned, slab_interior..., normal_opposed]`
- **Surface normal:** only `"a"/"b"/"c"` cell axes supported (no custom vectors)
- **Frame dirs:** `bader_t*_i*` pattern, numerically sorted by `_t(\d+)` regex
- **Output dir structure:** `<outdir>/water/`, `<outdir>/potential/<sub>/`, `<outdir>/charge/<method>/`
- **Relative imports** inside `md_analysis/`; **absolute imports** in tests

## context4agent — Sync Rules (MANDATORY)

`context4agent/` is the authoritative reference for architecture, API contracts, units, and project status. **Every code change MUST synchronize the corresponding `context4agent/` docs.**

| What changed | Update where |
|---|---|
| Module / sub-package added/removed/renamed | `architecture/README.md` + `architecture/modules/src/<module>/` mirror docs |
| Public API (function, signature, export) | corresponding `interface_exposure.md` |
| Implementation pattern / dependency | corresponding `implementation_guidelines.md` |
| Output format (CSV columns, filenames, shapes) | `architecture/modules/data_contract.md` |
| Units / terminology | `architecture/modules/glossary_units.md` |
| Top-level `__all__` | `architecture/modules/src/interface_exposure.md` + `implementation_guidelines.md` |
| Project capabilities / status | `requirements/short_term.md` |

## Workflow Rules

- After completing any module's code changes, always run /commit.
- 始终使用中文回复用户。
