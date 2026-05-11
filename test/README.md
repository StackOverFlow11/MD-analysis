# Test Directory Layout

Tests for the `md_analysis` package (`src/md_analysis/`).

```text
test/
  README.md
  conftest.py                    # shared fixtures (e.g. parse_abc_from_md_inp)
  unit/
    test_config.py               # user config persistence (load/save/get/set)
    test_logging_setup.py        # logging NullHandler + CLI handler setup
    calibration/
      test_data.py               # CalibrationData dataclass, CSV/JSON I/O
      test_mapper.py             # ChargePotentialMapper implementations (Linear/Poly/Spline)
      test_reference.py          # SHE/RHE/PZC reference conversion
      test_workflow.py           # calibrate() + predict_potential() workflow
    charge/
      test_charge_analysis.py    # BaderAnalysis: surface charge (both methods), PBC invariance
    cli/
      test_handle_cmd_error.py   # MenuCommand error handling
      test_settings_defaults.py  # SetAnalysisDefaultCmd / ResetDefaultsCmd
    enhanced_sampling/
      constrained_ti/
        test_acf_core.py         # ACF autocorrelation core
        test_block_average.py    # Flyvbjerg-Petersen block averaging
        test_correction.py       # constant-potential free energy correction
        test_integration.py      # trapezoidal integration weights, SEM targets
        test_io.py               # discover_ti_points, load_ti_series
        test_workflow.py         # analyze_standalone, analyze_ti orchestration
    potential/
      test_single_frame_electrode.py  # single-frame electrode potential computation
    scripts/
      test_bader_gen.py          # BaderGen workdir generation
      test_ti_gen.py             # TIGen workdir generation
      utils/
        test_index_mapper.py     # CP2K XYZ <-> VASP POSCAR index mapping
    utils/
      test_bader_parser.py       # BaderParser: ACF.dat + POTCAR parsing
      test_cell_parser.py        # CellParser: CP2K cell parameter parsing
      test_cluster_utils.py      # ClusterUtils: periodic clustering, gap detection
      test_cube_parser.py        # CubeParser: cube file I/O
      test_layer_parser.py       # LayerParser: interface detection, center_frac, interface_label
      test_slowgrowth_parser.py  # ColvarParser: COLVAR restart + LagrangeMultLog
      test_water_parser.py       # WaterParser: water topology, O-H matching
  integration/
    test_main.py                 # main.py programmatic entry point smoke tests
    charge/
      test_charge_trajectory.py       # trajectory surface charge, CSV output
    enhanced_sampling/
      test_constrained_ti_regression.py  # constrained TI end-to-end regression
      test_slowgrowth_plot.py            # slowgrowth analysis plot output
    potential/
      test_center_potential.py         # center slab potential analysis
      test_phi_z_profile.py            # phi(z) plane-averaged profile
    utils/
      test_water_layer_pipeline.py     # layer + water detection end-to-end
      plot_all_frames_auto_c_window_theta_pdf.py          # (visual sanity check)
      plot_all_frames_weighted_avg_water_z_distributions.py  # (visual sanity check)
      plot_last_frame_water_z_distributions.py             # (visual sanity check)
    water/
      test_adsorbed_water_theta_distribution.py
      test_water_density_orientation_plots_angstrom_only.py
      test_water_density_plot_from_potential.py
      test_water_orientation_plot_from_potential.py
      test_water_three_panel_plot.py
```

## Running tests

```bash
pip install .          # install package first (not editable)
pytest test/           # run all tests
pytest test/unit/      # unit tests only
pytest test/unit/charge/test_charge_analysis.py -v  # single file
```

## Conventions

- `unit/`: small, deterministic tests for single functions.
- `integration/`: frame-level workflows that read real data from `data_example/`.
- Test data lives in `data_example/` (potential cube files, Bader POSCAR/ACF/POTCAR, etc.).
- Integration tests under `integration/utils/` prefixed with `plot_` are visual sanity-check scripts, not automated assertions.
