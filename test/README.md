# Test Directory Layout

Tests for the `md_analysis` package (`src/md_analysis/`).

```text
test/
  README.md
  conftest.py              # shared fixtures (e.g. parse_abc_from_md_inp)
  unit/
    utils/
      test_layer_parser.py       # LayerParser: interface detection, center_frac, interface_label
      test_water_parser.py       # WaterParser: water topology, O-H matching
      test_cluster_utils.py      # ClusterUtils: periodic clustering, gap detection
      test_bader_parser.py       # BaderParser: ACF.dat + POTCAR parsing
    charge/
      test_charge_analysis.py    # BaderAnalysis: surface charge (both methods), PBC invariance
    potential/
      test_single_frame_electrode.py  # single-frame electrode potential computation
  integration/
    utils/
      test_water_layer_pipeline.py    # layer + water detection end-to-end
    charge/
      test_charge_trajectory.py       # trajectory surface charge, CSV output
    potential/
      test_center_potential.py        # center slab potential analysis
      test_phi_z_profile.py           # φ(z) plane-averaged profile
    water/
      test_water_three_panel_plot.py  # three-panel water analysis plot
      test_adsorbed_water_theta_distribution.py
      test_water_density_plot_from_potential.py
      test_water_orientation_plot_from_potential.py
      test_water_density_orientation_plots_angstrom_only.py
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
