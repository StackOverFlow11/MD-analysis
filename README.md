# MD Analysis

Lightweight analysis utilities for periodic metal-water interfaces from CP2K MD simulations.

## Features

**Water analysis** — density and orientation profiles along the interface-to-midpoint direction:

- Water mass-density profile (g/cm³)
- Orientation-weighted density profile (g/cm³)
- Adsorbed-layer auto-detection and orientation-angle distribution
- Integrated three-panel PNG visualization

**Potential analysis** — Hartree potential and electrode potential from cube files:

- Center-slab Hartree potential time series
- Fermi energy extraction from `md.out`
- Electrode potential U vs SHE (computational SHE)
- Plane-averaged φ(z) profile overlay across frames
- Thickness sensitivity sweep (dual-axis: U vs SHE + spatial std of φ(z))

**Charge analysis** — Bader surface charge density from VASP post-processing:

- Two surface charge methods:
  - `counterion` — only non-water, non-metal species (counterions, adsorbates) contribute to σ
  - `layer` — sum of interface-layer metal atom net charges / area
- Single-side surface charge + potential extrapolation via calibration
- Time-series CSV with cumulative average + PNG visualization
- Per-atom indexed charge extraction and counterion charge tracking across trajectory

**Charge-potential calibration** — σ→φ mapping for constant-charge simulations:

- Fit calibration from CSV or manual (φ, σ) data points (linear, polynomial, spline, differential capacitance)
- Predict electrode potential from surface charge density
- Reference scale conversion (SHE / RHE / PZC)

**Enhanced sampling** — slow-growth and constrained thermodynamic integration:

- Slow-growth TI: quick and publication-quality free energy plots
- Constrained TI: 4-step convergence diagnostics (ACF, Flyvbjerg-Petersen block average, running average drift, Geweke stationarity)
- Free-energy integration with error propagation
- Constant-potential correction (Norskov) for constant-charge TI

**Scripts & Tools** — automation for simulation workflows:

- VASP Bader work directory generation (single frame + batch from trajectory)
- Bijective CP2K XYZ ↔ VASP POSCAR index mapping (encoded in POSCAR comment line)
- CP2K constrained-MD work directory generation for TI (single target + batch)
- Persistent configuration for VASP/CP2K submission script paths

## Quick Start

### Install

```bash
pip install numpy matplotlib ase pytest tqdm scipy
pip install -e .   # editable install — recommended for development;
                   # without -e, tests may import a stale site-packages copy.
```

If you skip the editable install (or use a non-editable `pip install .`),
prefix every test command with `PYTHONPATH=src` so pytest imports the
current working tree instead of any older copy installed in your
environment.

#### Offline / HPC install

`cp2k-input-tools` is pinned as a git direct dependency
(`cp2k-input-tools @ git+https://github.com/cp2k/cp2k-input-tools.git@393e87b`),
so a plain `pip install .` needs network access to GitHub. For HPC
compute nodes without internet, pre-download all wheels on a networked
machine and install from the local directory:

```bash
# on a machine with network access (same Python/platform as the target):
pip download . -d wheels/

# copy the wheels/ directory to the HPC login node, then offline:
pip install --no-index --find-links=wheels/ .
```

The `pip download` step is the only part that requires connectivity
(it clones the git-pinned dependency and builds its wheel); the install
step afterwards runs fully offline.

### CLI (interactive menu)

Run `md-analysis` with no arguments to enter the interactive menu:

```bash
md-analysis
```

The VASPKIT-style numbered menu guides you through:

- **1) Water Analysis** — density (101), orientation (102), adsorbed-layer (103/104), three-panel (105)
- **2) Electrochemical Analysis**
  - **21) Potential** — center slab (211), Fermi (212), U vs SHE (213), phi(z) (214), thickness sweep (215), full (216)
  - **22) Charge** — counterion (221), layer (222), full (223), single-side + potential (224), tracked atoms (225), counterion tracking (226)
  - **23) Calibration** — calibrate from CSV (231), manual input (232), predict potential (233)
- **3) Enhanced Sampling**
  - **30) Slow-Growth** — quick plot (301), publication plot (302)
  - **31) Constrained TI** — single-point diagnostics (311), full TI analysis (312), constant-potential correction (313)
- **4) Scripts / Tools**
  - **41) Bader** — single frame (411), batch (412)
  - **42) TI** — single target (421), batch (422)
  - **43) SP Potential** — single (431), batch (432)
  - **44) DeePMD SP** — single (441), batch (442)
- **9) Settings**
  - **90) Show / Reset** — show config (900), reset all (909)
  - **91) Script Paths** — VASP submit (911), CP2K submit (912), SP inp template (913), DP SP inp template (914)
  - **92) Analysis Defaults** — layer tol (921), Z bin (922), theta bin (923), water O-H cutoff (924)
  - **93) Potential Output** — reference scale / pH / T / φ_PZC (931)

Each option prompts for required inputs, then offers an optional "Modify advanced parameters?" gate.

> **Chinese user docs (CLI guide / workflows / settings / pitfalls):** see [`docs/`](docs/README.md).

### Python API

```python
from md_analysis.water import plot_water_three_panel_analysis

plot_water_three_panel_analysis(
    xyz_path="data_example/potential/dense/md-pos-1.xyz",
    md_inp_path="data_example/potential/dense/md.inp",
    output_dir="output/",
)
```

Or use the programmatic entry points (canonical module:
`md_analysis.workflows`; the same names are re-exported by
`md_analysis.main` as a thin facade):

```python
from md_analysis.workflows import (
    run_water_three_panel, run_potential_full,
    run_surface_charge, run_tracked_charge, run_counterion_charge,
    run_interface_analysis,  # composite: water + potential
)

# Each run_* returns a WorkflowResult; iterate result.artifacts:
result = run_water_three_panel(
    xyz_path="data_example/potential/dense/md-pos-1.xyz",
    md_inp_path="data_example/potential/dense/md.inp",
    output_dir="output/water/",
)
for name, path in result.artifacts.items():
    print(name, path)
```

> The legacy `run_*_analysis` / `run_all` names that previously lived
> in `md_analysis.main` were removed in the entrance refactor
> (Phase 7a) without an alias — switch to the names above and read
> `WorkflowResult.artifacts` instead of the old `dict[str, Path]`.

Generate VASP Bader work directories from an MD trajectory:

```python
from md_analysis.scripts import generate_bader_workdir, batch_generate_bader_workdirs

# Single frame
from ase.io import read
atoms = read("trajectory.xyz", index=0)
atoms.set_cell([a, b, c])
atoms.set_pbc(True)
generate_bader_workdir(atoms, output_dir="bader_out/", script_path="run.sh")

# Batch across trajectory (frame 0, 10, 20, ...)
batch_generate_bader_workdirs(
    xyz_path="trajectory.xyz",
    cell_abc=(a, b, c),
    output_dir="bader_out/",
    frame_start=0, frame_end=100, frame_step=10,
    script_path="run.sh",
)
```

### Example input data

- `data_example/potential/dense/md-pos-1.xyz` — trajectory frames (continuous mode)
- `data_example/potential/dense/md.inp` — cell parameters (`ABC [angstrom] a b c`)
- `data_example/potential/dense/md.out` — CP2K output with Fermi energies
- `data_example/potential/dense/md-POTENTIAL-v_hartree-1_*.cube` — Hartree potential cube snapshots
- `data_example/potential/distributed/potential_t*_i*/` — per-frame SP subdirs (distributed mode)
- `data_example/bader/bader_work_dir/` — POSCAR, ACF.dat, POTCAR for Bader charge tests
- `data_example/sg/` — CP2K COLVAR restart + LagrangeMultLog files (4 scenarios: angle, distance, combinedCV, more_constrain)

## Project Layout

```
src/md_analysis/
├── __init__.py             # re-exports sub-packages + MDAnalysisError; __version__ = "0.1.0"
├── exceptions.py           # MDAnalysisError — base class for all 17 domain exceptions
├── config.py               # persistent user configuration (~/.config/md_analysis/)
├── main.py                 # programmatic entry points
├── workflows/              # canonical programmatic entry facades — run_* return WorkflowResult
│   ├── __init__.py         #   re-exports all run_* entry points by domain
│   ├── models.py           #   WorkflowResult / MissingArtifactError / require_artifacts_exist
│   ├── water.py            #   run_water_* entry points (density / orientation / three-panel)
│   ├── potential.py        #   run_potential_* entry points (center / Fermi / U vs SHE / phi(z) / sweep)
│   ├── charge.py           #   run_surface_charge / run_tracked_charge / run_counterion_charge
│   ├── calibration.py      #   run_calibration_fit / run_calibration_predict
│   ├── enhanced_sampling.py #   run_slowgrowth_* / run_ti_* entry points
│   ├── scripts.py          #   run_bader_* / run_ti_* / run_potential_* / run_sp_* workdir generators
│   └── composite.py        #   run_interface_analysis (composite: water three-panel + potential full)
├── cli/                    # VASPKIT-style interactive CLI (md-analysis console script)
│   ├── __init__.py         #   main() entry point, banner, top menu, build_menu_tree()
│   ├── _framework.py       #   MenuNode, MenuGroup, MenuCommand, lazy_import()
│   ├── _params.py          #   K key constants, ParamCollector ABC, generic param classes
│   ├── _prompt.py          #   reusable input prompt helpers
│   ├── _water.py           #   water sub-menu (101-105)
│   ├── _potential.py       #   potential sub-menu (211-216)
│   ├── _charge.py          #   charge sub-menu (221-226)
│   ├── _calibration.py     #   calibration sub-menu (231-233)
│   ├── _enhanced_sampling.py #   slow-growth sub-menu (301-302)
│   ├── _constrained_ti.py  #   constrained TI sub-menu (311-313)
│   ├── _scripts.py         #   scripts/tools sub-menu (411-412, 421-422, 431-432, 441-442)
│   └── _settings.py        #   settings sub-menu (900-931)
├── utils/                  # single-frame low-level tools (split into 3 layers)
│   ├── constants.py        #   physical constants, unit conversions, cSHE parameters, defaults
│   ├── formats/            #   single-file parsers, grouped by engine family
│   │   ├── common/         #     engine-neutral formats
│   │   │   └── cube.py     #       cube file I/O, plane-averaged φ(z), slab-averaged potential
│   │   ├── cp2k/           #     CP2K single-file parsers
│   │   │   ├── cell.py     #       cell parameter parsing (.restart + md.inp)
│   │   │   ├── colvar.py   #       COLVAR restart + LagrangeMultLog (ConstraintMetadata / LambdaSeries)
│   │   │   ├── inp.py      #       CP2K .inp parse/generate via cp2k-input-tools (cit tree)
│   │   │   ├── stdout.py   #       md.out / sp.out Fermi-level + step/time parsing
│   │   │   └── xyz.py      #       xyz step-comment parser + per-step Atoms streaming
│   │   ├── vasp/           #     VASP placeholders (NotImplementedError, NOT auto-registered)
│   │   │   ├── report.py   #       VASP REPORT (constraint MD)
│   │   │   ├── outcar.py   #       VASP OUTCAR (Fermi series)
│   │   │   └── locpot.py   #       VASP LOCPOT (planar-averaged potential)
│   │   └── bader/          #     Bader code outputs
│   │       ├── _errors.py  #       BaderParseError (shared)
│   │       ├── acf.py      #       ACF.dat parsing + load_bader_atoms
│   │       └── potcar.py   #       POTCAR ZVAL extraction
│   ├── structure/          #   geometric / chemical-semantics helpers
│   │   ├── cluster.py      #     1D periodic clustering + gap detection
│   │   ├── layer.py        #     metal layer detection, interface identification
│   │   └── water.py        #     water topology, density/orientation/angle profiles
│   └── io/                 #   generic I/O scaffolding (engine-neutral only)
│       ├── _frame_discovery.py  # private: frame directory discovery (bader_t*_i* / potential_t*_i*)
│       └── _io_helpers.py       # private shared I/O helpers (_cumulative_average, _write_csv)
├── engines/                # CP2K / VASP engine facade + neutral dataclasses
│   ├── __init__.py         #   public facade (Protocol, registry, models, CP2K read_*)
│   ├── protocols.py        #   ConstraintMDParser Protocol + parser registry
│   ├── models.py           #   ConstraintMetadata / LambdaSeries / PotentialFrame / FermiRecord
│   ├── cp2k.py             #   CP2KParser + read_constraint_metadata / read_lambda_series /
│   │                       #     read_fermi_series / read_continuous|distributed_potential_frames
│   └── vasp.py             #   VASPParser placeholder (NotImplementedError, NOT auto-registered)
├── water/                  # multi-frame water analysis workflows
│   ├── config.py           #   water analysis defaults + output filename constants
│   ├── Water.py            #   plot_water_three_panel_analysis() — primary entry point
│   ├── _plot.py            #   matplotlib three-panel plotting
│   └── WaterAnalysis/      #   density, orientation, adsorbed-layer sub-workflows
│       ├── _common.py      #     trajectory I/O, per-frame computation, ensemble averaging
│       ├── WaterDensity.py
│       ├── WaterOrientation.py
│       └── AdWaterOrientation.py
├── electrochemical/        # electrochemical analysis grouping package
│   ├── potential/          #   Hartree potential + electrode potential workflows
│   │   ├── CenterPotential.py  # center-slab, Fermi, electrode potential, thickness sweep
│   │   ├── PhiZProfile.py      # phi(z) overlay visualization
│   │   └── _plot.py            # matplotlib plotting helpers
│   ├── charge/             #   Bader charge analysis workflows
│   │   └── Bader/          #     BaderData, SurfaceCharge, AtomCharges, _frame_utils, _plot
│   └── calibration/        #   sigma->phi calibration mapping
│       └── CalibrationWorkflow.py  # calibrate(), predict_potential(), convert_reference()
├── enhanced_sampling/      # enhanced sampling workflows
│   ├── slowgrowth/         #   slow-growth TI (data classes + plots + CSV)
│   └── constrained_ti/     #   constrained TI diagnostics + free energy integration + correction
└── scripts/                # automation scripts
    ├── BaderGen.py         #   VASP Bader work directory generation
    ├── TIGen.py            #   CP2K constrained-MD work directory generation
    ├── PotentialGen.py     #   CP2K SP potential work directory generation (Hartree analysis)
    ├── SpGen.py            #   CP2K SP work directory generation for DeePMD training data
    ├── _frame_selector.py  #   private: shared frame selection for batch Gen scripts
    ├── _inp_tree_edit.py   #   private: workflow-level CP2K inp tree edits (parse -> edit -> generate)
    ├── _potential_template.py  #   private: built-in potential template constants + tree validators
    ├── template/           #   bundled templates (VASP INCAR, KPOINTS; CP2K sp_potential.inp)
    └── utils/
        └── IndexMapper.py  #   CP2K XYZ <-> VASP POSCAR bijective index mapping

test/
├── conftest.py             # shared fixtures
├── unit/
│   ├── utils/              # tests for formats.{common.cube, bader.{acf,potcar}, cp2k.{cell,colvar,stdout,xyz}}
│   │                       #   + structure.{layer,water,cluster} + io._frame_discovery
│   ├── engines/            # engines facade: parser registry, CP2K read_* contract tests
│   ├── cli/                # MenuCommand error handling, settings defaults
│   ├── calibration/        # CalibrationData, Mapper, reference conversion, workflow
│   ├── charge/             # surface charge, atom charges, tracked charges
│   ├── potential/          # single-frame electrode potential pipeline
│   ├── scripts/            # BaderGen, TIGen, IndexMapper
│   ├── enhanced_sampling/  # constrained_ti: ACF, block avg, correction, integration, I/O, workflow
│   ├── workflows/          # run_* facades: WorkflowResult contract tests per domain
│   ├── water/              # empty package placeholder (water tests live under integration/)
│   ├── composite/          # empty package placeholder (composite tests live under workflows/)
│   ├── test_config.py      # persistent user config
│   └── test_logging_setup.py  # NullHandler / StreamHandler setup
└── integration/
    ├── utils/              # water-layer pipeline, full-frame distribution scripts
    ├── water/              # density, orientation, theta, three-panel end-to-end
    ├── potential/          # center slab potential, phi(z) profile with real cube files
    ├── charge/             # trajectory surface charge, surface_charge_analysis end-to-end
    ├── enhanced_sampling/  # constrained TI regression, slow-growth plot
    └── test_main.py        # programmatic entry points (run_* functions)

data_example/               # minimal reproducible input data
├── potential/
│   ├── dense/              #   continuous-mode fixture: md.inp / md.out /
│   │                       #   md-pos-1.xyz / md-POTENTIAL-v_hartree-1_*.cube
│   └── distributed/        #   distributed-mode fixture: per-frame potential_t*_i*/
│                           #   sub-directories with cube + sp.out
├── bader/                  #   bader_work_dir/ with POSCAR, ACF.dat, POTCAR
├── sg/                     #   COLVAR restart + LagrangeMultLog test data (4 scenarios)
├── ti/                     #   8 ti_target_* constraint-point directories for TI regression tests
├── ti_with_correction/     #   9 ti_target_* dirs (restart + LagrangeMultLog + bader/) for TI correction tests
├── calibration/            #   calibration (phi, sigma) CSV data + fitted JSON (counterion / layer)
└── water/                  #   md-pos-1.xyz + md.inp for water analysis
context4agent/              # architecture contracts, decisions, requirements
```

## Running Tests

```bash
# All tests — prefix with PYTHONPATH=src if you used a non-editable
# install, otherwise pytest may import a stale site-packages copy:
PYTHONPATH=src pytest test/

# Single unit test file
PYTHONPATH=src pytest test/unit/utils/test_water_parser.py

# Specific module tests
PYTHONPATH=src pytest test/unit/utils/test_slowgrowth_parser.py   # cp2k.colvar tests
PYTHONPATH=src pytest test/unit/charge/test_charge_analysis.py     # Bader charge tests
PYTHONPATH=src pytest test/integration/                             # all integration tests

# With an editable install (`pip install -e .`) the PYTHONPATH=src
# prefix is unnecessary.
```

## Architecture Notes

- **Three-layer design**: `utils` (single-frame primitives) → `water` / `electrochemical` / `enhanced_sampling` (multi-frame workflows) → `main` / `cli` (integration entry points)
- **Interface detection**: 1D periodic clustering on metal z-coordinates, largest gap = water region, two bounding layers = interfaces
- **Water profiles**: computed from selected interface toward the cell midpoint, normalized to [0, 1], then ensemble-averaged across frames
- **Electrode potential**: U = -E_Fermi + φ_center + ΔΨ_a(H₃O⁺/w) - μ(H⁺,g⁰) - ΔE_ZP (computational SHE)
- **Surface charge**: two methods — `counterion` (non-water non-metal species only) and `layer` (interface-layer metal atom net charge sum)
- **Frame selection**: `frame_start` / `frame_end` / `frame_step` supported across all workflows (prompted via "advanced parameters" in the CLI)
- **Bader workdir generation**: reads CP2K XYZ trajectory → writes VASP POSCAR (element-grouped order) with bijective index map in comment line → copies bundled INCAR/KPOINTS templates → optionally generates POTCAR via vaspkit. Batch mode names directories `bader_t{time_fs}_i{step}` from frame metadata

## For Contributors

Architecture contracts and implementation rules:

- `context4agent/architecture/modules/data_contract.md` — output shapes, units, CSV headers
- `context4agent/architecture/modules/glossary_units.md` — terminology and unit definitions
- `context4agent/architecture/modules/src/**/interface_exposure.md` — public API per module
- `context4agent/architecture/modules/src/**/implementation_guidelines.md` — implementation rules
