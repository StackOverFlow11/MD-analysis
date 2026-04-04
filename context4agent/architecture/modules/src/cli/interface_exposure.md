# `md_analysis.cli` interface exposure

> Source: `src/md_analysis/cli/`

## Public interface

- `main()` -- VASPKIT-style interactive CLI entry point, registered as `md-analysis` console script via `pyproject.toml`.

## Internal modules (non-public)

- `_prompt.py` -- low-level prompt helpers (`prompt_str`, `prompt_str_required`, `prompt_int`, `prompt_float`, `prompt_choice`, `prompt_bool`, `set_input_source`)
- `_params.py` -- parameter collection layer (`K` key constants, `ParamCollector` ABC, generic param classes: `StrParam`, `FloatParam`, `IntParam`, `ChoiceParam`, `FixedParam`, `ConfigDefaultParam`; special: `CellAbcParam`, `MetalElementsParam`, `FrameSliceParam`)
- `_water.py` -- water analysis sub-menu (codes 101-105) + parameter collection + dispatch
- `_potential.py` -- potential analysis sub-menu (codes 211-216) + input_mode selection (continuous/distributed) + parameter collection + dispatch
- `_charge.py` -- charge analysis sub-menu (codes 221-226) + parameter collection + dispatch
- `_calibration.py` -- calibration sub-menu (codes 231-233): CSV/manual calibrate + predict
- `_enhanced_sampling.py` -- slow-growth sub-group 30 (codes 301-302) + file discovery helpers
- `_constrained_ti.py` -- constrained TI sub-group 31 (codes 311-313): single-point diagnostics + full TI analysis + constant-potential correction (lazy import of `enhanced_sampling.constrained_ti`)
- `_scripts.py` -- scripts/tools sub-menu; Bader sub-group 41 (codes 411-412), TI sub-group 42 (codes 421-422), SP Potential sub-group 43 (codes 431-432)
- `_settings.py` -- settings sub-menu (codes 901-910) + persistent config management (903-906: configurable analysis defaults, 907: reset all defaults, 909: potential output reference, 910: SP inp template path)

## Menu structure

- Top menu: 1=Water, 2=Electrochemical, 3=Enhanced Sampling, 4=Scripts/Tools, 9=Settings, 0=Exit
- Electrochemical sub-groups: 21=Potential, 22=Charge, 23=Calibration
- Enhanced Sampling sub-groups: 30=Slow-Growth, 31=Constrained TI Analysis
- Scripts sub-groups: 41=Bader Charge Preparation, 42=Thermodynamic Integration Preparation
- Leaf codes use numbered codes: 1xx=Water, 21x=Potential, 22x=Charge, 23x=Calibration, 30x=Slow-Growth, 31x=Constrained TI, 41x=Bader, 42x=TI, 9xx=Settings
- Flat index allows direct jump to any leaf code from the root menu
- After one analysis completes, control returns to the parent menu (MenuGroup uses a while-True loop)

## Entry point registration

```toml
[project.scripts]
md-analysis = "md_analysis.cli:main"
```
