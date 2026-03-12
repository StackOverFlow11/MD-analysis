# `md_analysis.cli` interface exposure

> Source: `src/md_analysis/cli/`

## Public interface

- `main()` -- VASPKIT-style interactive CLI entry point, registered as `md-analysis` console script via `pyproject.toml`.

## Internal modules (non-public)

- `_prompt.py` -- reusable input-prompt helpers (`_get_effective_default`, `_prompt_str`, `_prompt_int`, `_prompt_float`, `_prompt_choice`, `_prompt_bool`, `_parse_metal_elements`, `_prompt_global_params`)
- `_water.py` -- water analysis sub-menu (codes 101-105) + parameter collection + dispatch
- `_potential.py` -- potential analysis sub-menu (codes 211-216) + parameter collection + dispatch
- `_charge.py` -- charge analysis sub-menu (codes 221-223) + parameter collection + dispatch
- `_enhanced_sampling.py` -- enhanced sampling sub-menu (codes 301-302) + slow-growth plot commands (lazy import of `enhanced_sampling.slowgrowth`)
- `_scripts.py` -- scripts/tools sub-menu (codes 401-402) + single-frame and batch Bader work directory generation
- `_settings.py` -- settings sub-menu (codes 901-907) + persistent config management (903-906: configurable analysis defaults, 907: reset all defaults)

## Menu structure

- Top menu: 1=Water, 2=Electrochemical, 3=Enhanced Sampling, 4=Scripts/Tools, 9=Settings, 0=Exit
- Electrochemical sub-groups: 21=Potential, 22=Charge
- Leaf codes use numbered codes: 1xx=Water, 21x=Potential, 22x=Charge, 3xx=Enhanced Sampling, 4xx=Scripts, 9xx=Settings
- Flat index allows direct jump to any leaf code from the root menu
- After one analysis completes, the program exits (no loop)

## Entry point registration

```toml
[project.scripts]
md-analysis = "md_analysis.cli:main"
```
