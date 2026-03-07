# `md_analysis.cli` interface exposure

> Source: `src/md_analysis/cli/`

## Public interface

- `main()` -- VASPKIT-style interactive CLI entry point, registered as `md-analysis` console script via `pyproject.toml`.

## Internal modules (non-public)

- `_prompt.py` -- reusable input-prompt helpers (`_prompt_str`, `_prompt_int`, `_prompt_float`, `_prompt_choice`, `_prompt_bool`, `_parse_metal_elements`, `_prompt_global_params`)
- `_water.py` -- water analysis sub-menu (codes 101-105) + parameter collection + dispatch
- `_potential.py` -- potential analysis sub-menu (codes 201-206) + parameter collection + dispatch
- `_charge.py` -- charge analysis sub-menu (codes 301-303) + parameter collection + dispatch
- `_scripts.py` -- scripts/tools sub-menu (code 401) + Bader work directory generation
- `_settings.py` -- settings sub-menu (codes 901-902) + persistent config management

## Menu structure

- Top menu: 1=Water, 2=Potential, 3=Charge, 4=Scripts/Tools, 9=Settings, 0=Exit
- Sub-menus use numbered codes: 1xx=Water, 2xx=Potential, 3xx=Charge, 4xx=Scripts, 9xx=Settings
- After one analysis completes, the program exits (no loop)

## Entry point registration

```toml
[project.scripts]
md-analysis = "md_analysis.cli:main"
```
