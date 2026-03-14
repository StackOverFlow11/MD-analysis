# md_analysis.scripts — Implementation Guidelines

## Responsibilities

- Provide automation scripts for MD-analysis workflows (Bader work directory generation, TI constrained-MD setup).
- `BaderGen.py`: generate VASP single-point work directories for Bader charge analysis (single frame and batch from trajectory).
- `TIGen.py`: generate CP2K constrained-MD work directories for thermodynamic integration (single target + batch).

## Dependencies

### BaderGen
- `md_analysis.scripts.utils.IndexMapper`: `compute_index_map`, `write_poscar_with_map`
- `md_analysis.config`: `get_config`, `KEY_VASP_SCRIPT_PATH` (persistent user configuration)
- `ase.io.iread`: trajectory reading for batch generation
- `importlib.resources`: template file access (`INCAR`, `KPOINTS`)
- `shutil`: file copy, `which()` for vaspkit detection
- `subprocess`: vaspkit invocation for POTCAR generation
- `tqdm`: optional progress bar for batch generation

### TIGen
- `md_analysis.utils.RestartParser.ColvarParser`: `parse_colvar_restart`, `ColvarRestart`
- `md_analysis.utils.config`: `AU_TIME_TO_FS`
- `ase.io`: `iread` (trajectory reading), `write` (init.xyz output)
- `tqdm`: optional progress bar for batch generation

## Package Structure

```
scripts/
  __init__.py             # re-exports Bader + TI public symbols
  BaderGen.py             # Bader work directory generation (single + batch)
  TIGen.py                # TI constrained-MD work directory generation (single + batch)
  template/
    __init__.py           # empty (importlib.resources requirement)
    INCAR                 # VASP INCAR template (single-point, Bader settings)
    KPOINTS               # VASP KPOINTS template (Gamma-only)
  utils/
    __init__.py           # re-exports all IndexMapper public symbols
    IndexMapper.py        # bijective index mapping
```

## Key Design Decisions

1. **Single-frame + Batch**: `generate_bader_workdir` handles one frame; `batch_generate_bader_workdirs` iterates over a CP2K XYZ trajectory and calls the single-frame function per frame. Batch accepts `cell_abc` tuple (decoupled from cell source) and uses XYZ comment-line metadata (`atoms.info['i']`, `atoms.info['time']`) for directory naming (`bader_t{time}_i{step}`).
2. **Template packaging**: INCAR/KPOINTS are bundled as package data via `importlib.resources`, accessed with `as_file` context manager.
3. **POTCAR via vaspkit**: Uses `subprocess.run(["vaspkit"], input="103\n")` with timeout. Requires vaspkit in PATH and VASP pseudopotential directory configured.
4. **Config fallback**: Script path falls back to persistent config (`~/.config/md_analysis/config.json`) when not explicitly provided.

## TIGen Key Design Decisions

1. **Frame snapping**: For a given target CV value, finds the trajectory frame whose reconstructed CV (from `compute_target_series`) is closest; uses that frame's actual CV as the snapped TARGET.
2. **inp modification via regex**: Modifies PROJECT, STEPS (in `&MD` only), TARGET/TARGET_GROWTH (in `&COLLECTIVE`, strips `[unit]` annotations), and ensures `&TOPOLOGY` has `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ`.
3. **Units**: All CV values in atomic units (CP2K default). No custom unit conversion supported — avoids complex dimensionality issues with exotic CVs.
4. **Batch efficiency**: Pre-loads trajectory and parses restart once; reuses across all target points.
5. **Two batch modes**: Numeric (explicit a.u. values) or time-based (linspace over fs range → map to CV values).

## Sync Rules

Changes to this module must update:
- This file (`implementation_guidelines.md`)
- `interface_exposure.md` in the same directory
- `CLAUDE.md` project map
