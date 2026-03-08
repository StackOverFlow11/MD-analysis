# md_analysis.scripts — Implementation Guidelines

## Responsibilities

- Provide automation scripts for MD-analysis workflows (Bader work directory generation).
- `BaderGen.py`: generate VASP single-point work directories for Bader charge analysis (single frame and batch from trajectory).

## Dependencies

- `md_analysis.scripts.utils.IndexMapper`: `compute_index_map`, `write_poscar_with_map`
- `md_analysis.config`: `get_config`, `KEY_VASP_SCRIPT_PATH` (persistent user configuration)
- `ase.io.iread`: trajectory reading for batch generation
- `importlib.resources`: template file access (`INCAR`, `KPOINTS`)
- `shutil`: file copy, `which()` for vaspkit detection
- `subprocess`: vaspkit invocation for POTCAR generation
- `tqdm`: optional progress bar for batch generation

## Package Structure

```
scripts/
  __init__.py             # re-exports BaderGenError, generate_bader_workdir, batch_generate_bader_workdirs
  BaderGen.py             # core implementation (single + batch)
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

## Sync Rules

Changes to this module must update:
- This file (`implementation_guidelines.md`)
- `interface_exposure.md` in the same directory
- `CLAUDE.md` project map
