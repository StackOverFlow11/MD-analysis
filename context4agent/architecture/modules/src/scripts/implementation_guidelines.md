# md_analysis.scripts — Implementation Guidelines

## Responsibilities

- Provide automation scripts for MD-analysis workflows (Bader work directory generation, TI constrained-MD setup, SP potential preparation, DP training data SP preparation).
- `BaderGen.py`: generate VASP single-point work directories for Bader charge analysis (single frame and batch from trajectory).
- `TIGen.py`: generate CP2K constrained-MD work directories for thermodynamic integration (single target + batch).
- `PotentialGen.py`: generate CP2K single-point work directories for Hartree potential analysis (single frame and batch from trajectory).
- `SpGen.py`: generate CP2K single-point work directories for DeePMD training data collection. Mirrors `PotentialGen` but uses distinct config key (`KEY_DP_SP_INP_TEMPLATE_PATH`) and directory prefix (`sp_t*_i*`) to avoid mixing with potential analysis. Exposed via the `workflows.scripts.run_sp_batch` facade.
- `_inp_tree_edit.py`: private shared workflow-level CP2K inp **tree** editing helpers (`ensure_topology_init_xyz`, `set_project`, `set_md_steps`, `set_collective_target`, `zero_target_growth`, `ensure_kind_blocks`). Operate in place on cit trees produced by `utils.formats.cp2k.inp.parse_inp`; used by `TIGen`, `PotentialGen` and `SpGen` ("parse → tree edit → generate" write side).
- `_potential_template.py`: built-in default potential template constants (`DEFAULT_BASIS_SET` / `DEFAULT_POTENTIAL` / `DEFAULT_TEMPLATE_RESOURCE`) + tree-driven template validators (`extract_kind_elements`, `has_smear`, `has_v_hartree_cube`), all operating on cit trees.
- `_frame_selector.py`: private shared trajectory slicing module. Provides `FrameSelection` dataclass, `iter_selected_frames` iterator, and `resolve_single_frame` helper. Used by **all three batch Gen functions** (Bader/Potential/SpGen) to unify index-based and time-based frame selection. TIGen is excluded because its "time mode" is target-based (linspace → CV → nearest frame), not range-based slicing.

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
- `md_analysis.engines.cp2k`: `read_constraint_metadata_from_restart` (restart metadata for CV reconstruction)
- `md_analysis.utils.formats.cp2k.inp`: `parse_inp` / `generate_inp` (cit tree parse + generate)
- `md_analysis.scripts._inp_tree_edit`: `set_project`, `set_md_steps`, `set_collective_target`, `zero_target_growth`, `ensure_topology_init_xyz`
- `md_analysis.utils.constants`: `AU_TIME_TO_FS`
- `ase.io`: `iread` (trajectory reading), `write` (init.xyz output)
- `tqdm`: optional progress bar for batch generation

## Package Structure

```
scripts/
  __init__.py             # re-exports Bader + TI + Potential + SpGen public symbols
  BaderGen.py             # Bader work directory generation (single + batch)
  TIGen.py                # TI constrained-MD work directory generation (single + batch)
  PotentialGen.py         # SP potential work directory generation (single + batch)
  SpGen.py                # DP training SP work directory generation (single + batch)
  _inp_tree_edit.py       # Private shared CP2K inp tree-edit helpers (workflow semantics)
  _potential_template.py  # Private: built-in potential template constants + tree validators
  _frame_selector.py      # Private shared trajectory slicing (FrameSelection + iter/resolve)
  template/
    __init__.py           # empty (importlib.resources requirement)
    INCAR                 # VASP INCAR template (single-point, Bader settings)
    KPOINTS               # VASP KPOINTS template (Gamma-only)
    sp_potential.inp      # Built-in default CP2K potential-SP template (parseable, no placeholders)
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
2. **inp modification via cit tree editing**: Parses the SG inp once into a cit tree (`inp.parse_inp`), then per target deep-copies the base tree and edits canonical keys — `project_name` → `cMD`, `motion.md.steps` → user value (tree path touches only `&MD`, never `MAX_SCF`-style same-name keywords), `target` → snapped a.u. as a bare internal-unit value in the `&COLLECTIVE` block selected by `colvar_id` (`None` = primary / first block, matching reader-side `colvars.primary`; pre-refactor behavior wrote all blocks), `target_growth` → `0` in every `&COLLECTIVE` block — then regenerates via `inp.generate_inp`. Explicit `[unit]` tokens in the template are converted to internal units by the inp.py preprocessing (W2).
3. **Units**: All CV values in atomic units (CP2K default). No custom unit conversion supported — avoids complex dimensionality issues with exotic CVs.
4. **Batch efficiency**: Pre-loads trajectory and parses restart once; reuses across all target points.
5. **Two batch modes**: Numeric (explicit a.u. values) or time-based (linspace over fs range → map to CV values).

## PotentialGen Key Design Decisions

1. **Three-level template resolution**: explicit `inp_template_path` argument > `KEY_SP_INP_TEMPLATE_PATH` persistent config > built-in default potential template bundled at `scripts/template/sp_potential.inp` (plain parseable CP2K input, no placeholders). User templates are parsed to a cit tree and validated against the potential-SP hard requirements (`&V_HARTREE_CUBE` present, `&SMEAR` present, `&KIND` blocks cover every frame element — validators in `_potential_template.py`, tree-driven). The built-in template instead gets its `&KIND` blocks upserted from the frame's element composition via `ensure_kind_blocks` (defaults `DZVP-MOLOPT-SR-GTH` / `GTH-PBE`), so full coverage holds by construction.
2. **CELL ABC replacement via tree edit**: `inp.set_cell_abc` writes `force_eval.subsys.cell.abc` on the parsed tree; cell values come from restart/md.inp parsing. Vector-form (A/B/C) templates are not modified — warning logged, no-op (legacy semantics preserved).
3. **TOPOLOGY enforcement**: Same pattern as TIGen — ensures `COORD_FILE_NAME init.xyz` + `COORD_FILE_FORMAT XYZ` (section created when missing).
4. **Batch efficiency (generate-once)**: Parses the template once, applies the frame-independent edits once (cell + topology) and generates the sp.inp text once; the same text is written to every subdirectory (cell is constant across frames).
5. **Directory naming**: `potential_t{time}_i{step}` matches `_frame_source._SP_DIR_RE` for seamless downstream analysis with `input_mode="distributed"`.
6. **Shared helpers**: Tree edits delegated to `_inp_tree_edit.py` (`ensure_topology_init_xyz`, `ensure_kind_blocks`); parse/generate via `utils.formats.cp2k.inp`.

## SpGen Key Design Decisions

1. **Separate config key**: Uses `KEY_DP_SP_INP_TEMPLATE_PATH` (not `KEY_SP_INP_TEMPLATE_PATH`), so users can persist both a Hartree-potential template (containing `V_HARTREE CUBE`) and a DP training template (containing `PRINT FORCES`) and select them via independent CLI menus (913 vs 914).
2. **Distinct directory prefix**: `sp_t{time}_i{step}` instead of `potential_t{time}_i{step}`. This is critical because `electrochemical/potential/_frame_source.py:discover_distributed_frames()` defaults to `dir_pattern="potential_t*_i*"`; if SpGen reused that prefix, a mixed directory would cause accidental consumption of DP training dirs during potential analysis.
3. **`_frame_discovery.FRAME_DIR_STEP_TIME_RE` compatibility**: The generic `_t\d+_i\d+` regex matches any prefix, so `sp_t*_i*` still works with downstream generic frame discovery utilities — only the electrochemical default pattern is `potential_t*_i*`.
4. **Code reuse via `_inp_tree_edit.py` + `utils.formats.cp2k.inp`**: The same tree-edit helpers (`ensure_topology_init_xyz`, `inp.set_cell_abc`) are shared with `PotentialGen`; each module keeps its own workdir orchestration (single + batch) but the "parse → tree edit → generate" write path is factored out. Batch callers parse once, edit once, generate once and reuse the text (frame-independent edits).
5. **Workflow facade backend**: `SpGen` / `BaderGen` / `TIGen` batch wrappers are invoked through `workflows.scripts.run_sp_batch` / `run_bader_batch` / `run_ti_batch`; `PotentialGen` is CLI-only. Workdir paths are returned on `WorkflowResult.artifacts`.

## Unified Frame Selection Key Design Decisions (`_frame_selector.py`)

1. **Explicit `mode` discriminator**: rather than inferring mode from "which params are None" (implicit), the module uses an explicit `mode: Literal["index", "time"]` field. This keeps schema enum representation precise (`{"mode": {"enum": ["index", "time"]}}`) instead of requiring callers to guess which parameter combinations are valid.
2. **Python default `mode="index"`, CLI default `"time"`**: Python API keeps `"index"` for backward compatibility with existing code/tests. CLI `frame_mode = ChoiceParam(..., default="time")` reflects the user's preference for time-based slicing as the common case in DP training workflows.
3. **Strict time mode validation**: `FrameSelection.__post_init__` requires all three `time_start_fs/end_fs/step_fs` to be set together when `mode="time"`. Partial configuration raises `FrameSelectionError` immediately (fail-fast). This avoids silent fallbacks that could produce unexpected frame selections.
4. **Greedy time matching**: single-pass iteration. At each target `t_k`, yield the first frame with `time >= t_k`, advance `next_target = time + t_step`. Simpler than closest-neighbor matching and requires no look-ahead. For DP training data collection, greedy behavior is predictable enough.
5. **Nearest-neighbor for single frame**: `resolve_single_frame` in time mode uses a different algorithm (closest by `|time - target|` with early break on monotonic time). This is because for single-frame commands, "exactly time t" is semantically clearer than "first frame after t". Differences from the requested time are reported as warnings (CLI: logger.warning + stdout).
6. **Missing `atoms.info["time"]` → early fail**: Time mode requires CP2K-format XYZ comment lines. Missing metadata raises `FrameSelectionError` rather than falling back to index mode, because silent fallback could mask data corruption.
7. **TIGen explicitly excluded**: TIGen's time mode maps `linspace(t_init, t_final, n_points)` to nearest CV values (target-based point selection), not range slicing. Forcing a common abstraction would over-engineer both modules.
8. **Private module with public-like semantics**: `_frame_selector.py` is prefixed with underscore (not re-exported from `scripts/__init__.py`), but the `FrameSelection` / `iter_selected_frames` / `resolve_single_frame` API is stable and exposed via batch function signatures.

## Sync Rules

Changes to this module must update:
- This file (`implementation_guidelines.md`)
- `interface_exposure.md` in the same directory
- `AGENTS.md` project map
