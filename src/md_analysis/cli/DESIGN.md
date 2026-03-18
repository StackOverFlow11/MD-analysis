# CLI Composite + Command Refactoring Design

> **NOTE (2026-03-18)**: 本文档中的 `output_subdir` 硬编码模式已被 `output_name` 树推导机制取代。
> 菜单编号已更新：Charge 221-225、Calibration 231-233、Bader 411-412、TI 421-422。
> 最新约定参见 `cli/CLAUDE.md`。

## Design Goals

1. **Composite** tree structure: one `group.add(Cmd(...))` = menu text + dispatch + registration
2. **Command** pattern: leaf nodes encapsulate param collection + execution
3. Eliminate duplicated parameter collection across water/potential/charge/scripts
4. Unified error handling in base class, not per-handler decorators
5. VASPKIT-style interactive loop (stay in menu after command completes)
6. Support typing leaf code (e.g. `101`) directly at root menu

## File Structure

```
cli/
├── __init__.py       # main(), build_menu_tree(), logging setup, banner
├── _framework.py     # MenuNode, MenuGroup, MenuCommand base classes, lazy_import()
├── _prompt.py        # Low-level prompt functions (prompt_str, prompt_int, etc.)
├── _params.py        # ParamCollector hierarchy + predefined instances + key constants
├── _water.py         # 5 Water command classes
├── _potential.py     # 6 Potential command classes
├── _charge.py        # 3 Charge command classes
├── _scripts.py       # 2 Scripts command classes
└── _settings.py      # 7 Settings command classes
```

---

## Core Framework (`_framework.py`)

### MenuNode (ABC)

```python
class MenuNode(ABC):
    """Abstract base for all menu tree nodes."""

    def __init__(self, code: str, label: str) -> None:
        self.code = code
        self.label = label

    @abstractmethod
    def run(self) -> None: ...
```

`run()` returns `None` — **not** an exit code. This eliminates the "return value
semantics" ambiguity. `MenuGroup` loops internally; `MenuCommand` handles its own
errors. Neither propagates success/failure upward.

### MenuGroup

```python
class MenuGroup(MenuNode):
    """Non-leaf node: displays a menu and dispatches to children."""

    def __init__(self, code: str, label: str) -> None:
        super().__init__(code, label)
        self.children: list[MenuNode | str] = []   # str = separator line
        self._flat_index: dict[str, MenuCommand] | None = None

    def add(self, *nodes: MenuNode | str) -> None:
        for node in nodes:
            self.children.append(node)

    def run(self) -> None:
        while True:
            print(self._render_menu())
            choice = input(" Input: ").strip()

            if choice == "0":
                return

            # 1) Match direct child by code
            for child in self.children:
                if isinstance(child, MenuNode) and child.code == choice:
                    child.run()
                    break
            else:
                # 2) Fallback: flat index for leaf shortcut (root only)
                if self._flat_index and choice in self._flat_index:
                    self._flat_index[choice].run()
                else:
                    print(f"\n Invalid choice: {choice!r}")

    def build_flat_index(self) -> None:
        """Recursively index all descendant MenuCommand nodes by code.
        Call once on root after tree is fully assembled."""
        self._flat_index = {}
        self._collect_commands(self._flat_index)

    def _collect_commands(self, index: dict[str, MenuCommand]) -> None:
        for child in self.children:
            if isinstance(child, MenuCommand):
                index[child.code] = child
            elif isinstance(child, MenuGroup):
                child._collect_commands(index)

    def _render_menu(self) -> str:
        lines = [f"\n ---------- {self.label} ----------\n"]
        for child in self.children:
            if isinstance(child, str):
                if child:
                    lines.append(f"\n --- {child} ---")
                else:
                    lines.append("")                   # blank separator
            else:
                lines.append(f" {child.code}) {child.label}")
        lines.append("\n   0) Back / Exit\n")
        return "\n".join(lines)
```

Key behaviors:
- **Loops** until user enters "0". After any child command completes (success or
  failure), control returns to the loop and the menu is re-displayed.
- **Flat index**: built once via `build_flat_index()` on root. Enables typing "101"
  at the root menu to jump directly to a leaf command. Sub-menus only match their
  direct children (no flat fallback).
- **Separators**: `str` items in children render as visual dividers (empty string =
  blank line, non-empty = labeled separator).

### MenuCommand

```python
class MenuCommand(MenuNode):
    """Leaf node: collects parameters and executes an analysis command."""

    params: tuple[ParamCollector, ...] = ()
    advanced_params: tuple[ParamCollector, ...] = ()
    output_subdir: str = ""        # e.g. "water", "potential/center"

    def run(self) -> None:
        try:
            ctx = self._collect_all_params()
            if self.output_subdir:
                outdir = Path(ctx.get(K.OUTDIR, "analysis")) / self.output_subdir
                outdir.mkdir(parents=True, exist_ok=True)
                ctx[K.OUTDIR_RESOLVED] = outdir
            self.execute(ctx)
        except (MDAnalysisError, FileNotFoundError, ValueError, RuntimeError) as exc:
            print(f"\n  Error: {exc}")
        except Exception as exc:
            logger.error("Unexpected error in %s: %s", self.label, exc, exc_info=True)
            print(f"\n  Unexpected error ({type(exc).__name__}): {exc}")

    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        for p in self.params:
            p.collect(ctx)
        if self.advanced_params:
            if prompt_bool("Modify advanced parameters?", default=False):
                for p in self.advanced_params:
                    p.collect(ctx)
            else:
                for p in self.advanced_params:
                    p.apply_default(ctx)
        return ctx

    @abstractmethod
    def execute(self, ctx: dict) -> None:
        """Subclasses implement pure business logic here."""
        ...
```

Key behaviors:
- **Error handling** in `run()` replaces `@_handle_cmd_error` decorator. All
  commands get it for free.
- **Output dir creation** in `run()` replaces per-handler boilerplate.
- **`params` and `advanced_params`** are class-level **tuples** (immutable, no
  shared-mutable-state hazard). Subclasses that need per-instance customization
  set `self.params` in `__init__()`, which shadows the class attribute.
- **`_collect_all_params()`** can be overridden for commands with non-standard
  collection flows (e.g. scripts that display trajectory info mid-collection).
- **`execute(ctx)`** receives a dict and does only business logic — no error
  handling, no dir creation, no prompting (normally).

### lazy_import()

```python
def lazy_import(module_path: str, name: str):
    """Lazily import a callable from a dotted module path."""
    from importlib import import_module
    return getattr(import_module(module_path), name)
```

Provides a consistent pattern for delayed imports in `execute()`:

```python
def execute(self, ctx):
    analyze = lazy_import("md_analysis.water", "water_mass_density_z_distribution_analysis")
    csv = analyze(xyz_path=ctx[K.XYZ], ...)
```

Benefits: explicit full module path (greppable), consistent style across all
commands, avoids scattered `from ..water import ...` with relative paths.

---

## Parameter Collection (`_params.py`)

### Key Constants

```python
class K:
    """Parameter dict key constants. Prevents string-key typos."""
    XYZ = "xyz"
    CELL_ABC = "cell_abc"
    DZ_A = "dz_A"
    LAYER_TOL = "layer_tol"
    OUTDIR = "outdir"
    OUTDIR_RESOLVED = "outdir_resolved"
    FRAME_START = "frame_start"
    FRAME_END = "frame_end"
    FRAME_STEP = "frame_step"
    CUBE_PATTERN = "cube_pattern"
    THICKNESS = "thickness"
    CENTER_MODE = "center_mode"
    METAL_ELEMENTS = "metal_elements"
    MD_OUT = "md_out"
    FERMI_UNIT = "fermi_unit"
    MAX_CURVES = "max_curves"
    THICKNESS_END = "thickness_end"
    METHOD = "method"
    ROOT_DIR = "root_dir"
    DIR_PATTERN = "dir_pattern"
    NORMAL = "normal"
    N_SURFACE_LAYERS = "n_surface_layers"
    # scripts-specific
    FRAME = "frame"
    WORKDIR_NAME = "workdir_name"
    SCRIPT_PATH = "script_path"
    GEN_POTCAR = "gen_potcar"
```

All `execute()` methods and param collectors reference `K.XYZ` instead of `"xyz"`.
IDE autocomplete + find-references + typo detection at no runtime cost.

### ParamCollector (ABC)

```python
class ParamCollector(ABC):
    """One reusable unit of parameter collection."""

    @abstractmethod
    def collect(self, ctx: dict) -> None:
        """Prompt user and store result(s) in ctx."""
        ...

    @abstractmethod
    def apply_default(self, ctx: dict) -> None:
        """Store default value(s) in ctx without prompting."""
        ...
```

### Generic Param Classes (~80% of params)

```python
class StrParam(ParamCollector):
    """Single string prompt with optional default."""
    def __init__(self, key: str, label: str, *, default: str | None = None,
                 required: bool = False): ...

class FloatParam(ParamCollector):
    def __init__(self, key: str, label: str, *, default: float): ...

class IntParam(ParamCollector):
    def __init__(self, key: str, label: str, *, default: int | None = None): ...

class ChoiceParam(ParamCollector):
    def __init__(self, key: str, label: str, choices: list[str], *,
                 default: str): ...

class FixedParam(ParamCollector):
    """Injects a fixed value without prompting. Both collect() and
    apply_default() store the same value."""
    def __init__(self, key: str, value: Any): ...

class ConfigDefaultParam(ParamCollector):
    """Reads default from persistent user config (~/.config/md_analysis/).
    collect() prompts with config-aware default; apply_default() uses it silently."""
    def __init__(self, key: str, label: str, *, config_key: str): ...
```

### Special Param Classes (~20% of params)

```python
class CellAbcParam(ParamCollector):
    """Cell parameter acquisition with .restart/md.inp selection + retry."""

class MetalElementsParam(ParamCollector):
    """Comma-separated input -> set[str] | None."""

class FrameSliceParam(ParamCollector):
    """Collects frame_start, frame_end, frame_step as a group."""
    # collect() prompts 3 values; apply_default() sets all to None
```

### Predefined Instances

```python
# Required params (always in `params`, never in `advanced_params`)
xyz_path     = StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True)
xyz_path_opt = StrParam(K.XYZ, "XYZ trajectory file", default="md-pos-1.xyz")
cell_abc     = CellAbcParam()
cube_pattern = StrParam(K.CUBE_PATTERN, "Cube file glob pattern",
                        default="md-POTENTIAL-v_hartree-1_*.cube")
md_out_path  = StrParam(K.MD_OUT, "CP2K md.out file", default="md.out")

# Typed params
thickness     = FloatParam(K.THICKNESS, "Slab averaging thickness (A)", default=7.5)
center_mode   = ChoiceParam(K.CENTER_MODE, "Slab center mode",
                            ["interface", "cell"], default="interface")
fermi_unit    = ChoiceParam(K.FERMI_UNIT, "Fermi energy unit in md.out",
                            ["au", "ev"], default="au")
max_curves    = IntParam(K.MAX_CURVES, "Max curves on overlay (0=all)", default=0)
thickness_end = FloatParam(K.THICKNESS_END, "Thickness sweep upper limit (A)", default=15.0)
normal_axis   = ChoiceParam(K.NORMAL, "Surface normal axis",
                            ["a", "b", "c"], default="c")
charge_method = ChoiceParam(K.METHOD, "Charge method",
                            ["counterion", "layer"], default="counterion")
n_surface_layers = IntParam(K.N_SURFACE_LAYERS,
                            "Number of surface layers per interface", default=1)

# Config-backed params
dz_bin    = ConfigDefaultParam(K.DZ_A, "Z-axis bin width (A)", config_key=KEY_Z_BIN_WIDTH_A)
layer_tol = ConfigDefaultParam(K.LAYER_TOL, "Layer clustering tolerance (A)",
                               config_key=KEY_LAYER_TOL_A)

# Structural params
metal_elements = MetalElementsParam()
root_dir       = StrParam(K.ROOT_DIR, "Root directory", default=".")
dir_pattern    = StrParam(K.DIR_PATTERN, "Frame subdirectory pattern", default="bader_t*_i*")
outdir         = StrParam(K.OUTDIR, "Output directory", default="analysis")
frame_slice    = FrameSliceParam()
```

Total: **6 generic classes + 3 special classes + ~18 predefined instances**.
Replaces duplicated prompt code scattered across 4 files.

---

## Input Layer (`_prompt.py`)

Low-level prompt functions. Similar to current, with two changes:

### 1. Drop underscore prefix

These are internal to the `cli` package, but used directly by `_params.py` and
command `execute()` methods. Rename for readability:

- `_prompt_str` -> `prompt_str`
- `_prompt_int` -> `prompt_int`
- `_prompt_float` -> `prompt_float`
- `_prompt_choice` -> `prompt_choice`
- `_prompt_bool` -> `prompt_bool`
- `_prompt_str_required` -> `prompt_str_required`

### 2. Testable input source

```python
_input_fn: Callable[[str], str] = input

def set_input_source(fn: Callable[[str], str]) -> None:
    """Replace the input function (for testing)."""
    global _input_fn
    _input_fn = fn

def _read(prompt_text: str) -> str:
    """All prompt functions call this instead of input() directly."""
    return _input_fn(prompt_text)
```

All `prompt_*` functions and `MenuGroup.run()` use `_read()` internally. Tests call
`set_input_source()` to inject scripted responses without patching `builtins.input`.

### Removed from `_prompt.py`

- `_handle_cmd_error` -> replaced by `MenuCommand.run()` error handling
- `_get_effective_default` -> absorbed into `ConfigDefaultParam`
- `_parse_metal_elements` -> absorbed into `MetalElementsParam`
- `_prompt_global_params` -> replaced by `FrameSliceParam` + `outdir`
- `_prompt_cell_abc` -> absorbed into `CellAbcParam`

---

## Command Definitions

### Water (`_water.py`)

5 commands, all sharing the same param signature:

```python
_WATER_PARAMS = (xyz_path, cell_abc, dz_bin)
_WATER_ADVANCED = (layer_tol, outdir, frame_slice)

class WaterDensityCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = "water"

    def execute(self, ctx):
        analyze = lazy_import("md_analysis.water",
                              "water_mass_density_z_distribution_analysis")
        csv = analyze(
            xyz_path=ctx[K.XYZ], cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED], dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START], frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print(f"\n Analysis complete. Output:\n   density_csv: {csv}")
```

`WaterOrientationCmd`, `AdWaterOrientationCmd`, `AdWaterThetaCmd` follow the same
pattern, differing only in `execute()`. `WaterThreePanelCmd` calls
`run_water_analysis` from `md_analysis.main`.

### Potential (`_potential.py`)

Each command declares its own param subset (replaces conditional `_collect_params_for(code)`):

```python
class CenterPotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/center"

class FermiEnergyCmd(MenuCommand):
    params = (md_out_path, fermi_unit)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/fermi"

class ElectrodePotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/electrode"

class PhiZProfileCmd(MenuCommand):
    params = (cube_pattern, max_curves)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/phi_z"

class ThicknessSensitivityCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit, thickness_end)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/thickness_sensitivity"

class FullPotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              max_curves, thickness_end)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential"
```

What was `if code in ("211", "213", "215", "216"):` is now each class listing
exactly what it needs.

### Charge (`_charge.py`)

One class, 3 instances with different `method`:

```python
class SurfaceChargeCmd(MenuCommand):
    advanced_params = (root_dir, dir_pattern, normal_axis, metal_elements,
                       layer_tol, n_surface_layers, outdir, frame_slice)

    def __init__(self, code: str, label: str, *, method: str | None = None):
        super().__init__(code, label)
        if method is not None:
            self.params = (FixedParam(K.METHOD, method),)
        else:
            self.params = (charge_method,)
        # output_subdir set dynamically based on method
        ...

    def execute(self, ctx):
        analyze = lazy_import("md_analysis.main", "run_charge_analysis")
        results = analyze(
            output_dir=ctx[K.OUTDIR_RESOLVED], method=ctx[K.METHOD], ...
        )
        # print results + ensemble summary
```

Tree assembly:
```python
SurfaceChargeCmd("221", "Surface Charge (Counterion)", method="counterion")
SurfaceChargeCmd("222", "Surface Charge (Layer)", method="layer")
SurfaceChargeCmd("223", "Full Charge Analysis with Plots", method=None)  # prompts
```

### Scripts (`_scripts.py`)

Non-standard param flow -> override `_collect_all_params()`:

```python
class BaderSingleCmd(MenuCommand):
    output_subdir = ""   # handled in execute()

    def _collect_all_params(self) -> dict:
        """Custom flow: show trajectory info between prompts."""
        ctx: dict = {}
        ctx[K.XYZ] = prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
        _print_trajectory_info(ctx[K.XYZ])
        ctx[K.CELL_ABC] = CellAbcParam().collect_and_return()  # or inline
        ctx[K.FRAME] = prompt_int("Frame number (0-based)", default=0) or 0
        ctx[K.OUTDIR] = prompt_str("Output directory", default=".") or "."
        ctx[K.WORKDIR_NAME] = prompt_str("Work directory name", default="bader") or "bader"
        ctx[K.SCRIPT_PATH] = _resolve_script_path()
        ctx[K.GEN_POTCAR] = prompt_bool("Generate POTCAR via vaspkit?", default=True)
        return ctx

    def execute(self, ctx):
        # read frame, call generate_bader_workdir
        ...
```

Overriding `_collect_all_params()` gives full control over the prompt sequence while
still inheriting `run()`'s error handling and other base behavior.

### Settings (`_settings.py`)

No params, each command overrides `execute()` directly:

```python
class SetAnalysisDefaultCmd(MenuCommand):
    """Generic command for setting one configurable analysis default.
    Instantiated 4 times with different config_key."""

    def __init__(self, code: str, label: str, *, config_key: str):
        super().__init__(code, label)
        self.config_key = config_key

    def execute(self, ctx):
        entry = CONFIGURABLE_DEFAULTS[self.config_key]
        current = get_config(self.config_key, entry["default"])
        value = prompt_float(entry["label"], default=current)
        if value <= 0:
            print("  Value must be > 0, skipping.")
            return
        set_config(self.config_key, value)
        print(f"  Saved: {self.config_key} = {value}")
```

Tree assembly (903-906 are 1 class x 4 instances):
```python
SetAnalysisDefaultCmd("903", "Layer Clustering Tolerance (A)", config_key=KEY_LAYER_TOL_A)
SetAnalysisDefaultCmd("904", "Z-axis Bin Width (A)", config_key=KEY_Z_BIN_WIDTH_A)
SetAnalysisDefaultCmd("905", "Theta Bin Width (deg)", config_key=KEY_THETA_BIN_DEG)
SetAnalysisDefaultCmd("906", "Water O-H Cutoff (A)", config_key=KEY_WATER_OH_CUTOFF_A)
```

---

## Tree Assembly (`__init__.py`)

```python
_BANNER = f"""\
================================================================================
            MD-Analysis: Metal-Water Interface Analysis Toolkit
                              Version {__version__}
================================================================================"""


def build_menu_tree() -> MenuGroup:
    root = MenuGroup("0", "MD-Analysis")

    # --- Water ---
    water = MenuGroup("1", "Water Analysis")
    water.add(
        WaterDensityCmd("101", "Water Mass Density Profile"),
        WaterOrientationCmd("102", "Water Orientation-Weighted Density Profile"),
        AdWaterOrientationCmd("103", "Adsorbed-Water Orientation Profile"),
        AdWaterThetaCmd("104", "Adsorbed-Water Theta Distribution"),
        WaterThreePanelCmd("105", "Full Water Three-Panel Analysis  (includes 101-104)"),
    )

    # --- Electrochemical ---
    electrochemical = MenuGroup("2", "Electrochemical Analysis")

    potential = MenuGroup("21", "Potential Analysis")
    potential.add(
        CenterPotentialCmd("211", "Center Slab Potential (phi_center)"),
        FermiEnergyCmd("212", "Fermi Energy Time Series"),
        ElectrodePotentialCmd("213", "Electrode Potential (U vs SHE)  (includes 211+212)"),
        PhiZProfileCmd("214", "Phi(z) Plane-Averaged Profile"),
        ThicknessSensitivityCmd("215", "Thickness Sensitivity Sweep"),
        FullPotentialCmd("216", "Full Potential Analysis  (includes 211-215)"),
    )

    charge = MenuGroup("22", "Charge Analysis")
    charge.add(
        SurfaceChargeCmd("221", "Surface Charge (Counterion)", method="counterion"),
        SurfaceChargeCmd("222", "Surface Charge (Layer)", method="layer"),
        SurfaceChargeCmd("223", "Full Charge Analysis with Plots", method=None),
    )

    electrochemical.add(potential, charge)

    # --- Enhanced Sampling ---
    enhanced = MenuGroup("3", "Enhanced Sampling")
    enhanced.add(
        SGQuickPlotCmd("301", "Slow-Growth Quick Plot"),
        SGPublicationPlotCmd("302", "Slow-Growth Publication Plot"),
    )

    # --- Scripts ---
    scripts = MenuGroup("4", "Scripts / Tools")
    scripts.add(
        BaderSingleCmd("401", "Generate Bader Work Directory (single frame)"),
        BaderBatchCmd("402", "Batch Generate Bader Work Directories"),
    )

    # --- Settings ---
    settings = MenuGroup("9", "Settings")
    settings.add(
        SetVaspScriptCmd("901", "Set VASP Submission Script Path"),
        ShowConfigCmd("902", "Show Current Configuration"),
        "Analysis Defaults",
        SetAnalysisDefaultCmd("903", "Layer Clustering Tolerance (A)",
                              config_key=KEY_LAYER_TOL_A),
        SetAnalysisDefaultCmd("904", "Z-axis Bin Width (A)",
                              config_key=KEY_Z_BIN_WIDTH_A),
        SetAnalysisDefaultCmd("905", "Theta Bin Width (deg)",
                              config_key=KEY_THETA_BIN_DEG),
        SetAnalysisDefaultCmd("906", "Water O-H Cutoff (A)",
                              config_key=KEY_WATER_OH_CUTOFF_A),
        ResetDefaultsCmd("907", "Reset All Defaults"),
    )

    root.add(water, electrochemical, enhanced, scripts, "", settings)
    root.build_flat_index()
    return root


def main() -> int:
    """Console script entry point."""
    import logging

    _logger = logging.getLogger("md_analysis")
    if not any(not isinstance(h, logging.NullHandler) for h in _logger.handlers):
        _handler = logging.StreamHandler()
        _handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        _logger.addHandler(_handler)
        _logger.setLevel(logging.INFO)

    try:
        print(_BANNER)
        tree = build_menu_tree()
        tree.run()          # loops until "0"
        print("\n Bye!")
        return 0
    except (KeyboardInterrupt, EOFError):
        print("\n\n Bye!")
        return 0
```

Control flow with loops:

```
main()
  print(banner)
  root.run()                         # LOOP
    print(root menu)
    "1" -> water.run()               # LOOP
             print(water menu)
             "101" -> WaterDensityCmd.run()
                        collect params -> execute -> print result
                      (returns to water loop, menu re-displayed)
             "0"   -> break (returns to root loop)
    "101" -> WaterDensityCmd.run()   # flat index shortcut
               (returns to root loop)
    "0"   -> break
  print("Bye!")
  return 0
```

---

## Test Strategy

### What changes

| Old path | New path / approach |
|---|---|
| `test_handle_cmd_error.py` testing `_handle_cmd_error` decorator | Test `MenuCommand.run()` error handling directly |
| `test_settings_defaults.py` testing `_cmd_903()` etc. | Test `SetAnalysisDefaultCmd.execute()` etc. |
| `test_logging_setup.py` referencing `cli._prompt` logger | Reference new logger location |
| `monkeypatch("builtins.input", ...)` | Use `set_input_source()` from `_prompt.py` |

### What stays the same

- Integration tests (`test/integration/`) don't import CLI internals -> no changes
- All analysis module tests are unaffected
- `pyproject.toml` entry point `md_analysis.cli:main` unchanged

### New test opportunities

- `MenuGroup._render_menu()` can be tested in isolation (pure string output)
- `ParamCollector.collect()` / `apply_default()` testable with `set_input_source()`
- `build_menu_tree()` can be tested for structural correctness (all codes unique,
  flat index complete, no orphan nodes)

---

## Design Decisions

### Why `run()` returns `None`, not `int`

In the old design, `run() -> int` was ambiguous: exit code for `MenuGroup`,
success/failure for `MenuCommand`, but parent `MenuGroup` passed it through
via `return child.run()`. With the loop design, `MenuGroup` **must not** exit on
child failure — it continues the loop. Making the return type `None` eliminates
the temptation to propagate error codes. `main()` always returns 0; non-zero exit
only happens if we add explicit `sys.exit()` calls in the future.

### Why tuples for `params` / `advanced_params`

Class-level mutable lists are a classic Python hazard — all instances share the
same list object, and accidental mutation (e.g. `self.params.append(...)`) affects
every instance. Tuples are immutable. Per-instance customization works by setting
`self.params` in `__init__()`, which shadows the class attribute safely.

### Why flat index on root only

Sub-menu flat lookup (e.g. typing "211" inside the water menu) would be confusing —
the user expects water codes (101-105) in the water menu. Cross-menu shortcuts
make sense only at the top level, matching VASPKIT behavior.

### Why `lazy_import()` instead of inline `from .. import`

1. Full dotted module path is explicit and greppable
2. Consistent pattern across all commands (easy to audit)
3. Relative imports in `execute()` methods would break if file moves
4. Optional: could be extended to cache imports in the future

### Why `set_input_source()` instead of patching `builtins.input`

`monkeypatch("builtins.input", ...)` works but is fragile — it affects ALL code in
the process, including third-party libraries. `set_input_source()` scopes the
override to CLI prompt functions only. Tests can also compose input sequences:

```python
def scripted_input(responses):
    it = iter(responses)
    return lambda _: next(it)

set_input_source(scripted_input(["md-pos-1.xyz", ".restart", "cp2k.restart", ...]))
```

### What was NOT adopted

- **State pattern**: The menu nodes don't change behavior based on state. Each
  node's display and dispatch logic is fixed. "Current position in the tree" is
  just the call stack — no explicit state object needed.
- **TypedDict for ctx**: Added complexity for marginal benefit. The `K` class
  constants provide key safety. If type checking becomes important later, a
  TypedDict can be introduced without changing the framework.
- **Auto-mapping ctx -> function kwargs**: Each command's `execute()` explicitly
  maps `ctx[K.XYZ]` to `xyz_path=`. This is verbose but transparent — easy to
  debug, easy to read, no hidden magic.

---

## Clarifications (Potential Ambiguities)

### 1. `output_subdir` 的两种模式

命令存在两种输出目录模式，由是否调用 `main.py` 中的 `run_*_analysis()` 决定：

| 模式 | output_subdir | execute() 传什么 | 谁创建子目录 |
|---|---|---|---|
| **直接调用模块函数** (101-104, 211-215) | `"water"`, `"potential/center"` 等 | `ctx[K.OUTDIR_RESOLVED]` (已含子路径) | 基类 `run()` |
| **调用 `run_*_analysis` 包装** (105, 216, 223) | `""` (空) | `Path(ctx[K.OUTDIR])` (原始路径) | `run_*_analysis()` 内部 |

第二种模式下 `run_*_analysis()` 自己创建 `output_dir/water/`、`output_dir/charge/<method>/`
等子目录。基类 `run()` 检测到 `output_subdir` 为空时跳过 outdir 创建，也不设置
`ctx[K.OUTDIR_RESOLVED]`。

**Charge 的动态子目录**：`SurfaceChargeCmd` 使用模式二（`output_subdir = ""`），调用
`run_charge_analysis(output_dir=Path(ctx[K.OUTDIR]), method=ctx[K.METHOD], ...)`，由
`run_charge_analysis` 内部创建 `output_dir/charge/<method>/`。不需要在基类层面处理动态
`output_subdir`。

### 2. Scripts 命令中 CellAbcParam 的使用方式

DESIGN.md 中 `BaderSingleCmd._collect_all_params()` 写了 `CellAbcParam().collect_and_return()`，
但 `ParamCollector` 的 API 只有 `collect(ctx)` 和 `apply_default(ctx)`，没有
`collect_and_return()`。

正确写法应为：

```python
def _collect_all_params(self) -> dict:
    ctx: dict = {}
    ctx[K.XYZ] = prompt_str_required("XYZ trajectory file")
    _print_trajectory_info(ctx[K.XYZ])
    cell_abc.collect(ctx)          # CellAbcParam 实例，写入 ctx[K.CELL_ABC]
    ctx[K.FRAME] = prompt_int("Frame number (0-based)", default=0) or 0
    ...
    return ctx
```

即直接复用 `_params.py` 中的预定义 `cell_abc` 实例，调用 `collect(ctx)` 即可。
不引入额外的 API。

### 3. `_read()` 的作用范围

`_prompt.py` 中的 `_read()` 函数（封装 `_input_fn`）不仅被所有 `prompt_*` 函数使用，
**`MenuGroup.run()` 中的 `input(" Input: ")` 也必须改为 `_read(" Input: ")`**。

这样 `set_input_source()` 才能覆盖所有用户输入点，包括菜单选择和参数收集。实现时
`_framework.py` 从 `_prompt.py` 导入 `_read`：

```python
# _framework.py
from ._prompt import _read

class MenuGroup(MenuNode):
    def run(self) -> None:
        while True:
            print(self._render_menu())
            choice = _read(" Input: ").strip()
            ...
```

### 4. `FixedParam` 的行为

`FixedParam` 用于 `SurfaceChargeCmd(method="counterion")` 等场景，**collect() 和
apply_default() 行为相同**——都是无条件写入固定值，不提示用户：

```python
class FixedParam(ParamCollector):
    def __init__(self, key: str, value: Any):
        self.key = key
        self.value = value

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = self.value

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.value
```

这意味着当 `FixedParam` 出现在 `params` 元组中时，用户不会被提示输入——值在命令实例化
时就已确定。与 `ChoiceParam`（提示用户选择）形成对比。
