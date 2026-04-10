"""Parameter collection hierarchy and predefined instances."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any

from ..config import (
    KEY_CP2K_SCRIPT_PATH,
    KEY_DP_SP_INP_TEMPLATE_PATH,
    KEY_LAYER_TOL_A,
    KEY_SP_INP_TEMPLATE_PATH,
    KEY_VASP_SCRIPT_PATH,
    KEY_Z_BIN_WIDTH_A,
)
from ._prompt import prompt_bool, prompt_choice, prompt_float, prompt_int, prompt_str, prompt_str_required


class K:
    """Parameter dict key constants. Prevents string-key typos."""

    XYZ = "xyz"
    CELL_ABC = "cell_abc"
    DZ_A = "dz_A"
    LAYER_TOL = "layer_tol"
    OUTDIR = "outdir"
    # Trajectory frame selection (shared by Bader/Potential/SpGen batch commands)
    FRAME_MODE = "frame_mode"
    TIME_START_FS = "time_start_fs"
    TIME_END_FS = "time_end_fs"
    TIME_STEP_FS = "time_step_fs"
    SINGLE_TIME_FS = "single_time_fs"
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
    TARGET_SIDE = "target_side"
    # scripts-specific
    FRAME = "frame"
    WORKDIR_NAME = "workdir_name"
    SCRIPT_PATH = "script_path"
    GEN_POTCAR = "gen_potcar"
    INP_TEMPLATE = "inp_template"
    # enhanced-sampling-specific
    RESTART_PATH = "restart_path"
    LOG_PATH = "log_path"
    INITIAL_STEP = "initial_step"
    FINAL_STEP = "final_step"
    COLVAR_ID = "colvar_id"
    # Atom charge tracking
    ATOM_INDICES_XYZ = "atom_indices_xyz"
    # TI-specific
    INP_PATH = "inp_path"
    TARGET_AU = "target_au"
    TARGETS_AU = "targets_au"
    TIME_INITIAL_FS = "time_initial_fs"
    TIME_FINAL_FS = "time_final_fs"
    N_POINTS = "n_points"
    STEPS = "steps"
    # Calibration-specific
    CALIBRATION_CSV = "calibration_csv"
    FITTING_METHOD = "fitting_method"
    POLY_DEGREE = "poly_degree"
    SIGMA_VALUE = "sigma_value"
    CALIBRATION_JSON = "calibration_json"
    POTENTIAL_REFERENCE = "potential_reference"
    TEMPERATURE_K = "temperature_k"
    PH = "ph"
    PHI_PZC = "phi_pzc"
    # Constrained-TI analysis
    EQUILIBRATION = "equilibration"
    SEM_TARGET = "sem_target"
    TI_ROOT_DIR = "ti_root_dir"
    TI_DIR_PATTERN = "ti_dir_pattern"  # distinct from DIR_PATTERN (Bader)
    EPSILON_TOL_EV = "epsilon_tol_ev"
    TI_REVERSE = "ti_reverse"
    AUTO_EQUILIBRATION = "auto_equilibration"
    # Potential input mode (continuous MD vs distributed SP)
    INPUT_MODE = "input_mode"
    SP_ROOT_DIR = "sp_root_dir"
    SP_DIR_PATTERN = "sp_dir_pattern"
    SP_CUBE_FILENAME = "sp_cube_filename"
    SP_OUT_FILENAME = "sp_out_filename"


# ---------------------------------------------------------------------------
# ParamCollector ABC
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Generic param classes (~80% of params)
# ---------------------------------------------------------------------------

class StrParam(ParamCollector):
    """Single string prompt with optional default."""

    def __init__(self, key: str, label: str, *, default: str | None = None,
                 required: bool = False):
        self.key = key
        self.label = label
        self.default = default
        self.required = required

    def collect(self, ctx: dict) -> None:
        if self.required:
            ctx[self.key] = prompt_str_required(self.label)
        else:
            ctx[self.key] = prompt_str(self.label, default=self.default) or self.default

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.default


class FloatParam(ParamCollector):
    def __init__(self, key: str, label: str, *, default: float):
        self.key = key
        self.label = label
        self.default = default

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_float(self.label, default=self.default)

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.default


class IntParam(ParamCollector):
    def __init__(self, key: str, label: str, *, default: int | None = None):
        self.key = key
        self.label = label
        self.default = default

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_int(self.label, default=self.default)

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.default


class ChoiceParam(ParamCollector):
    def __init__(self, key: str, label: str, choices: list[str], *,
                 default: str):
        self.key = key
        self.label = label
        self.choices = choices
        self.default = default

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_choice(self.label, self.choices, default=self.default)

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.default


class FixedParam(ParamCollector):
    """Injects a fixed value without prompting."""

    def __init__(self, key: str, value: Any):
        self.key = key
        self.value = value

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = self.value

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.value


class ConditionalParam(ParamCollector):
    """Wraps another ParamCollector; prompts only when ``predicate(ctx)`` is true.

    When the predicate is false, the inner collector's default value is applied
    silently. Requires any ctx keys the predicate reads to be collected earlier
    in the same ``params`` tuple (predicate runs at ``collect`` time).
    """

    def __init__(self, inner: ParamCollector, predicate):
        self.inner = inner
        self.predicate = predicate

    def collect(self, ctx: dict) -> None:
        if self.predicate(ctx):
            self.inner.collect(ctx)
        else:
            self.inner.apply_default(ctx)

    def apply_default(self, ctx: dict) -> None:
        self.inner.apply_default(ctx)


class BoolParam(ParamCollector):
    """Yes/no prompt stored as bool."""

    def __init__(self, key: str, label: str, *, default: bool = True) -> None:
        self.key = key
        self.label = label
        self.default = default

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_bool(self.label, default=self.default)

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self.default


class DisplayAction(ParamCollector):
    """Execute a side-effect (e.g. print info) during parameter collection.

    Stores no value in ``ctx``.  Skipped silently when advanced params are
    not prompted (``apply_default`` is a no-op).

    The *action* callable receives the current ``ctx`` dict and must not
    rely on keys that have not yet been collected (ordering matters).
    """

    def __init__(self, action: Callable[[dict], None]) -> None:
        self.action = action

    def collect(self, ctx: dict) -> None:
        self.action(ctx)

    def apply_default(self, ctx: dict) -> None:
        pass  # intentional no-op: display actions are skipped in non-interactive mode


class _ConfigBackedParam(ParamCollector):
    """Internal base for params whose default comes from persistent user config.

    Subclasses must implement ``collect()`` with the appropriate prompt function.
    Note: ``apply_default`` may return ``None`` when no config value is set;
    subclasses may override to provide a registry-backed fallback.
    """

    def __init__(self, key: str, label: str, *, config_key: str) -> None:
        self.key = key
        self.label = label
        self.config_key = config_key

    def _get_config_value(self) -> Any:
        """Read current value from user config; returns None if unset."""
        from ..config import get_config

        return get_config(self.config_key)

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self._get_config_value()


class ConfigDefaultParam(_ConfigBackedParam):
    """Float prompt with default from persistent config registry.

    Uses ``CONFIGURABLE_DEFAULTS[config_key]["default"]`` as fallback
    when the user has not set a custom value.
    """

    def _get_effective_default(self) -> float:
        from ..config import CONFIGURABLE_DEFAULTS, get_config

        entry = CONFIGURABLE_DEFAULTS[self.config_key]
        return get_config(self.config_key, entry["default"])

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_float(self.label, default=self._get_effective_default())

    def apply_default(self, ctx: dict) -> None:
        ctx[self.key] = self._get_effective_default()


class ConfigStrParam(_ConfigBackedParam):
    """String prompt with default from persistent user config.

    Unlike ``ConfigDefaultParam``, does not require a ``CONFIGURABLE_DEFAULTS``
    registry entry.  Returns ``None`` when no config value is set.
    """

    def collect(self, ctx: dict) -> None:
        ctx[self.key] = prompt_str(self.label, default=self._get_config_value())


# ---------------------------------------------------------------------------
# Special param classes (~20% of params)
# ---------------------------------------------------------------------------

class CellAbcParam(ParamCollector):
    """Cell parameter acquisition with .restart/md.inp selection + retry."""

    def collect(self, ctx: dict) -> None:
        from ..utils.RestartParser.CellParser import (
            CellParseError,
            parse_abc_from_md_inp,
            parse_abc_from_restart,
        )

        for attempt in range(2):
            source = prompt_choice("Cell source", [".restart", "md.inp"],
                                   default=".restart")
            try:
                if source == ".restart":
                    path = prompt_str_required("CP2K .restart file")
                    abc = parse_abc_from_restart(path)
                else:
                    path = prompt_str_required("CP2K input file (e.g. md.inp)")
                    abc = parse_abc_from_md_inp(path)
            except (CellParseError, FileNotFoundError) as exc:
                print(f"\n  Error: {exc}")
                if attempt == 0:
                    if prompt_bool("Retry with different file?", default=True):
                        continue
                raise CellParseError(str(exc)) from exc

            print(f"  Cell: a={abc[0]:.4f}, b={abc[1]:.4f}, c={abc[2]:.4f} A")
            ctx[K.CELL_ABC] = abc
            return

        raise CellParseError("Failed to parse cell parameters.")  # pragma: no cover

    def apply_default(self, ctx: dict) -> None:
        pass  # CellAbc has no meaningful default


class AtomIndicesParam(ParamCollector):
    """Comma-separated integer atom indices (XYZ 0-based)."""

    def collect(self, ctx: dict) -> None:
        raw = prompt_str_required("Atom indices (XYZ, 0-based, comma-separated, e.g. 0,5,10)")
        indices = [int(x.strip()) for x in raw.split(",") if x.strip()]
        if not indices:
            raise ValueError("At least one atom index is required")
        ctx[K.ATOM_INDICES_XYZ] = indices

    def apply_default(self, ctx: dict) -> None:
        ctx[K.ATOM_INDICES_XYZ] = []


class MetalElementsParam(ParamCollector):
    """Comma-separated input -> set[str] | None."""

    def collect(self, ctx: dict) -> None:
        raw = prompt_str("Metal elements (comma-separated, e.g. Cu,Ag)", default=None)
        if raw is None:
            ctx[K.METAL_ELEMENTS] = None
        else:
            elements = {e.strip() for e in raw.split(",") if e.strip()}
            ctx[K.METAL_ELEMENTS] = elements or None

    def apply_default(self, ctx: dict) -> None:
        ctx[K.METAL_ELEMENTS] = None


class FrameSliceParam(ParamCollector):
    """Collects frame_start, frame_end, frame_step as a group."""

    def collect(self, ctx: dict) -> None:
        ctx[K.FRAME_START] = prompt_int("Frame start (0-based)", default=None)
        ctx[K.FRAME_END] = prompt_int("Frame end (exclusive)", default=None)
        ctx[K.FRAME_STEP] = prompt_int("Frame step", default=None)

    def apply_default(self, ctx: dict) -> None:
        ctx[K.FRAME_START] = None
        ctx[K.FRAME_END] = None
        ctx[K.FRAME_STEP] = None


# ---------------------------------------------------------------------------
# Predefined instances
# ---------------------------------------------------------------------------

# Required params (always in `params`, never in `advanced_params`)
xyz_path = StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True)
xyz_path_opt = StrParam(K.XYZ, "XYZ trajectory file", default="md-pos-1.xyz")
cell_abc = CellAbcParam()
cube_pattern = StrParam(K.CUBE_PATTERN, "Cube file glob pattern",
                        default="md-POTENTIAL-v_hartree-1_*.cube")
md_out_path = StrParam(K.MD_OUT, "CP2K md.out file", default="md.out")

# Typed params
thickness = FloatParam(K.THICKNESS, "Slab averaging thickness (A)", default=7.5)
center_mode = ChoiceParam(K.CENTER_MODE, "Slab center mode",
                          ["interface", "cell"], default="interface")
fermi_unit = ChoiceParam(K.FERMI_UNIT, "Fermi energy unit in md.out",
                         ["au", "ev"], default="au")
max_curves = IntParam(K.MAX_CURVES, "Max curves on overlay (0=all)", default=0)
thickness_end = FloatParam(K.THICKNESS_END, "Thickness sweep upper limit (A)",
                           default=15.0)
normal_axis = ChoiceParam(K.NORMAL, "Surface normal axis",
                          ["a", "b", "c"], default="c")
charge_method = ChoiceParam(K.METHOD, "Charge method",
                            ["counterion", "layer"], default="counterion")
n_surface_layers = IntParam(K.N_SURFACE_LAYERS,
                            "Number of surface layers per interface", default=1)
target_side = ChoiceParam(K.TARGET_SIDE, "Target surface side",
                          ["aligned", "opposed"], default="aligned")

# Config-backed params
dz_bin = ConfigDefaultParam(K.DZ_A, "Z-axis bin width (A)",
                            config_key=KEY_Z_BIN_WIDTH_A)
layer_tol = ConfigDefaultParam(K.LAYER_TOL, "Layer clustering tolerance (A)",
                               config_key=KEY_LAYER_TOL_A)

# Structural params
metal_elements = MetalElementsParam()
root_dir = StrParam(K.ROOT_DIR, "Root directory", default=".")
dir_pattern = StrParam(K.DIR_PATTERN, "Frame subdirectory pattern",
                       default="bader_t*_i*")
outdir = StrParam(K.OUTDIR, "Output directory", default="analysis")
frame_slice = FrameSliceParam()
atom_indices_xyz = AtomIndicesParam()

# Calibration params
calibration_csv_path = StrParam(K.CALIBRATION_CSV, "Calibration CSV file path")
fitting_method = ChoiceParam(K.FITTING_METHOD, "Fitting method",
                             ["linear", "polynomial", "spline",
                              "differential_capacitance"], default="linear")
poly_degree = IntParam(K.POLY_DEGREE, "Polynomial degree", default=2)
sigma_value = FloatParam(K.SIGMA_VALUE,
                         "Surface charge density sigma (uC/cm^2)", default=0.0)
calibration_json = StrParam(K.CALIBRATION_JSON, "Calibration JSON file path")

# Potential input mode params
input_mode = ChoiceParam(K.INPUT_MODE, "Input mode",
                         ["continuous", "distributed"], default="continuous")
_sp_root_dir_base = StrParam(K.SP_ROOT_DIR,
                             "Root directory of SP calculations", default=".")
_sp_dir_pattern_base = StrParam(K.SP_DIR_PATTERN, "SP directory pattern",
                                default="potential_t*_i*")
_sp_cube_filename_base = StrParam(K.SP_CUBE_FILENAME, "Cube filename in SP dirs",
                                  default="sp_potential-v_hartree-1_0.cube")
_sp_out_filename_base = StrParam(K.SP_OUT_FILENAME, "Output filename in SP dirs",
                                 default="sp.out")


def _is_distributed_mode(ctx: dict) -> bool:
    return ctx.get(K.INPUT_MODE) == "distributed"


# Only prompt sp_* params when input_mode == "distributed".
# Requires ``input_mode`` to appear earlier in the params tuple.
sp_root_dir = ConditionalParam(_sp_root_dir_base, _is_distributed_mode)
sp_dir_pattern = ConditionalParam(_sp_dir_pattern_base, _is_distributed_mode)
sp_cube_filename = ConditionalParam(_sp_cube_filename_base, _is_distributed_mode)
sp_out_filename = ConditionalParam(_sp_out_filename_base, _is_distributed_mode)
potential_reference = ChoiceParam(K.POTENTIAL_REFERENCE,
                                  "Output potential reference",
                                  ["SHE", "RHE", "PZC"], default="SHE")
temperature_k = FloatParam(K.TEMPERATURE_K, "Temperature (K)", default=298.15)
ph_value = FloatParam(K.PH, "pH", default=0.0)
phi_pzc = FloatParam(K.PHI_PZC, "Potential of zero charge (V vs SHE)", default=0.0)

# Config-backed script / template params (used by _scripts.py declarative commands)
vasp_script = ConfigStrParam(K.SCRIPT_PATH, "Submission script path",
                             config_key=KEY_VASP_SCRIPT_PATH)
cp2k_script = ConfigStrParam(K.SCRIPT_PATH, "Submission script path",
                             config_key=KEY_CP2K_SCRIPT_PATH)
sp_inp_template = ConfigStrParam(K.INP_TEMPLATE,
                                 "SP inp template path (e.g. sp.inp)",
                                 config_key=KEY_SP_INP_TEMPLATE_PATH)
dp_sp_inp_template = ConfigStrParam(K.INP_TEMPLATE,
                                    "DP SP inp template path (e.g. sp.inp)",
                                    config_key=KEY_DP_SP_INP_TEMPLATE_PATH)
gen_potcar = BoolParam(K.GEN_POTCAR, "Generate POTCAR via vaspkit?", default=True)

# Trajectory frame selection (shared by Bader/Potential/SpGen Batch + Single commands)
frame_mode = ChoiceParam(K.FRAME_MODE, "Frame selection mode",
                         choices=["index", "time"], default="time")
time_start_fs = FloatParam(K.TIME_START_FS, "Start time (fs)", default=0.0)
time_end_fs = FloatParam(K.TIME_END_FS, "End time (fs)", default=1000.0)
time_step_fs = FloatParam(K.TIME_STEP_FS, "Time step (fs)", default=10.0)
single_time_fs = FloatParam(K.SINGLE_TIME_FS, "Target time (fs)", default=0.0)
