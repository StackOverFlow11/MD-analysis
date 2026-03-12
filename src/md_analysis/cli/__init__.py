"""VASPKIT-style interactive CLI for md-analysis.

Entry point registered as ``md-analysis`` console script.
"""

from __future__ import annotations

from .. import __version__
from ..config import (
    KEY_LAYER_TOL_A,
    KEY_THETA_BIN_DEG,
    KEY_WATER_OH_CUTOFF_A,
    KEY_Z_BIN_WIDTH_A,
)
from ._charge import SurfaceChargeCmd
from ._enhanced_sampling import SGPublicationPlotCmd, SGQuickPlotCmd
from ._framework import MenuGroup
from ._potential import (
    CenterPotentialCmd,
    ElectrodePotentialCmd,
    FermiEnergyCmd,
    FullPotentialCmd,
    PhiZProfileCmd,
    ThicknessSensitivityCmd,
)
from ._scripts import BaderBatchCmd, BaderSingleCmd
from ._settings import (
    ResetDefaultsCmd,
    SetAnalysisDefaultCmd,
    SetVaspScriptCmd,
    ShowConfigCmd,
)
from ._water import (
    AdWaterOrientationCmd,
    AdWaterThetaCmd,
    WaterDensityCmd,
    WaterOrientationCmd,
    WaterThreePanelCmd,
)

_BANNER = f"""\
================================================================================
            MD-Analysis: Metal-Water Interface Analysis Toolkit
                              Version {__version__}
================================================================================"""


def build_menu_tree() -> MenuGroup:
    """Assemble the full CLI menu tree."""
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
        CenterPotentialCmd("201", "Center Slab Potential (phi_center)"),
        FermiEnergyCmd("202", "Fermi Energy Time Series"),
        ElectrodePotentialCmd("203", "Electrode Potential (U vs SHE)  (includes 201+202)"),
        PhiZProfileCmd("204", "Phi(z) Plane-Averaged Profile"),
        ThicknessSensitivityCmd("205", "Thickness Sensitivity Sweep"),
        FullPotentialCmd("206", "Full Potential Analysis  (includes 201-205)"),
    )

    charge = MenuGroup("22", "Charge Analysis")
    charge.add(
        SurfaceChargeCmd("301", "Surface Charge (Counterion)", method="counterion"),
        SurfaceChargeCmd("302", "Surface Charge (Layer)", method="layer"),
        SurfaceChargeCmd("303", "Full Charge Analysis with Plots", method=None),
    )

    electrochemical.add(potential, charge)

    # --- Enhanced Sampling ---
    enhanced = MenuGroup("3", "Enhanced Sampling")
    enhanced.add(
        SGQuickPlotCmd("501", "Slow-Growth Quick Plot"),
        SGPublicationPlotCmd("502", "Slow-Growth Publication Plot"),
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
    if not _logger.handlers or all(
        isinstance(h, logging.NullHandler) for h in _logger.handlers
    ):
        _handler = logging.StreamHandler()
        _handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        _logger.addHandler(_handler)
        _logger.setLevel(logging.INFO)

    try:
        print(_BANNER)
        tree = build_menu_tree()
        tree.run()
        print("\n Bye!")
        return 0
    except (KeyboardInterrupt, EOFError):
        print("\n\n Bye!")
        return 0
