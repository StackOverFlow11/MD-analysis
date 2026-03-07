"""VASPKIT-style interactive CLI for md-analysis.

Entry point registered as ``md-analysis`` console script.
"""

from __future__ import annotations

import sys

from .. import __version__


_BANNER = f"""\
================================================================================
            MD-Analysis: Metal-Water Interface Analysis Toolkit
                              Version {__version__}
================================================================================"""

_MENU = """\

 1) Water Analysis
 2) Potential Analysis
 3) Charge Analysis
 4) Scripts / Tools

 9) Settings
 0) Exit

================================================================================"""


def main() -> int:
    """Interactive top-level menu."""
    try:
        print(_BANNER)
        print(_MENU)
        choice = input(" Input: ").strip()

        if choice == "1":
            from ._water import water_menu
            return water_menu()
        elif choice == "2":
            from ._potential import potential_menu
            return potential_menu()
        elif choice == "3":
            from ._charge import charge_menu
            return charge_menu()
        elif choice == "4":
            from ._scripts import scripts_menu
            return scripts_menu()
        elif choice == "9":
            from ._settings import settings_menu
            return settings_menu()
        elif choice == "0":
            print("\n Bye!")
            return 0
        else:
            print(f"\n Invalid choice: {choice!r}")
            return 1

    except (KeyboardInterrupt, EOFError):
        print("\n\n Bye!")
        return 0
