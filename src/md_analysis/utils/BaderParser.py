"""Parse VASP Bader charge analysis output and attach results to ASE Atoms."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read as ase_read


from ..exceptions import MDAnalysisError


class BaderParseError(MDAnalysisError):
    """Raised when ACF file format is invalid or atom counts mismatch."""


def _read_acf(path: Path) -> tuple[np.ndarray, dict]:
    """Parse an ACF.dat file produced by the Bader charge analysis program.

    Parameters
    ----------
    path : Path
        Path to the ACF.dat file.

    Returns
    -------
    data : np.ndarray, shape (natoms, 4)
        Columns: x, y, z, charge.
    footer : dict
        Keys: ``vacuum_charge``, ``vacuum_volume``, ``number_of_electrons``.
    """
    lines = Path(path).read_text().splitlines()

    data_rows: list[list[float]] = []
    footer: dict[str, float] = {}

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("-"):
            continue

        upper = stripped.upper()
        if upper.startswith("VACUUM CHARGE"):
            footer["vacuum_charge"] = float(stripped.split(":")[-1])
        elif upper.startswith("VACUUM VOLUME"):
            footer["vacuum_volume"] = float(stripped.split(":")[-1])
        elif upper.startswith("NUMBER OF ELECTRONS"):
            footer["number_of_electrons"] = float(stripped.split(":")[-1])
        else:
            parts = stripped.split()
            if len(parts) < 5:
                continue
            try:
                # columns: index, x, y, z, charge, min_dist, atomic_vol
                _idx = int(parts[0])
                x, y, z, charge = (float(p) for p in parts[1:5])
                data_rows.append([x, y, z, charge])
            except ValueError:
                continue

    if not data_rows:
        raise BaderParseError(f"No atom data found in {path}")

    return np.array(data_rows), footer


def _read_potcar_zval(path: Path) -> list[tuple[str, float]]:
    """Extract (element_symbol, zval) pairs from a VASP POTCAR file.

    Parameters
    ----------
    path : Path
        Path to the POTCAR file.

    Returns
    -------
    list of (str, float)
        Element symbol and valence electron count for each species block,
        in the order they appear in POTCAR (matching POSCAR element order).
    """
    lines = Path(path).read_text().splitlines()
    results: list[tuple[str, float]] = []

    element: str | None = None
    for line in lines:
        stripped = line.strip()
        # The first line of each element block starts with PAW_PBE/PAW_LDA/US:
        #   PAW_PBE Cu_pv 06Sep2000
        # Avoid matching TITEL lines (e.g. "TITEL  = PAW_PBE Cu ...").
        if stripped.startswith(("PAW_PBE", "PAW_LDA", "PAW_GGA", "US ")):
            token = stripped.split()[1]
            # Strip suffixes like _pv, _sv, _GW, etc.
            element = token.split("_")[0]

        if "ZVAL" in line and element is not None:
            # Line format: "   POMASS =   63.546; ZVAL   =   11.000    mass and valenz"
            after_zval = line.split("ZVAL")[1]
            # Extract the number after '='
            num_str = after_zval.split("=")[1].split()[0]
            zval = float(num_str)
            results.append((element, zval))
            element = None  # reset for next block

    if not results:
        raise BaderParseError(f"No ZVAL entries found in {path}")

    return results


def load_bader_atoms(
    structure_path: Path,
    acf_path: Path,
    potcar_path: Path,
) -> Atoms:
    """Load structure and Bader charges, returning an augmented Atoms object.

    Parameters
    ----------
    structure_path : Path
        Path to POSCAR / CONTCAR (or any ASE-readable structure file).
    acf_path : Path
        Path to ACF.dat from Bader analysis.
    potcar_path : Path
        Path to POTCAR used in the calculation.

    Returns
    -------
    atoms : ase.Atoms
        The structure with two extra per-atom arrays:
        - ``"bader_charge"`` : raw electron count from ACF.dat
        - ``"bader_net_charge"`` : ZVAL - bader_charge (positive = lost electrons)
    """
    atoms = ase_read(str(structure_path))
    data, _footer = _read_acf(acf_path)

    if len(data) != len(atoms):
        raise BaderParseError(
            f"Atom count mismatch: ACF.dat has {len(data)} atoms, "
            f"structure has {len(atoms)} atoms"
        )

    zval_list = _read_potcar_zval(potcar_path)

    # Build per-atom ZVAL array using POSCAR element order + counts
    symbols = atoms.get_chemical_symbols()
    # Map element -> zval from POTCAR (preserve order for lookup)
    zval_map: dict[str, float] = {}
    for elem, zval in zval_list:
        zval_map[elem] = zval

    zval_per_atom = np.empty(len(atoms))
    for i, sym in enumerate(symbols):
        if sym not in zval_map:
            raise BaderParseError(
                f"Element '{sym}' in structure not found in POTCAR"
            )
        zval_per_atom[i] = zval_map[sym]

    bader_charge = data[:, 3]
    bader_net_charge = zval_per_atom - bader_charge

    atoms.arrays["bader_charge"] = bader_charge
    atoms.arrays["bader_net_charge"] = bader_net_charge

    return atoms
