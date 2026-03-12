"""Tests for md_analysis.utils.BaderParser."""

from pathlib import Path

import numpy as np
import pytest

from md_analysis.utils.BaderParser import (
    BaderParseError,
    _read_acf,
    _read_potcar_zval,
    load_bader_atoms,
)

DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "bader" / "bader_work_dir"


class TestReadAcf:
    """Tests for _read_acf."""

    def test_shape_and_values(self):
        data, footer = _read_acf(DATA_DIR / "ACF.dat")
        # 274 atoms based on POSCAR: 62 Cu + 2 Ag + 70 O + 140 H
        assert data.shape == (274, 4)
        # First atom: x=9.719248, y=0.071976, z=25.760657, charge=11.006087
        np.testing.assert_allclose(data[0, 0], 9.719248, atol=1e-5)
        np.testing.assert_allclose(data[0, 3], 11.006087, atol=1e-5)

    def test_footer(self):
        _data, footer = _read_acf(DATA_DIR / "ACF.dat")
        assert footer["vacuum_charge"] == pytest.approx(0.0)
        assert footer["number_of_electrons"] == pytest.approx(1264.0)

    def test_empty_file_raises(self, tmp_path):
        acf = tmp_path / "ACF.dat"
        acf.write_text("# header\n---\n")
        with pytest.raises(BaderParseError, match="No atom data"):
            _read_acf(acf)


class TestReadPotcarZval:
    """Tests for _read_potcar_zval."""

    def test_elements_and_zval(self):
        result = _read_potcar_zval(DATA_DIR / "POTCAR")
        assert len(result) == 4
        elements = [e for e, _ in result]
        zvals = [z for _, z in result]
        assert elements == ["Cu", "Ag", "O", "H"]
        assert zvals == pytest.approx([11.0, 11.0, 6.0, 1.0])

    def test_empty_potcar_raises(self, tmp_path):
        potcar = tmp_path / "POTCAR"
        potcar.write_text("some random content\n")
        with pytest.raises(BaderParseError, match="No ZVAL"):
            _read_potcar_zval(potcar)


class TestLoadBaderAtoms:
    """Tests for load_bader_atoms."""

    def test_full_pipeline(self):
        atoms = load_bader_atoms(
            DATA_DIR / "POSCAR", DATA_DIR / "ACF.dat", DATA_DIR / "POTCAR"
        )
        assert len(atoms) == 274
        assert "bader_charge" in atoms.arrays
        assert "bader_net_charge" in atoms.arrays
        assert atoms.arrays["bader_charge"].shape == (274,)
        assert atoms.arrays["bader_net_charge"].shape == (274,)

    def test_net_charge_values(self):
        atoms = load_bader_atoms(
            DATA_DIR / "POSCAR", DATA_DIR / "ACF.dat", DATA_DIR / "POTCAR"
        )
        # For Cu (ZVAL=11), net_charge = 11 - bader_charge
        # First atom is Cu with bader_charge ~11.006 → net_charge ~ -0.006
        net = atoms.arrays["bader_net_charge"]
        bc = atoms.arrays["bader_charge"]
        np.testing.assert_allclose(net[0], 11.0 - bc[0])
        # H atoms (last 140) have ZVAL=1
        np.testing.assert_allclose(net[-1], 1.0 - bc[-1])

    def test_charge_conservation(self):
        atoms = load_bader_atoms(
            DATA_DIR / "POSCAR", DATA_DIR / "ACF.dat", DATA_DIR / "POTCAR"
        )
        # Sum of bader_charge should equal NUMBER OF ELECTRONS (1264)
        total = atoms.arrays["bader_charge"].sum()
        assert total == pytest.approx(1264.0, abs=0.01)

    def test_atom_count_mismatch(self, tmp_path):
        # Create a fake ACF.dat with only 2 atoms
        acf = tmp_path / "ACF.dat"
        acf.write_text(
            "    #         X           Y           Z       CHARGE\n"
            " ----\n"
            "    1  0.0  0.0  0.0  1.0  0.5  1.0\n"
            "    2  1.0  1.0  1.0  1.0  0.5  1.0\n"
            " ----\n"
            "    VACUUM CHARGE:  0.0\n"
            "    VACUUM VOLUME:  0.0\n"
            "    NUMBER OF ELECTRONS:  2.0\n"
        )
        with pytest.raises(BaderParseError, match="Atom count mismatch"):
            load_bader_atoms(DATA_DIR / "POSCAR", acf, DATA_DIR / "POTCAR")
