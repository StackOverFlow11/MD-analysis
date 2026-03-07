"""Tests for md_analysis.scripts.BaderGen — Bader work directory generation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from md_analysis.scripts.BaderGen import (
    BaderGenError,
    DEFAULT_WORKDIR_NAME,
    generate_bader_workdir,
)
from md_analysis.scripts.utils.IndexMapper import read_index_map_from_poscar


def _make_test_atoms() -> Atoms:
    """Create a small Cu+O+H slab for testing."""
    symbols = ["Cu", "Cu", "O", "H", "H", "O", "H", "H"]
    positions = [
        [0.0, 0.0, 0.0],
        [1.8, 1.8, 0.0],
        [0.9, 0.9, 2.5],
        [0.9, 0.2, 3.1],
        [0.9, 1.6, 3.1],
        [2.7, 0.9, 2.5],
        [2.7, 0.2, 3.1],
        [2.7, 1.6, 3.1],
    ]
    atoms = Atoms(symbols=symbols, positions=positions,
                  cell=[3.6, 3.6, 10.0], pbc=True)
    return atoms


class TestGenerateBaderWorkdir:
    """Tests for generate_bader_workdir."""

    def test_creates_directory_structure(self, tmp_path):
        """Workdir is created with POSCAR, INCAR, KPOINTS."""
        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, generate_potcar=False,
        )
        assert workdir == tmp_path / DEFAULT_WORKDIR_NAME
        assert (workdir / "POSCAR").is_file()
        assert (workdir / "INCAR").is_file()
        assert (workdir / "KPOINTS").is_file()

    def test_poscar_has_index_map(self, tmp_path):
        """POSCAR first line contains a valid IndexMap encoding."""
        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, frame=5, source="test.xyz",
            generate_potcar=False,
        )
        imap = read_index_map_from_poscar(workdir / "POSCAR")
        assert imap.frame == 5
        assert imap.source == "test.xyz"
        assert imap.n_atoms == len(atoms)

    def test_incar_matches_template(self, tmp_path):
        """INCAR content matches the bundled template."""
        from importlib.resources import as_file, files

        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, generate_potcar=False,
        )
        template_pkg = files("md_analysis.scripts.template")
        with as_file(template_pkg / "INCAR") as tpl:
            expected = tpl.read_text()
        assert (workdir / "INCAR").read_text() == expected

    def test_kpoints_matches_template(self, tmp_path):
        """KPOINTS content matches the bundled template."""
        from importlib.resources import as_file, files

        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, generate_potcar=False,
        )
        template_pkg = files("md_analysis.scripts.template")
        with as_file(template_pkg / "KPOINTS") as tpl:
            expected = tpl.read_text()
        assert (workdir / "KPOINTS").read_text() == expected

    def test_script_copied(self, tmp_path):
        """script.sh is copied when script_path is given."""
        atoms = _make_test_atoms()
        script = tmp_path / "my_script.sh"
        script.write_text("#!/bin/bash\necho hello\n")
        workdir = generate_bader_workdir(
            atoms, tmp_path / "out", script_path=script,
            generate_potcar=False,
        )
        assert (workdir / "script.sh").is_file()
        assert (workdir / "script.sh").read_text() == script.read_text()

    def test_script_missing_raises(self, tmp_path):
        """FileNotFoundError raised when script_path does not exist."""
        atoms = _make_test_atoms()
        with pytest.raises(FileNotFoundError, match="not found"):
            generate_bader_workdir(
                atoms, tmp_path,
                script_path=tmp_path / "nonexistent.sh",
                generate_potcar=False,
            )

    def test_skip_potcar_when_disabled(self, tmp_path):
        """generate_potcar=False does not attempt vaspkit."""
        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, generate_potcar=False,
        )
        assert not (workdir / "POTCAR").exists()

    def test_vaspkit_not_found_raises(self, tmp_path, monkeypatch):
        """BaderGenError raised when vaspkit is not in PATH."""
        import shutil as _shutil

        monkeypatch.setattr(_shutil, "which", lambda cmd: None)
        atoms = _make_test_atoms()
        with pytest.raises(BaderGenError, match="vaspkit not found"):
            generate_bader_workdir(
                atoms, tmp_path, generate_potcar=True,
            )

    def test_custom_workdir_name(self, tmp_path):
        """Custom workdir_name is used."""
        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, workdir_name="my_bader",
            generate_potcar=False,
        )
        assert workdir == tmp_path / "my_bader"
        assert workdir.is_dir()

    def test_script_from_config(self, tmp_path, monkeypatch):
        """script_path=None reads from persistent config."""
        from md_analysis import config as _config

        script = tmp_path / "configured_script.sh"
        script.write_text("#!/bin/bash\necho configured\n")
        cfg_path = tmp_path / "test_config.json"
        _config.set_config(_config.KEY_VASP_SCRIPT_PATH, str(script),
                           config_path=cfg_path)

        monkeypatch.setattr(
            "md_analysis.scripts.BaderGen.get_config",
            lambda key, **kw: _config.get_config(key, config_path=cfg_path),
        )

        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path / "out", script_path=None,
            generate_potcar=False,
        )
        assert (workdir / "script.sh").is_file()
        assert (workdir / "script.sh").read_text() == script.read_text()

    def test_element_order(self, tmp_path):
        """Custom element_order is reflected in the IndexMap."""
        atoms = _make_test_atoms()
        workdir = generate_bader_workdir(
            atoms, tmp_path, element_order=("H", "O", "Cu"),
            generate_potcar=False,
        )
        imap = read_index_map_from_poscar(workdir / "POSCAR")
        assert imap.element_order == ("H", "O", "Cu")
