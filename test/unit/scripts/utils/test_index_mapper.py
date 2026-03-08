"""Tests for md_analysis.scripts.utils.IndexMapper."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

from md_analysis.scripts.utils.IndexMapper import (
    IndexMap,
    IndexMapError,
    IndexMapParseError,
    compute_index_map,
    decode_comment_line,
    encode_comment_line,
    read_index_map_from_poscar,
    remap_array,
    write_poscar_with_map,
)


def _make_atoms(symbols: list[str]) -> Atoms:
    """Create a minimal Atoms object with random positions in a cubic cell."""
    rng = np.random.default_rng(42)
    n = len(symbols)
    positions = rng.random((n, 3)) * 10.0
    return Atoms(symbols=symbols, positions=positions, cell=[10, 10, 26], pbc=True)


# ---------------------------------------------------------------------------
# compute_index_map
# ---------------------------------------------------------------------------


class TestComputeIndexMap:
    def test_compute_default_order(self):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(atoms, frame=0, source="test.xyz")
        assert m.element_order == ("Cu", "O", "H")
        assert m.n_atoms == 5
        # POSCAR order should group: Cu(0), Cu(2), O(1), H(3), H(4)
        expected_p2x = np.array([0, 2, 1, 3, 4])
        np.testing.assert_array_equal(m.poscar_to_xyz, expected_p2x)

    def test_compute_custom_order(self):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(
            atoms, frame=1, source="t.xyz", element_order=("H", "O", "Cu")
        )
        assert m.element_order == ("H", "O", "Cu")
        # POSCAR order: H(3), H(4), O(1), Cu(0), Cu(2)
        expected_p2x = np.array([3, 4, 1, 0, 2])
        np.testing.assert_array_equal(m.poscar_to_xyz, expected_p2x)

    def test_bijection_property(self):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        n = m.n_atoms
        np.testing.assert_array_equal(
            m.poscar_to_xyz[m.xyz_to_poscar], np.arange(n)
        )
        np.testing.assert_array_equal(
            m.xyz_to_poscar[m.poscar_to_xyz], np.arange(n)
        )

    def test_intragroup_stable_order(self):
        atoms = _make_atoms(["O", "Cu", "O", "Cu", "H", "O"])
        m = compute_index_map(atoms, frame=0, source="")
        # Default order = (O, Cu, H) by first occurrence
        # O group: indices 0, 2, 5 → should appear in this order in POSCAR
        # Cu group: indices 1, 3
        # H group: index 4
        p2x = m.poscar_to_xyz
        o_indices = p2x[:3].tolist()
        cu_indices = p2x[3:5].tolist()
        assert o_indices == [0, 2, 5]
        assert cu_indices == [1, 3]

    def test_mismatched_elements_raises(self):
        atoms = _make_atoms(["Cu", "O", "H"])
        # Missing element in element_order
        with pytest.raises(ValueError, match="not found in element_order"):
            compute_index_map(atoms, element_order=("Cu", "O"))
        # Extra element in element_order
        with pytest.raises(ValueError, match="not present in atoms"):
            compute_index_map(atoms, element_order=("Cu", "O", "H", "Ag"))

    def test_single_element_identity(self):
        atoms = _make_atoms(["Cu", "Cu", "Cu"])
        m = compute_index_map(atoms, frame=0, source="")
        np.testing.assert_array_equal(m.poscar_to_xyz, np.arange(3))
        np.testing.assert_array_equal(m.xyz_to_poscar, np.arange(3))


# ---------------------------------------------------------------------------
# encode / decode comment line
# ---------------------------------------------------------------------------


class TestCommentLine:
    def _roundtrip_map(self) -> IndexMap:
        atoms = _make_atoms(["Cu", "Ag", "O", "H", "Cu", "H"])
        return compute_index_map(atoms, frame=42, source="./md-pos-1.xyz")

    def test_comment_roundtrip(self):
        m = self._roundtrip_map()
        line = encode_comment_line(m)
        m2 = decode_comment_line(line)
        assert m2.frame == m.frame
        assert m2.source == m.source
        assert m2.n_atoms == m.n_atoms
        assert m2.element_order == m.element_order
        np.testing.assert_array_equal(m2.poscar_to_xyz, m.poscar_to_xyz)
        np.testing.assert_array_equal(m2.xyz_to_poscar, m.xyz_to_poscar)

    def test_comment_special_chars(self):
        atoms = _make_atoms(["Cu", "O"])
        m = compute_index_map(atoms, frame=0, source="path with spaces/file(1).xyz")
        line = encode_comment_line(m)
        m2 = decode_comment_line(line)
        assert m2.source == "path with spaces/file(1).xyz"

    def test_comment_invalid_tag_raises(self):
        with pytest.raises(IndexMapParseError, match="does not start with"):
            decode_comment_line("some random POSCAR comment")

    def test_comment_missing_field_raises(self):
        with pytest.raises(IndexMapParseError, match="Cannot parse"):
            decode_comment_line("md_analysis::v1 frame=0")

    def test_comment_invalid_permutation_raises(self):
        atoms = _make_atoms(["Cu", "O", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        line = encode_comment_line(m)
        # Corrupt the base64 payload to create an invalid permutation
        # Replace p2x value with one encoding [0, 0, 0]
        import base64

        bad_p2x = base64.b64encode(
            np.array([0, 0, 0], dtype=np.uint16).tobytes()
        ).decode()
        line = line.rsplit("p2x=", 1)[0] + f"p2x={bad_p2x}"
        with pytest.raises(IndexMapParseError, match="not a valid permutation"):
            decode_comment_line(line)


# ---------------------------------------------------------------------------
# write / read POSCAR roundtrip
# ---------------------------------------------------------------------------


class TestPoscarRoundtrip:
    def test_write_read_poscar_roundtrip(self, tmp_path):
        atoms = _make_atoms(["Cu", "Ag", "O", "H", "H", "Cu", "O", "Ag"])
        m = compute_index_map(atoms, frame=7, source="traj.xyz")
        out = tmp_path / "POSCAR"
        write_poscar_with_map(atoms, out, m)

        m2 = read_index_map_from_poscar(out)
        np.testing.assert_array_equal(m2.poscar_to_xyz, m.poscar_to_xyz)
        np.testing.assert_array_equal(m2.xyz_to_poscar, m.xyz_to_poscar)
        assert m2.frame == 7
        assert m2.source == "traj.xyz"
        assert m2.element_order == m.element_order

    def test_positions_match_after_remap(self, tmp_path):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        out = tmp_path / "POSCAR"
        write_poscar_with_map(atoms, out, m)

        from ase.io import read as ase_read

        atoms_poscar = ase_read(str(out), format="vasp")
        # Remap POSCAR positions back to XYZ order
        pos_remapped = remap_array(
            atoms_poscar.get_positions(), m, "poscar_to_xyz"
        )
        np.testing.assert_allclose(pos_remapped, atoms.get_positions(), atol=1e-4)

    def test_ase_readable(self, tmp_path):
        """ASE can read the written POSCAR without errors."""
        atoms = _make_atoms(["Cu", "O", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        out = tmp_path / "POSCAR"
        write_poscar_with_map(atoms, out, m)

        from ase.io import read as ase_read

        atoms_read = ase_read(str(out), format="vasp")
        assert len(atoms_read) == len(atoms)


# ---------------------------------------------------------------------------
# remap_array
# ---------------------------------------------------------------------------


class TestRemapArray:
    def test_remap_double_trip_identity(self):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        data = np.arange(5, dtype=float)

        to_poscar = remap_array(data, m, "xyz_to_poscar")
        back = remap_array(to_poscar, m, "poscar_to_xyz")
        np.testing.assert_array_equal(back, data)

    def test_remap_multidimensional(self):
        atoms = _make_atoms(["Cu", "O", "Cu", "H", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        positions = atoms.get_positions()  # shape (5, 3)

        remapped = remap_array(positions, m, "xyz_to_poscar")
        assert remapped.shape == (5, 3)
        # Each row in remapped should come from the correct XYZ atom
        for j in range(5):
            xyz_idx = m.poscar_to_xyz[j]
            np.testing.assert_array_equal(remapped[j], positions[xyz_idx])

    def test_remap_length_mismatch_raises(self):
        atoms = _make_atoms(["Cu", "O", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        with pytest.raises(IndexMapError, match="data length"):
            remap_array(np.arange(5), m, "xyz_to_poscar")

    def test_remap_invalid_direction_raises(self):
        atoms = _make_atoms(["Cu", "O", "H"])
        m = compute_index_map(atoms, frame=0, source="")
        with pytest.raises(ValueError, match="Unknown direction"):
            remap_array(np.arange(3), m, "bad_direction")
