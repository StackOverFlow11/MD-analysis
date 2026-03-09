"""Unit tests for md_analysis.utils.CubeParser."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.utils.CubeParser import (
    CubeHeader,
    discover_cube_files,
    extract_step_from_cube_filename,
    plane_avg_phi_z_ev,
    read_cube_header_and_values,
    slab_average_potential_ev,
    z_coords_ang,
)
from md_analysis.utils.config import BOHR_TO_ANG, HA_TO_EV


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_header(**overrides) -> CubeHeader:
    """Build a CubeHeader with sensible defaults (all in bohr)."""
    defaults = dict(
        natoms=1,
        origin_bohr=np.array([0.0, 0.0, 0.0]),
        nx=2,
        ny=2,
        nz=4,
        vx_bohr=np.array([1.0, 0.0, 0.0]),
        vy_bohr=np.array([0.0, 1.0, 0.0]),
        vz_bohr=np.array([0.0, 0.0, 1.0]),
    )
    defaults.update(overrides)
    return CubeHeader(**defaults)


def _write_cube(
    path: Path,
    *,
    natoms: int = 1,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    nx: int = 2,
    ny: int = 2,
    nz: int = 4,
    vx: tuple[float, float, float] = (1.0, 0.0, 0.0),
    vy: tuple[float, float, float] = (0.0, 1.0, 0.0),
    vz: tuple[float, float, float] = (0.0, 0.0, 1.0),
    values: np.ndarray | None = None,
    fortran_d: bool = False,
) -> Path:
    """Write a minimal Gaussian cube file for testing."""
    if values is None:
        values = np.ones(nx * ny * nz, dtype=float)
    lines: list[str] = []
    lines.append("comment line 1\n")
    lines.append("comment line 2\n")
    lines.append(f"  {natoms}  {origin[0]:.6f}  {origin[1]:.6f}  {origin[2]:.6f}\n")
    lines.append(f"  {nx}  {vx[0]:.6f}  {vx[1]:.6f}  {vx[2]:.6f}\n")
    lines.append(f"  {ny}  {vy[0]:.6f}  {vy[1]:.6f}  {vy[2]:.6f}\n")
    lines.append(f"  {nz}  {vz[0]:.6f}  {vz[1]:.6f}  {vz[2]:.6f}\n")
    for _ in range(natoms):
        lines.append("  6  6.000000  0.000000  0.000000  0.000000\n")

    # Write data values
    val_strs: list[str] = []
    for v in values:
        s = f"{v:13.5E}"
        if fortran_d:
            s = s.replace("E", "D")
        val_strs.append(s)
    # Write 6 values per line (like CP2K)
    for i in range(0, len(val_strs), 6):
        lines.append("  ".join(val_strs[i : i + 6]) + "\n")

    path.write_text("".join(lines))
    return path


# ===========================================================================
# TestReadCubeHeaderAndValues
# ===========================================================================


class TestReadCubeHeaderAndValues:

    def test_happy_path(self, tmp_path: Path):
        cube = _write_cube(tmp_path / "test.cube")
        header, values = read_cube_header_and_values(cube)

        assert header.natoms == 1
        assert header.nx == 2
        assert header.ny == 2
        assert header.nz == 4
        assert values.shape == (2 * 2 * 4,)
        np.testing.assert_allclose(values, 1.0)

    def test_origin_and_voxels(self, tmp_path: Path):
        cube = _write_cube(
            tmp_path / "test.cube",
            origin=(1.5, 2.5, 3.5),
            vx=(0.5, 0.0, 0.0),
            vy=(0.0, 0.7, 0.0),
            vz=(0.0, 0.0, 0.3),
        )
        header, _ = read_cube_header_and_values(cube)

        np.testing.assert_allclose(header.origin_bohr, [1.5, 2.5, 3.5])
        np.testing.assert_allclose(header.vx_bohr, [0.5, 0.0, 0.0])
        np.testing.assert_allclose(header.vy_bohr, [0.0, 0.7, 0.0])
        np.testing.assert_allclose(header.vz_bohr, [0.0, 0.0, 0.3])

    def test_fortran_d_exponent(self, tmp_path: Path):
        vals = np.array([1.5, 2.0, 3.0, 0.001] * 4)  # 16 values for 2x2x4
        cube = _write_cube(
            tmp_path / "test.cube", values=vals, fortran_d=True,
        )
        _, parsed = read_cube_header_and_values(cube)
        np.testing.assert_allclose(parsed, vals, rtol=1e-4)

    def test_multi_atom(self, tmp_path: Path):
        vals = np.arange(16, dtype=float)
        cube = _write_cube(
            tmp_path / "test.cube", natoms=3, values=vals,
        )
        header, parsed = read_cube_header_and_values(cube)

        assert header.natoms == 3
        np.testing.assert_allclose(parsed, vals, atol=1e-4)

    def test_header_too_short_raises(self, tmp_path: Path):
        cube = tmp_path / "bad.cube"
        cube.write_text(
            "comment 1\n"
            "comment 2\n"
            "1 0.0 0.0\n"  # only 3 tokens — need 4
        )
        with pytest.raises(ValueError, match="too short|index out of range"):
            read_cube_header_and_values(cube)

    def test_data_size_mismatch_raises(self, tmp_path: Path):
        vals = np.ones(10)  # expected 2*2*4=16
        cube = _write_cube(tmp_path / "test.cube", values=vals)
        with pytest.raises(ValueError, match="data size mismatch|size mismatch"):
            read_cube_header_and_values(cube)

    def test_frozen_header(self, tmp_path: Path):
        cube = _write_cube(tmp_path / "test.cube")
        header, _ = read_cube_header_and_values(cube)
        with pytest.raises(AttributeError):
            header.nx = 999  # type: ignore[misc]


# ===========================================================================
# TestSlabAveragePotentialEv
# ===========================================================================


class TestSlabAveragePotentialEv:

    def test_cell_center(self):
        header = _make_header()
        values = np.full(2 * 2 * 4, 0.1)  # uniform 0.1 Ha
        phi, info = slab_average_potential_ev(header, values, thickness_ang=2.0)

        assert phi == pytest.approx(0.1 * HA_TO_EV, rel=1e-10)
        assert info["z_center_source"] == "cell"

    def test_explicit_z_center(self):
        header = _make_header()
        values = np.full(2 * 2 * 4, 0.1)
        _, info = slab_average_potential_ev(
            header, values, thickness_ang=2.0, z_center_ang=1.0,
        )
        assert info["z_center_source"] == "interface"

    def test_uniform_zero_std(self):
        header = _make_header()
        values = np.full(2 * 2 * 4, 0.05)
        _, info = slab_average_potential_ev(header, values, thickness_ang=3.0)
        assert info["phi_z_std_ev"] == pytest.approx(0.0, abs=1e-14)

    def test_nonuniform_nonzero_std(self):
        header = _make_header(nz=8)
        # Build field where the z-plane average varies linearly
        field = np.zeros((2, 2, 8))
        for iz in range(8):
            field[:, :, iz] = iz * 0.1  # 0.0, 0.1, 0.2, ... Ha
        values = field.ravel()
        _, info = slab_average_potential_ev(header, values, thickness_ang=6.0)
        assert info["phi_z_std_ev"] > 0

    def test_n_slices_count(self):
        # dz = |vz| = 1.0 bohr = 0.529... Ang
        # thickness = 1.5 Ang → half = 0.75 Ang ≈ 1.417 bohr
        # z_center_rel = 0.5*4 = 2.0 bohr; z_coords = 0.5,1.5,2.5,3.5
        # dist from 2.0: 1.5, 0.5, 0.5, 1.5  (raw); PBC min: 1.5,0.5,0.5,1.5
        # half_thickness_bohr = 0.75/0.529... ≈ 1.417 → slices: dist<=1.417 → 2 slices (0.5,0.5)
        header = _make_header()
        values = np.full(16, 0.1)
        _, info = slab_average_potential_ev(header, values, thickness_ang=1.5)
        assert info["n_slices"] == 2

    def test_zero_dz_raises(self):
        header = _make_header(vz_bohr=np.array([0.0, 0.0, 0.0]))
        values = np.full(16, 0.1)
        with pytest.raises(ValueError, match="Invalid dz"):
            slab_average_potential_ev(header, values, thickness_ang=2.0)

    def test_zero_slices_raises(self):
        header = _make_header()
        values = np.full(16, 0.1)
        # dz = 1.0 bohr ≈ 0.529 Ang. Thickness = 0.01 Ang → 0 slices
        with pytest.raises(ValueError, match="selects 0"):
            slab_average_potential_ev(header, values, thickness_ang=0.01)


# ===========================================================================
# TestPlaneAvgPhiZEv
# ===========================================================================


class TestPlaneAvgPhiZEv:

    def test_shape_and_conversion(self):
        header = _make_header()
        val_ha = 0.2
        values = np.full(2 * 2 * 4, val_ha)

        phi_z = plane_avg_phi_z_ev(header, values)
        assert phi_z.shape == (4,)
        np.testing.assert_allclose(phi_z, val_ha * HA_TO_EV)

    def test_xy_averaging(self):
        header = _make_header()
        field = np.zeros((2, 2, 4))
        # Set different xy values for each z-slice
        field[0, 0, :] = 1.0
        field[0, 1, :] = 2.0
        field[1, 0, :] = 3.0
        field[1, 1, :] = 4.0
        values = field.ravel()

        phi_z = plane_avg_phi_z_ev(header, values)
        # xy average = (1+2+3+4)/4 = 2.5 Ha at each z
        expected = 2.5 * HA_TO_EV
        np.testing.assert_allclose(phi_z, expected)


# ===========================================================================
# TestZCoordsAng
# ===========================================================================


class TestZCoordsAng:

    def test_shape_and_values(self):
        header = _make_header()  # origin_z=0, dz=1.0 bohr, nz=4
        z = z_coords_ang(header)

        assert z.shape == (4,)
        dz_ang = 1.0 * BOHR_TO_ANG
        expected = np.array([0.5, 1.5, 2.5, 3.5]) * dz_ang
        np.testing.assert_allclose(z, expected)

    def test_nonzero_origin(self):
        origin_z_bohr = 5.0
        header = _make_header(origin_bohr=np.array([0.0, 0.0, origin_z_bohr]))
        z = z_coords_ang(header)

        origin_z_ang = origin_z_bohr * BOHR_TO_ANG
        dz_ang = 1.0 * BOHR_TO_ANG
        expected = origin_z_ang + (np.arange(4) + 0.5) * dz_ang
        np.testing.assert_allclose(z, expected)


# ===========================================================================
# TestExtractStepFromCubeFilename
# ===========================================================================


class TestExtractStepFromCubeFilename:

    def test_typical(self):
        assert extract_step_from_cube_filename(Path("md-POTENTIAL-v_hartree-1_1050.cube")) == 1050

    def test_zero_step(self):
        assert extract_step_from_cube_filename(Path("prefix_0.cube")) == 0

    def test_no_match(self):
        assert extract_step_from_cube_filename(Path("some_data.dat")) is None

    def test_multi_underscore(self):
        assert extract_step_from_cube_filename(Path("a_b_999.cube")) == 999


# ===========================================================================
# TestDiscoverCubeFiles
# ===========================================================================


class TestDiscoverCubeFiles:

    def test_happy_path(self, tmp_path: Path):
        """Sorted list of matching cube files is returned."""
        for name in ["c_3.cube", "a_1.cube", "b_2.cube"]:
            _write_cube(tmp_path / name)
        result = discover_cube_files("*.cube", workdir=tmp_path)
        assert [p.name for p in result] == ["a_1.cube", "b_2.cube", "c_3.cube"]

    def test_no_match_raises(self, tmp_path: Path):
        """FileNotFoundError when no files match."""
        with pytest.raises(FileNotFoundError, match="No cube files"):
            discover_cube_files("*.cube", workdir=tmp_path)

    def test_frame_slicing(self, tmp_path: Path):
        """frame_start/end/step correctly slice the sorted list."""
        for i in range(5):
            _write_cube(tmp_path / f"cube_{i}.cube")
        result = discover_cube_files(
            "*.cube", workdir=tmp_path,
            frame_start=1, frame_end=4, frame_step=2,
        )
        assert [p.name for p in result] == ["cube_1.cube", "cube_3.cube"]

    def test_default_workdir(self, tmp_path: Path, monkeypatch):
        """When workdir is None, cwd is used."""
        _write_cube(tmp_path / "test_0.cube")
        monkeypatch.chdir(tmp_path)
        result = discover_cube_files("*.cube")
        assert len(result) == 1
        assert result[0].name == "test_0.cube"
