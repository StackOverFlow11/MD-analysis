"""Tests for md_analysis.scripts._vasp_md_source — VASP MD frame source."""

from __future__ import annotations

from pathlib import Path

import pytest

from md_analysis.scripts._frame_selector import FrameSelection
from md_analysis.scripts._vasp_md_source import (
    VaspMdSourceError,
    iter_vasp_md_frames,
    read_element_symbols,
    select_trajectory_file,
)


# ---------------------------------------------------------------------------
# Fixture builders — minimal hand-written VASP files that ase.io accepts
# ---------------------------------------------------------------------------

# Orthorhombic test cell (Angstrom), row-vector convention.
CELL_ORTHO = ((8.0, 0.0, 0.0), (0.0, 8.0, 0.0), (0.0, 0.0, 10.0))
SPECIES = ("Pt", "O", "H")
COUNTS = (2, 1, 2)
N_ATOMS = sum(COUNTS)  # 5


def _write_poscar(
    path: Path,
    cell=CELL_ORTHO,
    species=SPECIES,
    counts=COUNTS,
    vasp5: bool = True,
) -> None:
    """Write a minimal POSCAR (Direct coords); ``vasp5=False`` drops the
    symbols line (VASP4-style header)."""
    lines = ["test slab", "1.0"]
    lines += ["  " + "  ".join(f"{v:.8f}" for v in vec) for vec in cell]
    if vasp5:
        lines.append("  " + "  ".join(f"{s:>3s}" for s in species))
    lines.append("  " + "  ".join(f"{n:3d}" for n in counts))
    lines.append("Direct")
    for i in range(sum(counts)):
        f = 0.1 * (i + 1)
        lines.append(f"  {f:.6f}  {f:.6f}  {f:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_xdatcar(
    path: Path,
    n_frames: int,
    cell=CELL_ORTHO,
    species=SPECIES,
    counts=COUNTS,
    vasp5: bool = True,
) -> None:
    """Write a minimal XDATCAR; ``vasp5=False`` drops the species line.

    Atom ``i`` of frame ``k`` sits at fractional ``0.1*(i+1) + 0.001*k``
    in every direction, so frames are distinguishable by position.
    """
    lines = ["test slab", "1.0"]
    lines += ["  " + "  ".join(f"{v:.8f}" for v in vec) for vec in cell]
    if vasp5:
        lines.append("  " + "  ".join(f"{s:>3s}" for s in species))
    lines.append("  " + "  ".join(f"{n:3d}" for n in counts))
    for frame in range(n_frames):
        lines.append(f"Direct configuration=  {frame + 1:>6d}")
        for i in range(sum(counts)):
            f = 0.1 * (i + 1) + 0.001 * frame
            lines.append(f"  {f:.6f}  {f:.6f}  {f:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_vasprun(path: Path, n_frames: int = 3, cell=CELL_ORTHO) -> None:
    """Write a minimal vasprun.xml accepted by ``ase.io.read(format='vasp-xml')``.

    ase's reader requires <atominfo>, <kpoints> and one <calculation> per
    frame (each with an <scstep> energy block and a <structure>). The
    first atom's fractional z distinguishes frames: ``0.10 + 0.01*k``.
    """
    basis = "\n".join(
        f"        <v>{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}</v>" for v in cell
    )

    def _calc(frame: int) -> str:
        z0 = 0.10 + 0.01 * frame
        return f"""  <calculation>
    <scstep>
      <energy>
        <i name="e_0_energy">-100.0</i>
        <i name="e_fr_energy">-100.0</i>
      </energy>
    </scstep>
    <structure>
      <crystal><varray name="basis">
{basis}
      </varray></crystal>
      <varray name="positions">
        <v>0.10 0.10 {z0:.6f}</v>
        <v>0.20 0.20 0.20</v>
        <v>0.30 0.30 0.30</v>
      </varray>
    </structure>
    <energy><i name="e_fr_energy">-100.0</i></energy>
  </calculation>"""

    text = f"""<?xml version="1.0" encoding="UTF-8"?>
<modeling>
  <atominfo>
    <array name="atoms"><set>
        <rc><c>Pt</c><c>Pt_1</c></rc>
        <rc><c>Pt</c><c>Pt_2</c></rc>
        <rc><c>O</c><c>O_1</c></rc>
    </set></array>
  </atominfo>
  <kpoints>
    <varray name="kpointlist"><v>0.0 0.0 0.0</v></varray>
    <varray name="weights"><v>1.0</v></varray>
  </kpoints>
  <structure name="initialpos">
    <crystal><varray name="basis">
{basis}
    </varray></crystal>
    <varray name="positions">
      <v>0.10 0.10 0.10</v>
      <v>0.20 0.20 0.20</v>
      <v>0.30 0.30 0.30</v>
    </varray>
  </structure>
{chr(10).join(_calc(k) for k in range(n_frames))}
</modeling>
"""
    path.write_text(text, encoding="utf-8")


def _make_xdatcar_dir(tmp_path: Path, n_frames: int = 3, **kwargs) -> Path:
    """tmp_path with an XDATCAR + matching POSCAR."""
    _write_xdatcar(tmp_path / "XDATCAR", n_frames, **kwargs)
    _write_poscar(tmp_path / "POSCAR", cell=kwargs.get("cell", CELL_ORTHO))
    return tmp_path


# ---------------------------------------------------------------------------
# select_trajectory_file
# ---------------------------------------------------------------------------


class TestSelectTrajectoryFile:
    def test_prefers_vasprun_over_xdatcar(self, tmp_path):
        (tmp_path / "vasprun.xml").write_text("dummy", encoding="utf-8")
        (tmp_path / "XDATCAR").write_text("dummy", encoding="utf-8")
        assert select_trajectory_file(tmp_path).name == "vasprun.xml"

    def test_falls_back_to_xdatcar(self, tmp_path):
        (tmp_path / "XDATCAR").write_text("dummy", encoding="utf-8")
        assert select_trajectory_file(tmp_path).name == "XDATCAR"

    def test_missing_dir_raises(self, tmp_path):
        with pytest.raises(VaspMdSourceError, match="not found"):
            select_trajectory_file(tmp_path / "nope")

    def test_no_trajectory_lists_candidates(self, tmp_path):
        (tmp_path / "POSCAR").write_text("dummy", encoding="utf-8")
        with pytest.raises(VaspMdSourceError) as excinfo:
            select_trajectory_file(tmp_path)
        msg = str(excinfo.value)
        assert "vasprun.xml" in msg
        assert "XDATCAR" in msg


# ---------------------------------------------------------------------------
# iter_vasp_md_frames — XDATCAR route
# ---------------------------------------------------------------------------


class TestXdatcarFrames:
    def test_frames_symbols_cell_pbc(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=3)
        frames = list(iter_vasp_md_frames(d, FrameSelection()))
        assert [idx for idx, _ in frames] == [0, 1, 2]
        atoms = frames[0][1]
        assert atoms.get_chemical_symbols() == ["Pt", "Pt", "O", "H", "H"]
        assert len(atoms) == N_ATOMS
        cell = atoms.get_cell()
        assert cell[0, 0] == pytest.approx(8.0)
        assert cell[1, 1] == pytest.approx(8.0)
        assert cell[2, 2] == pytest.approx(10.0)
        assert atoms.pbc.all()

    def test_cell_comes_from_xdatcar_not_poscar(self, tmp_path):
        """The yielded cell is the trajectory's, not the POSCAR's."""
        other_cell = ((12.0, 0.0, 0.0), (0.0, 12.0, 0.0), (0.0, 0.0, 14.0))
        _write_xdatcar(tmp_path / "XDATCAR", 2, cell=CELL_ORTHO)
        _write_poscar(tmp_path / "POSCAR", cell=other_cell)
        frames = list(iter_vasp_md_frames(tmp_path, FrameSelection()))
        assert frames[0][1].get_cell()[0, 0] == pytest.approx(8.0)

    def test_frame_positions_track_frames(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=3)
        frames = list(iter_vasp_md_frames(d, FrameSelection()))
        for k, (_, atoms) in enumerate(frames):
            expected_z = (0.1 + 0.001 * k) * 10.0
            assert atoms.get_positions()[0][2] == pytest.approx(expected_z)

    def test_index_slicing(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=5)
        sel = FrameSelection(frame_start=1, frame_end=4, frame_step=2)
        frames = list(iter_vasp_md_frames(d, sel))
        assert [idx for idx, _ in frames] == [1, 3]
        # Positions confirm the right frames were sliced out
        assert frames[0][1].get_positions()[0][2] == pytest.approx(
            (0.1 + 0.001 * 1) * 10.0
        )
        assert frames[1][1].get_positions()[0][2] == pytest.approx(
            (0.1 + 0.001 * 3) * 10.0
        )

    def test_step_only(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=5)
        sel = FrameSelection(frame_step=2)
        assert [idx for idx, _ in iter_vasp_md_frames(d, sel)] == [0, 2, 4]

    def test_vasp4_style_xdatcar_uses_poscar_elements(self, tmp_path):
        """No species line in XDATCAR → elements come from the POSCAR."""
        d = _make_xdatcar_dir(tmp_path, n_frames=2, vasp5=False)
        frames = list(iter_vasp_md_frames(d, FrameSelection()))
        assert [idx for idx, _ in frames] == [0, 1]
        assert frames[0][1].get_chemical_symbols() == ["Pt", "Pt", "O", "H", "H"]

    def test_missing_poscar_raises(self, tmp_path):
        _write_xdatcar(tmp_path / "XDATCAR", 2)
        with pytest.raises(VaspMdSourceError, match="POSCAR"):
            list(iter_vasp_md_frames(tmp_path, FrameSelection()))

    def test_vasp4_style_poscar_rejected(self, tmp_path):
        """POSCAR without a symbols line cannot anchor element order."""
        _write_xdatcar(tmp_path / "XDATCAR", 2)
        _write_poscar(tmp_path / "POSCAR", vasp5=False)
        with pytest.raises(VaspMdSourceError, match="symbol"):
            list(iter_vasp_md_frames(tmp_path, FrameSelection()))

    def test_species_mismatch_with_poscar_raises(self, tmp_path):
        _write_xdatcar(tmp_path / "XDATCAR", 2, species=("Au", "O", "H"))
        _write_poscar(tmp_path / "POSCAR")
        with pytest.raises(VaspMdSourceError, match="disagrees"):
            list(iter_vasp_md_frames(tmp_path, FrameSelection()))


# ---------------------------------------------------------------------------
# iter_vasp_md_frames — vasprun.xml route
# ---------------------------------------------------------------------------


class TestVasprunFrames:
    def test_frames_symbols_cell_pbc(self, tmp_path):
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=3)
        frames = list(iter_vasp_md_frames(tmp_path, FrameSelection()))
        assert [idx for idx, _ in frames] == [0, 1, 2]
        atoms = frames[0][1]
        assert atoms.get_chemical_symbols() == ["Pt", "Pt", "O"]
        cell = atoms.get_cell()
        assert cell[0, 0] == pytest.approx(8.0)
        assert cell[2, 2] == pytest.approx(10.0)
        assert atoms.pbc.all()

    def test_frame_positions_track_frames(self, tmp_path):
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=3)
        frames = list(iter_vasp_md_frames(tmp_path, FrameSelection()))
        for k, (_, atoms) in enumerate(frames):
            expected_z = (0.10 + 0.01 * k) * 10.0
            assert atoms.get_positions()[0][2] == pytest.approx(expected_z)

    def test_index_slicing(self, tmp_path):
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=5)
        sel = FrameSelection(frame_start=1, frame_end=4, frame_step=2)
        frames = list(iter_vasp_md_frames(tmp_path, sel))
        assert [idx for idx, _ in frames] == [1, 3]

    def test_no_poscar_needed(self, tmp_path):
        """The vasprun route has no POSCAR requirement."""
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=1)
        frames = list(iter_vasp_md_frames(tmp_path, FrameSelection()))
        assert len(frames) == 1

    def test_vasprun_takes_priority_over_xdatcar(self, tmp_path):
        """Both files present → vasprun wins; its cell is yielded.

        No POSCAR is provided, so reaching the XDATCAR route would raise
        here — proving the vasprun route was taken.
        """
        other_cell = ((12.0, 0.0, 0.0), (0.0, 12.0, 0.0), (0.0, 0.0, 14.0))
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=2)
        _write_xdatcar(tmp_path / "XDATCAR", 2, cell=other_cell)
        frames = list(iter_vasp_md_frames(tmp_path, FrameSelection()))
        assert len(frames) == 2
        assert frames[0][1].get_cell()[0, 0] == pytest.approx(8.0)
        assert frames[0][1].get_chemical_symbols() == ["Pt", "Pt", "O"]


# ---------------------------------------------------------------------------
# Orthogonality guard
# ---------------------------------------------------------------------------


class TestOrthogonality:
    def test_non_orthogonal_cell_raises(self, tmp_path):
        cell = ((8.0, 0.0, 0.0), (0.5, 8.0, 0.0), (0.0, 0.0, 10.0))
        d = _make_xdatcar_dir(tmp_path, n_frames=2, cell=cell)
        with pytest.raises(VaspMdSourceError) as excinfo:
            list(iter_vasp_md_frames(d, FrameSelection()))
        msg = str(excinfo.value)
        assert "orthorhombic" in msg
        assert "sp.inp" in msg

    def test_non_orthogonal_vasprun_raises(self, tmp_path):
        cell = ((8.0, 0.0, 0.0), (0.0, 8.0, 0.0), (0.1, 0.0, 10.0))
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=2, cell=cell)
        with pytest.raises(VaspMdSourceError, match="orthorhombic"):
            list(iter_vasp_md_frames(tmp_path, FrameSelection()))

    def test_off_diagonal_below_tolerance_passes(self, tmp_path):
        """Off-diagonals within the 1e-6 A tolerance are accepted."""
        cell = ((8.0, 1e-8, 0.0), (0.0, 8.0, 0.0), (0.0, 0.0, 10.0))
        d = _make_xdatcar_dir(tmp_path, n_frames=2, cell=cell)
        frames = list(iter_vasp_md_frames(d, FrameSelection()))
        assert len(frames) == 2

    def test_off_diagonal_above_tolerance_raises(self, tmp_path):
        cell = ((8.0, 1e-3, 0.0), (0.0, 8.0, 0.0), (0.0, 0.0, 10.0))
        d = _make_xdatcar_dir(tmp_path, n_frames=2, cell=cell)
        with pytest.raises(VaspMdSourceError, match="orthorhombic"):
            list(iter_vasp_md_frames(d, FrameSelection()))


# ---------------------------------------------------------------------------
# Time-mode rejection
# ---------------------------------------------------------------------------


class TestTimeMode:
    def test_time_mode_rejected_eagerly(self, tmp_path):
        """Time mode fails fast — even before any file resolution."""
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=100.0, time_step_fs=10.0,
        )
        # No list(): validation must raise at call time, on an empty dir
        with pytest.raises(VaspMdSourceError, match="time mode"):
            iter_vasp_md_frames(tmp_path, sel)

    def test_time_mode_message_points_to_index_mode(self, tmp_path):
        sel = FrameSelection(
            mode="time", time_start_fs=0.0, time_end_fs=100.0, time_step_fs=10.0,
        )
        with pytest.raises(VaspMdSourceError) as excinfo:
            iter_vasp_md_frames(tmp_path, sel)
        assert "frame_start" in str(excinfo.value)


# ---------------------------------------------------------------------------
# read_element_symbols
# ---------------------------------------------------------------------------


class TestReadElementSymbols:
    def test_vasprun_route(self, tmp_path):
        _write_vasprun(tmp_path / "vasprun.xml", n_frames=2)
        assert read_element_symbols(tmp_path) == ["Pt", "Pt", "O"]

    def test_xdatcar_route_reads_poscar(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=2)
        assert read_element_symbols(d) == ["Pt", "Pt", "O", "H", "H"]

    def test_xdatcar_route_vasp4_style(self, tmp_path):
        d = _make_xdatcar_dir(tmp_path, n_frames=2, vasp5=False)
        assert read_element_symbols(d) == ["Pt", "Pt", "O", "H", "H"]

    def test_xdatcar_route_requires_poscar(self, tmp_path):
        _write_xdatcar(tmp_path / "XDATCAR", 2)
        with pytest.raises(VaspMdSourceError, match="POSCAR"):
            read_element_symbols(tmp_path)

    def test_no_trajectory_raises(self, tmp_path):
        with pytest.raises(VaspMdSourceError, match="vasprun.xml"):
            read_element_symbols(tmp_path)

    def test_calculation_less_vasprun_raises(self, tmp_path):
        """Crashed run (no <calculation> blocks) → clear error, not a crash."""
        (tmp_path / "vasprun.xml").write_text(
            '<?xml version="1.0"?><modeling><atominfo>'
            "<array name='atoms'><set><rc><c>Pt</c><c>Pt_1</c></rc></set>"
            "</array></atominfo></modeling>",
            encoding="utf-8",
        )
        with pytest.raises(VaspMdSourceError, match="No frames"):
            read_element_symbols(tmp_path)
