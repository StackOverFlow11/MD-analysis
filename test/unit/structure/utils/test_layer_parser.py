"""Unit tests for layer parsing utilities."""

from __future__ import annotations

import numpy as np
import pytest

ase = pytest.importorskip("ase")
from ase import Atoms

from scripts.structure.utils.LayerParser import detect_interface_layers


@pytest.fixture
def simple_slab_with_environment() -> Atoms:
    """
    Build two Cu layers plus non-metal atoms on both sides along z.
    """
    symbols = [
        # Layer 1 (low z)
        "Cu",
        "Cu",
        "Cu",
        "Cu",
        # Layer 2 (high z)
        "Cu",
        "Cu",
        "Cu",
        "Cu",
        # Environment atoms near both sides
        "O",
        "O",
    ]
    positions = np.array(
        [
            [2.0, 2.0, 2.0],
            [2.0, 8.0, 2.0],
            [8.0, 2.0, 2.0],
            [8.0, 8.0, 2.0],
            [2.0, 2.0, 8.0],
            [2.0, 8.0, 8.0],
            [8.0, 2.0, 8.0],
            [8.0, 8.0, 8.0],
            [5.0, 5.0, 0.8],
            [5.0, 5.0, 9.2],
        ],
        dtype=float,
    )
    return Atoms(symbols=symbols, positions=positions, cell=[10.0, 10.0, 10.0], pbc=[True, True, True])


def test_detect_interface_layers_marks_two_sides(simple_slab_with_environment: Atoms) -> None:
    result = detect_interface_layers(
        simple_slab_with_environment,
        metal_symbols={"Cu"},
        normal="c",
        layer_tol_A=0.6,
        n_interface_layers=1,
    )
    assert len(result.metal_indices) == 8
    assert len(result.metal_layers_sorted) == 2

    interfaces = result.interface_layers()
    assert len(interfaces) == 2
    z_signs = {int(np.sign(layer.normal_unit[2])) for layer in interfaces if layer.normal_unit is not None}
    assert z_signs == {-1, 1}
