"""Geometry relaxation through ASE reaches the reference relaxed energy."""

import pytest
from ase.filters import UnitCellFilter
from ase.io import read
from ase.optimize import LBFGS

from conftest import requires_model


@requires_model
def test_relaxation(calculator):
    atoms = read("basin.xyz")
    atoms.calc = calculator

    mask = [1, 1, 1, 1, 1, 1]
    converged = LBFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=1.0e-3, steps=1000)

    assert converged
    assert atoms.get_potential_energy() == pytest.approx(-2205.071010)
