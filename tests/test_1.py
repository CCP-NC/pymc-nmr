"""Smoke test: ASE, janus-core and the MACE model load and give the reference energy."""

import pytest
from ase.io import read

from conftest import requires_model


@requires_model
def test_single_point_energy(calculator):
    atoms = read("basin.xyz")
    atoms.calc = calculator
    assert atoms.get_potential_energy() == pytest.approx(-2184.6450253232610)
