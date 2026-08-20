"""pymc-nmr's single-point energy agrees with ASE used directly.

This is a differential test against an oracle: the Field/Config/Species path
through pymc-nmr must reproduce what plain ASE computes for the same structure
and potential.
"""

import pytest
from ase.io import read

from conftest import requires_model


@requires_model
def test_pymc_energy_matches_ase(field, basin, calculator):
    energy_pymc = field.calculate_energy(basin.create_atoms_object(), False)

    atoms = read("basin.xyz")
    atoms.calc = calculator
    energy_ase = atoms.get_potential_energy()

    assert energy_ase == pytest.approx(energy_pymc)
