"""pymc-nmr's relaxed energy agrees with ASE relaxation used directly."""

import pytest
from ase.filters import UnitCellFilter
from ase.io import read
from ase.optimize import LBFGS

from conftest import requires_model


@requires_model
def test_pymc_relaxed_energy_matches_ase(field, basin, calculator):
    energy_pymc = field.calculate_energy_relax(
        basin.create_atoms_object(), "lbfgs", 1000, 1.0e-3, "conp", False
    )

    atoms = read("basin.xyz")
    atoms.calc = calculator
    mask = [1, 1, 1, 1, 1, 1]
    converged = LBFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=1.0e-3, steps=1000)
    assert converged
    energy_ase = atoms.get_potential_energy()

    assert energy_ase == pytest.approx(energy_pymc)
