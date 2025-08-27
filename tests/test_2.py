#tests whether some of the basics are ok - ASE, mace model, configuration
from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.io import read

from janus_core.helpers.mlip_calculators import choose_calculator

import numpy as np
import pytest


atoms = read("basin.xyz")

device = "cuda"
precsn = "float64"
arch = "mace_mp"
model = "./data/MACE-matpes-r2scan-omat-ft.model"
atoms.calc = choose_calculator(arch=arch, model=model, precision=precsn, device=device)

mask=[1,1,1,1,1,1]

flag = LBFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=1.0e-3, steps=1000)

assert flag == True

energy = atoms.get_potential_energy()

assert energy == pytest.approx(-2205.071010)

print(energy)