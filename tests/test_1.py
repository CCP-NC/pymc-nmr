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

energy = atoms.get_potential_energy()

assert energy == pytest.approx(-2184.6450253232610)

print(energy)