#same as test 1 except uses cpu
from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.io import read

from janus_core.helpers.mlip_calculators import choose_calculator
from mace.calculators import mace_mp
import numpy as np
import pytest


atoms = read("basin.xyz")

device = "cpu"
precsn = "float64"
arch = "mace_mp"
model = "./data/MACE-matpes-r2scan-omat-ft.model"
atoms.calc = choose_calculator(arch=arch, model=model, precision=precsn, device=device)

energy = atoms.get_potential_energy()

assert energy == pytest.approx(-2184.6450253232610)

print(energy)