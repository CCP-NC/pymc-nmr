#same as 10 on cpu but with dispersion included via janus and then directly through mace
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
atoms.calc = choose_calculator(architecture=arch, dispersion = True, model=model, precision=precsn, device=device, calc_kwargs={'dispersion' : True})
janus_energy = atoms.get_potential_energy()

atoms.calc = mace_mp(model=model, dispersion=False, default_dtype="float64", device='cpu')
mace_energy = atoms.get_potential_energy()

diff_energy = mace_energy - janus_energy

assert diff_energy == pytest.approx(0.0)

print(diff_energy)