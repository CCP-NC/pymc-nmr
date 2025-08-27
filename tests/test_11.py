#same as 10 on cpu but with dispersion included via janus and then directly through mace
from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.io import read

from janus_core.helpers.mlip_calculators import choose_calculator

from mace.calculators import mace_mp # directly get the mace calculator

import numpy as np
import pytest

try:
    import dftd3
except Exception as e:
            print(f"{e} whilst trying to import dft-d3. Check to see whether it has been installed! \n")
            exit()

atoms = read("basin.xyz")

device = "cpu"
precsn = "float64"
arch = "mace_mp"
model = "./data/MACE-matpes-r2scan-omat-ft.model"
atoms.calc = choose_calculator(arch=arch, dispersion=True, model=model, precision=precsn, device=device)
janus_energy = atoms.get_potential_energy()

atoms.calc = mace_mp(model=model, dispersion=True, default_dtype="float64", device='cpu')
mace_energy = atoms.get_potential_energy()

diff_energy = mace_energy - janus_energy

assert diff_energy == pytest.approx(0.0)

print(janus_energy, mace_energy, diff_energy)