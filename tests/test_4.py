#determines whether a relaxed energy calculated within pymc-nmr is the same as that from ase
#this checks some of the functionality on Field, Config, Species
import sys
sys.path.insert(1, '../src')

from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.io import read

from janus_core.helpers.mlip_calculators import choose_calculator

import numpy as np
import pytest

from config import Config
from JobControl import JobControl
from Species import Species
from Field import Field
from Energy import Energy



def pymc_energy():
    restart_iteration = 0
    md_time = 0.0

    job = JobControl()
    spec = Species()

    fld = Field()

    basisFileName = "basin.xyz"
    fieldFileName = "potentials"

    out_io = open("error.log", "w")
    

    #field file
    try:
        with open(fieldFileName, 'r') as in_stream:
            
            fld.readPotential(in_stream, out_io, spec)
            
    except FileNotFoundError:
        print("\n*** could not find potentials file\n")
        sys.exit(1)

    # Read data from basis file in xyz format
    basin = Config()
    restart_iteration = 0
    restart_time = 0.0 
    restart_energy = 0.0
    try:
        instream = open(basisFileName, "r")
        restart_iteration, restart_time, restart_energy = basin.read_config(instream)
        basin.setup_configuration(spec, out_io)
    except FileNotFoundError:
        print("\n*** could not find configuration file: basin.xyz \n")
        sys.exit(1)

    fld.setup()
    atoms = basin.create_atoms_object()
    energy = fld.calculate_energy_relax(atoms, "lbfgs", 1000, 1.0e-3, "conp", False)

    return energy.totalEnergy


atoms = read("basin.xyz")

device = "cuda"
precsn = "float64"
arch = "mace_mp"
model = "MACE-matpes-r2scan-omat-ft.model"
atoms.calc = choose_calculator(architecture=arch, model=model, precision=precsn, device=device)

mask=[1,1,1,1,1,1]

flag = LBFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=1.0e-3, steps=1000)

assert flag == True

energy_ase = atoms.get_potential_energy()

energy_pymc = pymc_energy()
print(energy_pymc)

assert energy_ase == pytest.approx(energy_pymc)

