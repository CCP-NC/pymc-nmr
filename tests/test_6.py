#determines energy of MC swap
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

    np.random.seed(0) # random seed to get same sequence

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
    old_energy = fld.calculate_energy_relax(atoms, "lbfgs", 1000, 1.0e-3, "conp", False)

    atm1 = basin.select_atom_of_type("Si")
    atm2 = basin.select_atom_of_type("Al")
    basin.swap_atom_positions(atm1, atm2)
    atoms = basin.create_atoms_object()
    new_energy = fld.calculate_energy_relax(atoms, "lbfgs", 1000, 1.0e-3, "conp", False)
    delta_v1 = old_energy.totalEnergy - new_energy.totalEnergy

    basin.swap_atom_positions(atm1, atm2) # swap atoms back around
    atoms = basin.create_atoms_object()
    new_energy = fld.calculate_energy_relax(atoms, "lbfgs", 1000, 1.0e-3, "conp", False)
    delta_v2 = old_energy.totalEnergy - new_energy.totalEnergy

    return delta_v1, delta_v2


delta_v1, delta_v2 = pymc_energy()
print(delta_v1, delta_v2)

assert delta_v1 == pytest.approx(-0.12121373476838926)
assert delta_v2 == pytest.approx(0.0, rel=1.0e-6, abs=1.0e-6)

