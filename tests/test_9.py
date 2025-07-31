#tests MC displacement of volume
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

    np.random.seed(0)
    
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
    old_energy = fld.calculate_energy(atoms, False)

    old_vec = basin.get_vectors()
    old_pos = basin.get_positions()
    bulks = np.zeros(6, dtype=np.float64)
    maxVol = np.zeros(6, dtype=np.float64)
    maxVol[:] = 0.1
    indx = int(6.0 * np.random.random())
    vol_new = basin.distort_cell(indx, bulks, maxVol)
    
    atoms = basin.create_atoms_object()
    new_energy = fld.calculate_energy(atoms, False)
    delta_v1 = old_energy.totalEnergy - new_energy.totalEnergy

    basin.set_positions(old_pos)
    basin.set_vectors(old_vec)
    atoms = basin.create_atoms_object()
    new_energy = fld.calculate_energy(atoms, False)
    delta_v2 = old_energy.totalEnergy - new_energy.totalEnergy

    return delta_v1, delta_v2


delta_v1, delta_v2 = pymc_energy()
print(delta_v1, delta_v2)

assert delta_v1 == pytest.approx(-0.5764823374752268)
assert delta_v2 == pytest.approx(0.0)

