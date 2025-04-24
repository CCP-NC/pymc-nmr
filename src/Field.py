import sys

import numpy as np

from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from janus_core.helpers.mlip_calculators import choose_calculator

from Energy import Energy
from Species import Species

class Field:

    def __init__(self):
        
        self.model = None
        self.struc = None
        self.cutoff = 0.0
        self.species = None
        #self.atoms = None

        self.device = "cuda"

        #self.atoms = Atoms()
        self.janCalc = None #janCalculator()

        self.first_setup = True

    def readPotential(self, in_stream, out_stream, spec: Species):
        #read the species and potentials

        while True:
            line = in_stream.readline()
            if not line:
                break
            
            #out_stream.write(f" input line: {line.strip()}\n")
            if line[0] == '#':
                continue

            words = line.split()

            if not words:
                continue

            if words[0].lower() == "close":
                return
            
            elif words[0] == "species":
                num = int(words[1])
                spec.load_species(in_stream, num)

            elif words[0].lower() == "device":
                self.device = words[1].lower()

            elif words[0].lower() == "model":
                self.model = words[1]

        
        
    def setup(self):

        if self.first_setup == False:
            return
        
        if self.model is None:
            print("a model name is required")
            exit(-1)
        
        self.janCalc = choose_calculator(architecture="mace_mp", model=self.model, precision="float64", device=self.device)


        self.first_setup = False

    def calculate_energy(self, atoms: Atoms):
        
        total_energy = Energy()
        
        atoms.set_calculator(self.janCalc)

        total_energy.totalEnergy = atoms.get_potential_energy()
        
        return total_energy

    def calculate_energy_relax(self, atoms: Atoms):
        
        total_energy = Energy()
        max_force = 1e-3
        max_steps = 1000
       
        atoms.calc = self.janCalc

        
        flag = LBFGS(UnitCellFilter(atoms, mask=[1,1,1,1,1,1])).run(fmax=max_force, steps=max_steps)

        if flag:
            total_energy.totalEnergy = atoms.get_potential_energy()
            print("field final energy ", atoms.get_potential_energy())
        else:
            total_energy.totalEnergy = 1.0e6
        
        return total_energy

    
