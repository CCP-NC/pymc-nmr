import sys

import numpy as np

from ase.filters import UnitCellFilter
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.md import VelocityVerlet, langevin
from ase.units import fs

from janus_core.helpers.mlip_calculators import choose_calculator

from Energy import Energy
from Species import Species

class Field:

    def __init__(self):
        
        self.model = None # full path to the MLIP parameters
        self.species = None # list of elements

        self.device = "cuda" # device where the calculations for MLIP are undertaken
        self.janCalc = None # Calculator

        self.arch = "mace_mp" # MLIP type or architecture
        self.precision = "float64" # precision used by MLIP

        self.first_setup = True # flag if multiple calls are carried out to prevent further setup+

    def readPotential(self, in_stream, out_stream, spec: Species):
        """
        This function read in the species and location of MLIP library

        Parameters
        ----------

        in_stream : io stream
            stream for reading data - should already be opened

        out_stream :  io stream 
            stream for writing - it is redundant in this instance but included for expansion etc
            
        spec : Species object
            list of elements

        """

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

            elif words[0].lower() == "precision": 
                self.precision = words[1].lower()

            elif words[0].lower() == "arch":
                self.arch = words[1].lower()

            elif words[0].lower() == "model":
                self.model = words[1]

        
        
    def setup(self):
        """
        The function sets up the fiels and creates an ASE calculator for use later
        """
        if self.first_setup == False:
            return
        
        if self.model is None:
            print("a model name is required")
            exit(-1)

        if self.species == None:
            print("the spcies must be provided")
            exit(-1)
        
        try:
            self.janCalc = choose_calculator(architecture=self.arch, model=self.model, precision=self.precision, device=self.device)
        except Exception as e:
            print(f"{e} whilst trying to load {self.model}\n")
            exit()


        self.first_setup = False

    def calculate_energy(self, atoms: Atoms, wrap):
        """
        The subroutine is a wrapper around ASE get_potential_energy function

        Parameters
        ----------
        atoms : Atoms
            An ASE atoms object with positions etc

        wrap : bool
            if true the cell is reconfigured so that all atoms are within the box

        returns the total energy

        Returns
        -------

        total_energy : float
            the final energy in eV
        """
        
        total_energy = Energy()

        if wrap:
            atoms.wrap()
        
        atoms.calc = self.janCalc

        total_energy.totalEnergy = atoms.get_potential_energy()
        
        return total_energy

    def calculate_energy_relax(self, atoms: Atoms, relmethod, relsteps, reltol, relstyle, wrap):
        """
        The subroutine is a wrapper around ASE energy/force minimisation routines

        Parameters
        ----------
        atoms : Atoms
            An ASE atoms object with positions etc

        relmethod : str
            string indetifying the relaxation method

        relsteps : int
            the max number of steps

        reltol : float
            accuracy of the minimisation

        selstyle : str
            indentifies which cell parameters should be relaxed

        wrap : bool
            if true the cell is reconfigured so that all atoms are within the box

        Returns
        -------

        total_energy : float
            the final energy in eV on a valid minimisation else a very high energy is sent back
        """
        
        total_energy = Energy()

        if relstyle == "conv":
            mask=[0,0,0,0,0,0]
        elif relstyle == "cona":
            mask=[1,1,1,0,0,0]
        else:
            mask=[1,1,1,1,1,1]

        if wrap:
            atoms.wrap()
       
        atoms.calc = self.janCalc
        
        if relmethod == "lbfgs":
            flag = LBFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=reltol, steps=relsteps)
        elif relmethod == "fire":
            flag = FIRE(UnitCellFilter(atoms, mask=mask)).run(fmax=reltol, steps=relsteps)
        elif relmethod == "bfgs":
            flag = BFGS(UnitCellFilter(atoms, mask=mask)).run(fmax=reltol, steps=relsteps)
        else:
            print("unrecognised relaxation method")
            exit()

        if flag:
            total_energy.totalEnergy = atoms.get_potential_energy()
            #print("field final energy ", atoms.get_potential_energy())
        else:
            total_energy.totalEnergy = 1.0e6
        
        return total_energy
    
    def run_md(self, atoms: Atoms, timestep, mdtemperature_K, mdfriction, mdsteps, wrap):
        """
        The subroutine is a wrapper around ASE constant volume MD. An MD of specified length
        is undertaken using the Langevin thormostat to shake the system.

        Parameters
        ----------
        atoms : Atoms
            An ASE atoms object with positions etc

        timestep : float
            the MD timestep in fs

        mdtemperature_K : float
            the temperature in Kelvin

        mdfriction : float
            friction parameter for the Langevin thermostat

        mdsteps : str
            the number of steps for the MD simulations

        wrap : bool
            if true the cell is reconfigured so that all atoms are within the box
        
        Returns
        -------

        total_energy : float
            the final energy in eV 
        """
        total_energy = Energy()

        if wrap:
            atoms.wrap()
        
        atoms.calc = self.janCalc

        #dyn = VelocityVerlet(atoms, timestep=2.0 * fs, temperature_K=1000)
        dyn = langevin.Langevin(atoms, timestep=timestep, temperature_K=mdtemperature_K, friction=mdfriction)
        dyn.run(mdsteps)

        total_energy.totalEnergy = atoms.get_potential_energy()
        
        return total_energy

    
