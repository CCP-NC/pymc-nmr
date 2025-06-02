import numpy as np
import math
from typing import List

from ase import Atoms
from ase.io import write

from Energy import Energy
from Field import Field
from Species import Species
from config import Config
from JobControl import JobControl
from Statistics import Statistics, TypeStatistics

BOLTZMANN = .00008617333262145 # in eV
EXIT_FAILURE = 1

class BasinHop:
    def __init__(self):

        self.numMCMoves = 0

        self.mcMoveList = []
        

        #atom swaps
        self.numSwaps = 0
        self.successfulDownSwaps = None
        self.attemptedSwaps = None
        self.successfulUpSwaps = None
        self.swapType1 = []
        self.swapType2 = []

        #semi grand or transmutational
        self.numTrans = 0
        self.transmuteChemPot = 0.0
        self.transType1 = []
        self.transType2 = []
        self.forwardMutations = 0
        self.attemptForwardMutations = 0
        self.backwardMutations = 0
        self.attemptBackwardMutations = 0

        #semi-Widom method
        self.numSemiWidom = 0
        self.chemPotAtom1 = None
        self.chemPotAtom2 = None
        self.semiWidomType1 = None
        self.semiWidomType2 = None

        
        self.md_runs = 0

    def _setupBasinHop(self, spec: Species, job: JobControl, out_stream):
        numSpec = spec.get_num_species()

        
        # Swap of atom positions
        if job.num_swap_atoms > 0:
            self.numSwaps = job.num_swap_atoms
            self.successfulDownSwaps = np.zeros(self.numSwaps)
            self.attemptedSwaps = np.zeros(self.numSwaps)
            self.successfulUpSwaps = np.zeros(self.numSwaps)
            for j in range(self.numSwaps):
                found = False
                ele1 = None
                ele2 = None

                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.swapType1[j]:
                        self.swapType1.append(job.swapType1[j])
                        found = True
                        ele1 = ele
                        break

                if not found:
                    out_stream.write(f"\n atom swap type 1 {job.swapType1[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

                found = False
                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.swapType2[j]:
                        self.swapType2.append(job.swapType2[j])
                        found = True
                        ele2 = ele
                        break

                if not found:
                    out_stream.write(f"\n atom swap type 2 {job.swapType2[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

                if ele1.charge != ele2.charge:
                    self.chargedSwap = True

        # Transmutation of atom positions
        if job.num_transmutate_atoms > 0:
            self.numTrans = job.num_transmutate_atoms
            self.transmuteChemPot = job.transmuteChemPot

            for j in range(self.numTrans):
                found = False
                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.mutateType1[j]:
                        self.transType1.append(job.mutateType1[j])
                        found = True
                        break

                if not found:
                    out_stream.write(f"\n atom transmutation type 1 {job.mutateType1[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

                found = False
                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.mutateType2[j]:
                        self.transType2.append(job.mutateType2[j])
                        found = True
                        break

                if not found:
                    out_stream.write(f"\n atom transmutation type 2 {job.mutateType2[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

    def _createMCMoves(self, job):
        self.numMCMoves = 0

        self.numMCMoves = (job.mdMoveFreq + job.swapFrequency)

        # Allocate mcMoveList array
        self.mcMoveList = np.zeros(self.numMCMoves, np.dtype('uint32'))

        j = 0
        for i in range(job.mdMoveFreq):
            self.mcMoveList[j] = 1
            j += 1

        for i in range(job.swapFrequency):
            self.mcMoveList[j] = 2
            j += 1

    def initialise(self, spec, job, out_stream):
        #setup the BasinHopWalker calculation and the moves
        self._setupBasinHop(spec, job, out_stream)

        self._createMCMoves(job)
    
    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin: Config, numSteps, cycle, 
            initialise, restart_iteration, restart_energy, out_stream):

        totalEnergy = Energy()
        checkEnergy = Energy()

        # Initiate the statistics
        stats.zero(1000, 0.0, False)
        type_stats.zero_types(1000, spec.get_num_species(), False) 

        beta = 1.0 / (job.temperature * BOLTZMANN)
    
        out_stream.write(f"\n beta (1/KT) {beta:.8f}\n")
        
        fld.setup()

        new_basin = basin.create_atoms_object()
        totalEnergy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol)
        basin.update_from_atoms(new_basin)
           
        totalEnergy.print_energy(1, out_stream)

        
        if job.restart == False:
            archive_io = open("archive.xyz", "w")
            basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=0)
            archive_io.close()
        
        numSteps = 1
        if job.restart:
            numSteps = restart_iteration

        while numSteps <= job.mcSteps:

            choice = int(self.numMCMoves * np.random.random())

            selection = self.mcMoveList[choice]

            if selection == 1:
                self.run_md(basin, fld, totalEnergy, job, beta, out_stream)

            elif selection == 2:
                self.swapAtoms_relax(basin, fld, totalEnergy, job, beta, out_stream)

            new_basin = basin.create_atoms_object()
            energy_new = fld.calculate_energy(new_basin)
            basin.update_from_atoms(new_basin)
            
            print("energy in main routine ", energy_new.totalEnergy)

            stats.sample(job.equilSteps, numSteps, totalEnergy, basin.get_volume(), basin.vectors.flatten(), out_stream)
            type_stats.sample_types(numSteps, job.equilSteps, basin, spec)

            if numSteps % job.printFreq == 0:
                stats.check_point(numSteps, job.equilSteps, 0.0, out_stream)
                type_stats.check_point_types(spec, out_stream)
                    
            if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                archive_io = open("archive.xyz", "a")
                basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                archive_io.close()

            if numSteps % job.sanityCheckFreq == 0: 
                restart_io = open("restart.xyz", "w")
                basin.write_config(restart_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                restart_io.close()
              
                new_basin = basin.create_atoms_object()
                checkEnergy = fld.calculate_energy(new_basin)
                basin.update_from_atoms(new_basin)

                eDiff = checkEnergy.get_total_energy() - totalEnergy.get_total_energy()

                if abs(eDiff) > 1.0e-6:
                    out_stream.write(f"\n sanity check failed on iteration {numSteps} !!!!!!!\n")
                    out_stream.write(f" total diff {eDiff.totalEnergy:.10e}\n")
                        
                totalEnergy.totalEnergy = checkEnergy.totalEnergy

            numSteps += 1

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(f" final energy of system containing {basin.natoms} atoms\n")
        out_stream.write(" *****************************************************************************************************\n")

        
        final_energy = Energy()
        new_basin = basin.create_atoms_object()
        final_energy = fld.calculate_energy(new_basin)
        basin.update_from_atoms(new_basin)

        final_energy.print_energy(1, out_stream)

        out_stream.write("\n final sanity check")
        checkEnergy = final_energy - totalEnergy
        checkEnergy.print_energy(cycle, out_stream)

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(" Summary of simulation\n")
        out_stream.write(" *****************************************************************************************************\n")

        out_stream.write("\n")
       
        for j in range(self.numSwaps):
            swapRatio = (self.successfulUpSwaps[j] + self.successfulDownSwaps[j]) / self.attemptedSwaps[j]
            successfulSwaps = self.successfulUpSwaps[j] + self.successfulDownSwaps[j]
            out_stream.write(f"\n swaps {job.swapType1[j]} {job.swapType2[j]}: attempted, successful and ratio {self.attemptedSwaps[j]} {successfulSwaps} {swapRatio:.10e}\n")
            out_stream.write(f"\n swaps {job.swapType1[j]} {job.swapType2[j]}: Downhill and Uphill {self.successfulDownSwaps[j]} {self.successfulUpSwaps[j]} \n")

        if self.numTrans > 0:
            forwardRatio = self.forwardMutations / self.attemptForwardMutations
            backwardRatio = self.backwardMutations / self.attemptBackwardMutations
            out_stream.write(f"\n forward mutations : attempted, successful and ratio {self.attemptForwardMutations} "
                            f"{self.forwardMutations} {forwardRatio:.10e}\n")
            out_stream.write(f" backward mutations : attempted, successful and ratio {self.attemptBackwardMutations} "
                        f"{self.backwardMutations} {backwardRatio:.10e}\n")

        if self.numSemiWidom > 0:
            out_stream.write(f"\n forward semi-widom mutations {self.forwardSemiWidom}\n")
            out_stream.write(f" backward semi-widom mutations {self.backwardSemiWidom}\n")

        if self.md_runs > 0:
            out_stream.write(f"\n the number of MD runs {self.md_runs}\n")

        restart_io = open("restart.xyz", "w")
        basin.write_config(restart_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
        restart_io.close()
        
        out_stream.flush()


    def swapAtoms_relax(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        

        j = int(np.random.random() * self.numSwaps)
        print("swap selection ",j," ", self.numSwaps,self.swapType1[j],self.swapType2[j])
        self.attemptedSwaps[j] += 1

        atm1 = basin.select_atom_of_type(self.swapType1[j])
        atm2 = basin.select_atom_of_type(self.swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy
        

        #print("swapping", atm1, atm2, basin.chem_symbols[atm1], basin.chem_symbols[atm2])
        print("swapping ", atm1, basin.symbol[atm1], atm2, basin.symbol[atm2])
        print("pos atm1 ", basin.symbol[atm1], basin.pos[atm1,:])
        print("pos atm2 ", basin.symbol[atm2], basin.pos[atm2,:])
        basin.swap_atom_positions(atm1, atm2)

        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol)
        
        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
        deltaVB = beta * deltaV
        print("swap ", old_energy.get_total_energy(), new_energy.get_total_energy()," ", deltaV," ", beta, " ", deltaVB)
        accept = False
        arg = np.random.random()
        if arg < np.exp(-deltaVB):
            accept = True
        
        if accept:
            totalEnergy.totalEnergy = new_energy.totalEnergy
            write(filename="accepted.xyz", images=new_basin, format="extxyz", append=True)
            if deltaV < 0.0:
                self.successfulDownSwaps[j] += 1
                if job.save_downhill:
                    write(filename="downhill.xyz", images=new_basin, format="extxyz", append=True)
            else:
                self.successfulUpSwaps[j] += 1
            print("swap accepted")
            basin.update_from_atoms(new_basin)
        else:
            basin.swap_atom_positions(atm1, atm2)
        
    def run_md(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        self.md_runs += 1

        j = int(np.random.random() * self.numSwaps)
        
        new_basin = basin.create_atoms_object()

        #run an md simulation
        fld.run_md(new_basin, job.timestep, job.mdtemperature_K, job.mdfriction, job.mdsteps)

        new_energy = Energy()
        #the energy needs to be relaxed to get the "new" energy
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol)
        
        totalEnergy.totalEnergy = new_energy.totalEnergy

        basin.update_from_atoms(new_basin)
        