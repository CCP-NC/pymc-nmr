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

class AirssStyle:
    def __init__(self):

        self.min_cfg:Config = None  # the lowest energy configuration
        self.sucessfull_airss = 0
        self.attempted_airss = 0
        


    def _setup_airss(self, spec: Species, job: JobControl, out_stream):
        numSpec = spec.get_num_species()

        
        # Swap of atom positions
        if job.num_swap_atoms > 0:
            self.numSwaps = job.num_swap_atoms
            for j in range(self.numSwaps):
                found = False
                ele1 = None
                ele2 = None

                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.swapType1[j]:
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
                        found = True
                        ele2 = ele
                        break

                if not found:
                    out_stream.write(f"\n atom swap type 2 {job.swapType2[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

    def initialise(self, spec, job, out_stream):
        
        #setup the BasinHopWalker calculation and the moves
        self._setup_airss(spec, job, out_stream)

    
    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin: Config, numSteps, cycle, 
            initialise, restart_iteration, restart_energy, out_stream):

        totalEnergy = Energy()  #stores the lowest energy so far
        checkEnergy = Energy()

        # Initiate the statistics
        stats.zero(1000, 0.0, False)
        type_stats.zero_types(1000, spec.get_num_species(), False) 
        
        fld.setup()

        if job.restart == False: #if there is restart then randomise the initial config
            old_pos = basin.get_positions()
            out_stream.write("\n creating a random configuration at start")
            for i in range(job.num_swap_atoms):
                basin.randomise(job.swapType1[i], job.swapType2[i], old_pos)

        new_basin = basin.create_atoms_object()
        totalEnergy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle)
        basin.update_from_atoms(new_basin)
           
        totalEnergy.print_energy(1, out_stream)

        
        numSteps = 1
        if job.restart:
            numSteps = restart_iteration

        while numSteps <= job.mcSteps:

            
            self.randomise_relax(basin, fld, totalEnergy, job, out_stream)

            new_basin = basin.create_atoms_object()
            energy_new = fld.calculate_energy(new_basin)
            basin.update_from_atoms(new_basin)
            
            print("energy in main routine ", energy_new.totalEnergy)
     
            if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                archive_io = open("archive.xyz", "a")
                basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                archive_io.close()
            
            stats.sample(job.equilSteps, numSteps, totalEnergy, basin.get_volume(), basin.cell_properties(), out_stream)
            type_stats.sample_types(numSteps, job.equilSteps, basin, spec)

            if numSteps % job.printFreq == 0:
                stats.check_point(numSteps, job.equilSteps, 0.0, out_stream)
                type_stats.check_point_types(spec, out_stream)

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
                        
                totalEnergy = checkEnergy

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
       
        if self.sucessfull_airss > 0:
            ratio = self.sucessfull_airss / self.attempted_airss
            out_stream.write(f"\n swaps : attempted, successful and ratio {self.attempted_airss} {self.sucessfull_airss} {ratio:.10e}\n")

        restart_io = open("restart.xyz", "w")
        basin.write_config(restart_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
        restart_io.close()
        
        out_stream.flush()


    def randomise_relax(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, out_stream):
        
        self.attempted_airss += 1

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy
        
        j = int(np.random.random() * float(len(job.swapType1)))

        old_pos = basin.get_positions()  # make a temprary config

        basin.randomise(job.swapType1[j], job.swapType2[j], old_pos)

        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        #write(filename="shuffled.xyz", images=new_basin, format="extxyz", append=False)
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle)
        
        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
    
        print("random ", old_energy.get_total_energy(), new_energy.get_total_energy()," ", deltaV)
    
        
        if deltaV < 0.0:
            print("airss new basin found")
            totalEnergy.totalEnergy = new_energy.totalEnergy
            self.sucessfull_airss += 1
            if job.save_downhill:
                 write(filename="downhill.xyz", images=new_basin, format="extxyz", append=True)
            
            basin.update_from_atoms(new_basin) # update the saved basin

        else:
            basin.set_positions(old_pos)
        