import numpy as np
import math
from typing import List

from ase import Atoms
from ase.io import write

from energy import Energy
from field import Field
from species import Species
from config import Config
from job_control import JobControl
from statistics import Statistics, TypeStatistics
from grid import Grid

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

        #combined atom swaps
        self.numCombiSwaps = 0
        self.successfulDownCombiSwaps = None
        self.attemptedCombiSwaps = None
        self.successfulUpCombiSwaps = None
        self.combi_swapType1 = []
        self.combi_swapType2 = []
        self.combi_swapType3 = []
        self.grd = None

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
        """
        sets up and checks the different types of swap
        
        Parameters
        ----------

        spec : Species oblect
            elemental data
        
        job : JobControl object
            control parameters

        out_stream : IO
            output file mc.log stream

        """
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
        #combined atom swaps
        if job.num_combi_swap_atoms > 0:
            self.numCombiSwaps = job.num_combi_swap_atoms
            self.successfulDownCombiSwaps = np.zeros(self.numCombiSwaps)
            self.attemptedCombiSwaps = np.zeros(self.numCombiSwaps)
            self.successfulUpCombiSwaps = np.zeros(self.numCombiSwaps)
            for j in range(self.numCombiSwaps):
                found = False

                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.combi_swapType1[j]:
                        self.combi_swapType1.append(job.combi_swapType1[j])
                        found = True
                        break

                if not found:
                    out_stream.write(f"\n combined atom swap type 1 {job.combi_swapType1[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

                found = False
                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.combi_swapType2[j]:
                        self.combi_swapType2.append(job.combi_swapType2[j])
                        found = True
                        break

                if not found:
                    out_stream.write(f"\n combined atom swap type 2 {job.combi_swapType2[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

                found = False
                for i in range(numSpec):
                    ele = spec.get_species(i)
                    if ele.name == job.combi_swapType3[j]:
                        self.combi_swapType3.append(job.combi_swapType3[j])
                        found = True
                        break

                if not found:
                    out_stream.write(f"\n combined atom swap type 3 {job.combi_swapType3[j]} not found in species list!\n")
                    out_stream.flush()
                    exit(EXIT_FAILURE)

            #initialise the grid
            self.grd = Grid(job.gridx, job.gridy, job.gridz, job.grid_cut)

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

        """
        The MC moves are created with the appropriate probability
        
        Parameters
        ----------

        job : JobControl object
            control parameters

        """
        self.numMCMoves = 0

        self.numMCMoves = (job.mdMoveFreq + job.swapFrequency + job.combi_swapFrequency)

        # Allocate mcMoveList array
        self.mcMoveList = np.zeros(self.numMCMoves, np.dtype('uint32'))

        j = 0
        for i in range(job.mdMoveFreq):
            self.mcMoveList[j] = 1
            j += 1

        for i in range(job.swapFrequency):
            self.mcMoveList[j] = 2
            j += 1

        for i in range(job.combi_swapFrequency):
            self.mcMoveList[j] = 3
            j += 1

    def initialise(self, spec, job, out_stream):
        """
        initialises the Basin Hopping functionality
        
        Parameters
        ----------

        spec : Species oblect
            elemental data
        
        job : JobControl object
            control parameters

        out_stream : IO
            output file mc.log stream

        """
        self._setupBasinHop(spec, job, out_stream)

        self._createMCMoves(job)

    def write_statistics(self, numSteps, total_energy, cell_properties, stats_io):
        stats_io.write(f" {numSteps} ")
        stats_io.write(f" {total_energy } ")
        stats_io.write(f" {cell_properties[0]} ")
        stats_io.write(f" {cell_properties[1]} ")
        stats_io.write(f" {cell_properties[2]} ")
        stats_io.write(f" {cell_properties[3]} ")
        stats_io.write(f" {cell_properties[4]} ")
        stats_io.write(f" {cell_properties[5]} \n")
        stats_io.flush()

    
    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin: Config, numSteps, cycle, 
            initialise, restart_iteration, restart_energy, out_stream):
        """
        The function that controls the Basin Hopping functionality and the flow of data output

        Parameters
        ----------

        spec : Species oblect
            elemental data
        
        fld : Field
            MLIP parameters and calculation of energies

        job : JobControl object
            control parameters

        stats : Statistics object
            maintains data on energies and unit cell

        type_stats: TypeStatistics
            maintains data on the number of types in the unit cell
        basin: Config
            initial configuration
             
        numSteps : int
            the number of steps within the current simulation
            
        cycle : int
             
        initialise : bool
            flag to indicate whether any initialisation is required
            
        restart_iteration : int
            the iteration read in from the configuration file
        
        restart_energy : float
            the energy read in from the configuration file
            
        out_stream : IO
            output file mc.log stream
        """

        totalEnergy = Energy()
        checkEnergy = Energy()

        # Initiate the statistics
        stats.zero(1000, 0.0, False)
        type_stats.zero_types(1000, spec.get_num_species(), False) 

        beta = 1.0 / (job.temperature * BOLTZMANN)
    
        out_stream.write(f"\n beta (1/KT) {beta:.8f}\n")
        
        fld.setup()

        new_basin = basin.create_atoms_object()
        totalEnergy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle, job.wrap)
        basin.update_from_atoms(new_basin)
           
        totalEnergy.print_energy(1, out_stream)

        numSteps = 1
        if job.restart:
            numSteps = restart_iteration

        if job.restart == False:
            archive_io = open("archive.xyz", "w")
            basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=0)
            archive_io.close()

            if job.writestats:
                stats_io = open("stats", "w")
                self.write_statistics(numSteps, totalEnergy.get_total_energy(), basin.cell_properties(), stats_io)
                stats_io.close()
        
        while numSteps <= job.mcSteps:

            choice = int(self.numMCMoves * np.random.random())

            selection = self.mcMoveList[choice]

            if selection == 1:
                self.run_md(basin, fld, totalEnergy, job, beta, out_stream)

            elif selection == 2:
                self.swapAtoms_relax(basin, fld, totalEnergy, job, beta, out_stream)

            elif selection == 3:
                self.grd.build_grid(basin) # before overey combi-swap is possibly over-kill
                                           #but I do not know how much they will move on relaxation
                self.combi_atom_swap_relax(basin, fld, totalEnergy, job, beta, out_stream)

            #new_basin = basin.create_atoms_object()
            #energy_new = fld.calculate_energy(new_basin)
            #basin.update_from_atoms(new_basin)
            
            #print("energy in main routine ", energy_new.totalEnergy)

            stats.sample(job.equilSteps, numSteps, totalEnergy, basin.get_volume(), basin.cell_properties(), out_stream)
            type_stats.sample_types(numSteps, job.equilSteps, basin, spec)

            if numSteps % job.printFreq == 0:
                stats.check_point(numSteps, job.equilSteps, 0.0, out_stream)
                type_stats.check_point_types(spec, out_stream)
                    
            if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                archive_io = open("archive.xyz", "a")
                basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                archive_io.close()

            if job.writestats and numSteps % job.writestats_freq == 0:
                stats_io = open("stats", "a")
                self.write_statistics(numSteps, totalEnergy.get_total_energy(), basin.cell_properties(), stats_io)
                stats_io.close()

            if numSteps % job.sanityCheckFreq == 0: 
                restart_io = open("restart.xyz", "w")
                basin.write_config(restart_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                restart_io.close()
              
                new_basin = basin.create_atoms_object()
                checkEnergy = fld.calculate_energy(new_basin, job.wrap)
                basin.update_from_atoms(new_basin)

                eDiff = checkEnergy.get_total_energy() - totalEnergy.get_total_energy()

                if abs(eDiff) > 1.0e-6:
                    out_stream.write(f"\n sanity check failed on iteration {numSteps} !!!!!!!\n")
                    out_stream.write(f" total diff {eDiff:.10e}\n")
                        
                totalEnergy.totalEnergy = checkEnergy.totalEnergy

            numSteps += 1

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(f" final energy of system containing {basin.natoms} atoms\n")
        out_stream.write(" *****************************************************************************************************\n")

        
        final_energy = Energy()
        new_basin = basin.create_atoms_object()
        final_energy = fld.calculate_energy(new_basin, job.wrap)
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

        for j in range(self.numCombiSwaps):
            swapRatio = (self.successfulUpCombiSwaps[j] + self.successfulDownCombiSwaps[j]) / self.attemptedCombiSwaps[j]
            successfulSwaps = self.successfulUpCombiSwaps[j] + self.successfulDownCombiSwaps[j]
            out_stream.write(f"\n combi-swaps {job.combi_swapType1[j]} {job.combi_swapType2[j]} {job.combi_swapType3[j]}: attempted, successful and ratio {self.attemptedCombiSwaps[j]} {successfulSwaps} {swapRatio:.10e}\n")
            out_stream.write(f"\n combi-swaps {job.combi_swapType1[j]} {job.combi_swapType2[j]} {job.combi_swapType3[j]}: Downhill and Uphill {self.successfulDownCombiSwaps[j]} {self.successfulUpCombiSwaps[j]} \n")

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
        """
        Swaps two atoms and then uses energy minimisation. The old and new energies are used within Boltzmann sampling

        Parameters
        ----------

        basin : Config
            The working configuration of atoms

        fld : Field
            The container for the energy minimisation using ASE calculator

        totalEnergy : Energy
            The working energy of the cell

        job : JobControl
            control parameters

        beta : float
            1 / kT

        out_stream : IO
            stream for writing out data

        """
        

        j = int(np.random.random() * self.numSwaps)
        #print("swap selection ",j," ", self.numSwaps,self.swapType1[j],self.swapType2[j])
        self.attemptedSwaps[j] += 1

        atm1 = basin.select_atom_of_type(self.swapType1[j])
        atm2 = basin.select_atom_of_type(self.swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy
        
        #print("swapping ", atm1, basin.symbol[atm1], atm2, basin.symbol[atm2])
        #print("pos atm1 ", basin.symbol[atm1], basin.pos[atm1,:])
        #print("pos atm2 ", basin.symbol[atm2], basin.pos[atm2,:])
        basin.swap_atom_positions(atm1, atm2)

        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle, job.wrap)
        
        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
        deltaVB = beta * deltaV
        #print("swap ", old_energy.get_total_energy(), new_energy.get_total_energy()," ", deltaV," ", beta, " ", deltaVB)
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
            #print("swap accepted")
            basin.update_from_atoms(new_basin)
        else:
            basin.swap_atom_positions(atm1, atm2)

    def combi_atom_swap_relax(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        """
        Swaps two atoms and then uses energy minimisation. A third atom is also moved to a suitable space in the structure/zeolite.
        The old and new energies are used within Boltzmann sampling - it does NOT obey detailed balance

        Parameters
        ----------

        basin : Config
            The working configuration of atoms

        fld : Field
            The container for the energy minimisation using ASE calculator

        totalEnergy : Energy
            The working energy of the cell

        job : JobControl
            control parameters

        beta : float
            1 / kT

        out_stream : IO
            stream for writing out data

        """
        j = int(np.random.random() * self.numCombiSwaps)
        #print("swap selection ",j," ", self.numSwaps,self.combi_swapType1[j],self.combi_swapType2[j], self.combi_swapType3)
        self.attemptedCombiSwaps[j] += 1

        atm1 = basin.select_atom_of_type(self.combi_swapType1[j])
        atm2 = basin.select_atom_of_type(self.combi_swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return
        
        atm3 = basin.find_closest_atom(atm2, self.combi_swapType3[j])
        #print("original atm2 pos",basin.symbol[atm2], basin.pos[atm2,0], basin.pos[atm2,1], basin.pos[atm2,2])
        #print("closest atm3 pos",basin.symbol[atm3], basin.pos[atm3,0], basin.pos[atm3,1], basin.pos[atm3,2])
        old_pos = basin.get_positions() # so we know where to put it back

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy
        
        #print("swapping ", atm1, basin.symbol[atm1], atm2, basin.symbol[atm2])
        #print("pos atm1 ", basin.symbol[atm1], basin.pos[atm1,:])
        #print("pos atm2 ", basin.symbol[atm2], basin.pos[atm2,:])
        self.grid_swap_atom_positions(basin, atm1, atm2, atm3, job.combi_swap_dist)
        
        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle, job.wrap)
        
        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
        deltaVB = beta * deltaV
        #print("swap ", old_energy.get_total_energy(), new_energy.get_total_energy()," ", deltaV," ", beta, " ", deltaVB)
        accept = False
        arg = np.random.random()
        if arg < np.exp(-deltaVB):
            accept = True
        
        if accept:
            totalEnergy.totalEnergy = new_energy.totalEnergy
            write(filename="accepted.xyz", images=new_basin, format="extxyz", append=True)
            if deltaV < 0.0:
                self.successfulDownCombiSwaps[j] += 1
                if job.save_downhill:
                    write(filename="downhill.xyz", images=new_basin, format="extxyz", append=True)
            else:
                self.successfulUpCombiSwaps[j] += 1
            #print("combi swap accepted")
            basin.update_from_atoms(new_basin)
        else:
            #basin.restore_grid_swap(atm1, atm2)
            basin.set_positions(old_pos)
        
    def run_md(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        """
        MD is used to shake the atoms into new configuration and then uses energy minimisation. 

        Parameters
        ----------

        basin : Config
            The working configuration of atoms

        fld : Field
            The container for the energy minimisation using ASE calculator

        totalEnergy : Energy
            The working energy of the cell

        job : JobControl
            control parameters

        beta : float
            1 / kT

        out_stream : IO
            stream for writing out data

        """
        
        self.md_runs += 1

        old_pos = basin.get_positions() 
        
        new_basin = basin.create_atoms_object()

        #run an md simulation
        fld.run_md(new_basin, job.mdtimestep, job.mdtemperature_K, job.mdfriction, job.mdsteps, job.wrap)

        new_energy = Energy()
        #the energy needs to be relaxed to get the "new" energy
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol, job.relstyle, job.wrap)
        
        if new_energy.totalEnergy > 1.0e5:   # ie it has failed to minimisa - this prevents it going into a stupid position

            basin.set_positions(old_pos)
        else:
            totalEnergy.totalEnergy = new_energy.totalEnergy

            basin.update_from_atoms(new_basin)

    def grid_swap_atom_positions(self, basin: Config, atm1, atm2, atm3, rcut):
        """
        Swaps the two atoms and then uses a cavity bias style approach to find a position for the thirs atom
        
        Parameters
        ----------

        basin : Config
            The working configuration of atoms
        
        atm1 : int
            index of the first atom

        atm2 : int
            index of the second atom

        atm3 : int
            index of the third atom where a space is required
        
        rcut : float
            the maximum distance for the grid point away from the second atom

        """

        basin.swap_atom_positions(atm1, atm2)

        #get grid positions near the new position of atom 2
        grd_list = self.grd.find_empty_grids(basin.pos[atm2,0], basin.pos[atm2,1], basin.pos[atm2,2], basin.get_vectors(), rcut)

        if len(grd_list) == 0:
            print("the value of combination swap distance is too small")
            exit()

        choice = int(len(grd_list) * np.random.random())
        atm = grd_list[choice]
        #print("grid choice", atm)
        basin.pos[atm3,:] = self.grd.grid_pos[atm,:]

        #print("atm2 pos",basin.symbol[atm2], basin.pos[atm2,0], basin.pos[atm2,1], basin.pos[atm2,2])
        #print("atm2 pos",basin.symbol[atm3], basin.pos[atm3,0], basin.pos[atm3,1], basin.pos[atm3,2])

        