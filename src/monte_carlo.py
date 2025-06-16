import numpy as np
import math
from typing import List

from ase import Atoms

from Energy import Energy
from Field import Field
from Species import Species
from JobControl import JobControl
from Statistics import Statistics, TypeStatistics
from config import Config

BOLTZMANN = .00008617333262145 # in eV
EXIT_FAILURE = 1

def _scale_displacement(rcp_vector, lat_vector, delta):

    xx = delta[0] * rcp_vector[0] + delta[1] * rcp_vector[3] + delta[2] * rcp_vector[6]
    yy = delta[0] * rcp_vector[1] + delta[1] * rcp_vector[4] + delta[2] * rcp_vector[7]
    zz = delta[0] * rcp_vector[2] + delta[1] * rcp_vector[5] + delta[2] * rcp_vector[8]

    delta[0] = xx * lat_vector[0] + yy * lat_vector[3] + zz * lat_vector[6];
    delta[1] = xx * lat_vector[1] + yy * lat_vector[4] + zz * lat_vector[7];
    delta[2] = xx * lat_vector[2] + yy * lat_vector[5] + zz * lat_vector[8]

class MonteCarlo:
    def __init__(self):

        self.numMCMoves = 0

        self.mcMoveList = []
        
        #displacement of atoms
        self.totalAtomMoves = 0
        self.successfulAtomMoves = 0
        self.distance_atom_max = None
        self.ratioAtomMove = None
        self.attemptedAtomMoves = None
        self.noAtomMoves = None
        self.noAtomMovers = None
        self.atomMoveTypes = None

        #volume moves
        self.maxVolChange = None
        self.successful_vol_change = None
        self.numVolChange = None
        self.attemptedVolChange = None
        self.totalVolChanges = None

        #atom swaps
        self.numSwaps = 0
        self.successfulSwaps = 0
        self.attemptedSwaps = 0
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

    
        

        #fragments
        self.noFragments = 0

    def _setupMonteCarlo(self, spec: Species, job: JobControl, out_stream):
        numSpec = spec.get_num_species()

        # Allocate distance and ratio vectors
        self.distance_atom_max = np.zeros(numSpec, dtype=np.float64) 
        self.distance_atom_max[:] = job.maxDistance
        self.ratioAtomMove = np.zeros(numSpec, dtype=np.float64) 
        self.attemptedAtomMoves = np.zeros(numSpec, dtype=np.float64) 
        self.noAtomMoves = np.zeros(numSpec, dtype=np.float64) 
        self.totalAtomMoves = 0
        self.successfulAtomMoves = 0

        self.noAtomMovers = len(job.moveTypes)
        self.atomMoveTypes = np.zeros(self.noAtomMovers, np.dtype('uint32'))

        # Initialize atom move types
        for j in range(self.noAtomMovers):
            found = False
            for i in range(numSpec):
                ele = spec.get_species(i)
                if ele.name == job.moveTypes[j]:
                    self.atomMoveTypes[j] = i
                    found = True
                    break

            if not found:
                out_stream.write(f"\n move atom type {job.moveTypes[j]} not found in species list!\n")
                out_stream.flush()
                exit(EXIT_FAILURE)

        # Setup volume moves
        self.maxVolChange = np.zeros(6, dtype=np.float64)
        self.maxVolChange[:] = job.maxVolDisplacement
        self.successful_vol_change = np.zeros(6, dtype=np.float64)
        self.numVolChange = np.zeros(6, dtype=np.float64)
        self.attemptedVolChange = np.zeros(6, dtype=np.float64)
        self.totalVolChanges = np.zeros(6, dtype=np.float64)

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

        self.numMCMoves = (job.atomMoveFreq + job.volMoveFreq + job.swapFrequency + job.transmutateFrequency)

        # Allocate mcMoveList array
        self.mcMoveList = np.zeros(self.numMCMoves, np.dtype('uint32'))

        j = 0
        for i in range(job.atomMoveFreq):
            self.mcMoveList[j] = 1
            j += 1

        for i in range(job.swapFrequency):
            self.mcMoveList[j] = 2
            j += 1

        for i in range(job.volMoveFreq):
            self.mcMoveList[j] = 3
            j += 1


    def initialise(self, spec, job, out_stream):
        #setup the MonteCarloWalker calculation and the moves
        self._setupMonteCarlo(spec, job, out_stream)
    
        self._createMCMoves(job)
        


    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin: Config, numSteps, cycle, initialise, out_stream):

       
        volRatio = np.zeros(6, dtype=np.float64)

        totalEnergy = Energy()
        checkEnergy = Energy()

        # Initiate the statistics
        stats.zero(1000, 0.0, False)
        type_stats.zero_types(1000, spec.get_num_species(), False) 

        beta = 1.0 / (job.temperature * BOLTZMANN)
    
        out_stream.write(f"\n beta (1/KT) {beta:.8f}\n")
        
        fld.setup()

        new_basin = basin.create_atoms_object()
        totalEnergy = fld.calculate_energy(new_basin)
           
        #totalEnergy.totalEnergy = energy_new.totalEnergy
        totalEnergy.print_energy(1, out_stream)

        if job.restart == False:
            archive_io = open("archive.xyz", "w")
            basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=0)
            archive_io.close()

        numSteps = 1

        while numSteps <= job.mcSteps:

            choice = int(self.numMCMoves * np.random.random())

            selection = self.mcMoveList[choice]

            if selection == 1:
                self.move_atom(basin, fld, totalEnergy, beta, out_stream)

            elif selection == 2:
                self.swapAtoms(basin, fld, totalEnergy, job, beta, out_stream)

            elif selection == 3:               
                 self.move_volume(basin, fld, totalEnergy, spec, job, beta, out_stream)

            stats.sample(job.equilSteps, numSteps, totalEnergy, basin.get_volume(), basin.cell_properties(), out_stream)
            type_stats.sample_types(numSteps, job.equilSteps, basin, spec)

            if numSteps % job.printFreq == 0:
                stats.check_point(numSteps, job.equilSteps, 0.0, out_stream)
                type_stats.check_point_types(spec, out_stream)
            
            if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                archive_io = open("archive.xyz", "a")
                basin.write_config(archive_io, total_energy=totalEnergy.get_total_energy(), iteration=numSteps)
                archive_io.close()

            
            if numSteps % job.accAtomMoveUpdate == 0 and numSteps > 0:
                for i in range(spec.number_of_elements):
                    if self.attemptedAtomMoves[i] > 0:
                        self.ratioAtomMove[i] = self.noAtomMoves[i] / self.attemptedAtomMoves[i]

                        if self.ratioAtomMove[i] > job.acceptAtomMoveRatio:
                            self.distance_atom_max[i] *= 1.05
                        else:
                            self.distance_atom_max[i] *= 0.95

                    self.noAtomMoves[i] = 0
                    self.attemptedAtomMoves[i] = 0

            if numSteps % job.acceptVolUpdate == 0 and numSteps > 0:
                for i in range(6):
                    if self.attemptedVolChange[i] > 0:
                        volRatio[i] = self.numVolChange[i] / self.attemptedVolChange[i]

                        if volRatio[i] > job.acceptVolMoveRatio:
                            self.maxVolChange[i] *= 1.05
                        else:
                            self.maxVolChange[i] *= 0.95

                    self.numVolChange[i] = 0.0
                    self.attemptedVolChange[i] = 0

            #if numSteps > job.equilSteps and job.sampleBasin and numSteps % job.sampleBasinFreq == 0:
            #    basin[cycle].samplePositions()

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
                else:
                    out_stream.write(f"\n sanity check passed  on iteration {numSteps} \n")
                        
                totalEnergy.totalEnergy = checkEnergy.totalEnergy


            numSteps += 1

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(f" final energy of system containing {basin.get_number_of_atoms()} atoms\n")
        out_stream.write(" *****************************************************************************************************\n")

        
        final_energy = Energy()
        new_basin = basin.create_atoms_object()
        final_energy = fld.calculate_energy(new_basin)

        final_energy.print_energy(0, out_stream)

        out_stream.write("\n final sanity check")
        checkEnergy = final_energy - totalEnergy
        checkEnergy.print_energy(cycle, out_stream)

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(" Summary of simulation\n")
        out_stream.write(" *****************************************************************************************************\n")

        out_stream.write("\n")
        if self.totalAtomMoves > 0:
            out_stream.write(f" total, successful atom moves and ratio {self.totalAtomMoves} {self.successfulAtomMoves} "
                            f"{self.successfulAtomMoves / self.totalAtomMoves:.10e}\n")
            out_stream.write("\n")
            for i in range(spec.number_of_elements):
                out_stream.write(f" maximum atom displacement {i} {self.distance_atom_max[i]}\n")

        out_stream.write("\n")
        for indx in range(6):
            if self.successful_vol_change[indx] > 0:
                out_stream.write(f" total, successful volume changes and ratio {self.totalVolChanges[indx]} "
                                f"{self.successful_vol_change[indx]} {self.successful_vol_change[indx] / self.totalVolChanges[indx]:.10e}\n")

        for indx in range(6):
            if self.successful_vol_change[indx] > 0:
                out_stream.write(f" maximum volume change parameter {self.maxVolChange[indx]}\n")

        if self.numSwaps > 0:
            swapRatio = self.successfulSwaps / self.attemptedSwaps
            out_stream.write(f"\n swaps : attempted, successful and ratio {self.attemptedSwaps} {self.successfulSwaps} {swapRatio:.10e}\n")


        out_stream.flush()

    

    def move_atom(self, basin: Config, fld: Field, total_energy: Energy, beta: float, out_stream):
        
        
        choice = int(np.random.random() * self.noAtomMovers)
        typ = self.atomMoveTypes[choice]
        
        atm = basin.select_atom()
        
        if atm < 0:
            return
        
        energy_old = Energy()
        energy_old.totalEnergy = total_energy.totalEnergy
        energy_new = Energy()
        
        
        
        self.totalAtomMoves += 1
        self.attemptedAtomMoves[typ] += 1

        print("before move ", atm, basin.pos[atm,:])
        old_pos = basin.make_atom_move(atm, self.distance_atom_max[typ])
        print("after move ", atm, basin.pos[atm,:])
        
        #create atom object and calculate new energy
        new_basin = basin.create_atoms_object()
        energy_new = fld.calculate_energy(new_basin)
        
        deltaV = energy_new.get_total_energy() - energy_old.get_total_energy()
        deltaVB = deltaV * beta
        print("smove ", energy_old.get_total_energy(), energy_new.get_total_energy(), deltaV, deltaVB)
        #energyDifference.print_energy(0, out_stream)
        accept = False
        arg = np.random.random()

        if deltaV < 0.0:
            accept = True
        else:
            try:
                if arg < np.exp(-deltaVB):
                    accept = True
            except Exception as e:
                out_stream.write(f"{e} whilst moving atom {atm} \n")
                accept = False
            
        if accept:
            #update the total energy (basin can remain the same)
            total_energy.totalEnergy = energy_new.totalEnergy
            
            self.noAtomMoves[typ] += 1
            self.successfulAtomMoves += 1
            print("accepted ")
            
        else:
            #revert basin back to its old state
            basin.reject_atom_move(atm, old_pos)
            print("rejected")


    def move_volume(self, basin: Config, fld: Field, totalEnergy: Energy, spec: Species, job: JobControl, beta: float, out_stream):
        indx = 0
        vol_new = 1.0
        natoms = basin.natoms

        betaInv = 1.0 / beta

        oldEnergy = Energy()
        oldEnergy.totalEnergy = totalEnergy.totalEnergy
        vol_old = basin.get_volume()
        old_vec = basin.get_vectors()
        old_pos = basin.get_positions()

        bulks = np.ones(3, dtype=np.float64)
        
        if job.volMoveSymmetry == 0:
            indx = 0
            maxVol = self.maxVolChange[0]
            vol_new = basin.expand_cell_cubic(bulks, maxVol)
        elif job.volMoveSymmetry == 1:
            indx = int(2.0 * np.random.random())
            vol_new = basin.expand_cell_tetragonal(indx, bulks, self.maxVolChange)
        elif job.volMoveSymmetry == 2:
            indx = int(3.0 * np.random.random())
            vol_new = basin.expand_cell_orthorhombic(indx, bulks, self.maxVolChange)
        else:
            pass  # Implement other volume changes if needed
        
        self.attemptedVolChange[indx] += 1
        self.totalVolChanges[indx] += 1

        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        new_energy = fld.calculate_energy(new_basin)

        deltav = new_energy.get_total_energy() - oldEnergy.get_total_energy()
        #print("energies", new_energy.get_total_energy(), oldEnergy.get_total_energy(), deltav)
        arg = beta * (deltav + job.extPressure * (vol_new - vol_old) - (natoms) * betaInv * math.log(vol_new / vol_old))
        #print ("old ", oldEnergy.get_total_energy(), " new ", new_energy.get_total_energy(), " diff ", deltav, "arg ", math.exp(-arg))
        rNum = np.random.random()
        
        if rNum < math.exp(-arg):
            totalEnergy.totalEnergy = new_energy.totalEnergy
            
            self.numVolChange[indx] += 1
            self.successful_vol_change[indx] += 1

        else:
            basin.set_positions(old_pos)
            basin.set_vectors(old_vec)

    def swapAtoms(self, basin: Config, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        self.attemptedSwaps += 1

        j = int(np.random.random() * self.numSwaps)
        
        atm1 = basin.select_atom_of_type(self.swapType1[j])
        atm2 = basin.select_atom_of_type(self.swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy

        #print("swapping", atm1, atm2, basin.chem_symbols[atm1], basin.chem_symbols[atm2])
        basin.swap_atom_positions(atm1, atm2)

        new_energy = Energy()
        new_basin = basin.create_atoms_object()
        new_energy = fld.calculate_energy(new_basin)

        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
        deltaVB = beta * deltaV
        print("swap ", old_energy.get_total_energy(), new_energy.get_total_energy(), deltaV, deltaVB)
        accept = False
        arg = np.random.random()
        if deltaV < 0.0:
            accept = True
        else:
            try:
                if arg < np.exp(-deltaVB):
                    accept = True
            except Exception as e:
                out_stream.write(f"{e} whilst swapping atoms {atm1} and {atm2}\n")
                accept = False
            
        if accept:
            totalEnergy.totalEnergy = new_energy.totalEnergy
            self.successfulSwaps += 1
            #basin.update_from_atoms(new_basin)

        else:
            basin.swap_atom_positions(atm1, atm2)

    
 

    