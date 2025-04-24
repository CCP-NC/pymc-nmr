import numpy as np
import math
from typing import List

from ase import Atoms

from Energy import Energy
from Field import Field
from Species import Species
from JobControl import JobControl
from Statistics import Statistics, TypeStatistics

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

        for i in range(job.transmutateFrequency):
            self.mcMoveList[j] = 3
            j += 1

        for i in range(job.volMoveFreq):
            self.mcMoveList[j] = 5
            j += 1


    def initialise(self, spec, job, out_stream):
        #setup the MonteCarloWalker calculation and the moves
        self._setupMonteCarlo(spec, job, out_stream)
    
        self._createMCMoves(job)
        


    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin, numSteps, cycle, initialise, out_stream):

       
        volRatio = np.zeros(6, dtype=np.float64)

        totalEnergy = []
        checkEnergy = Energy()

        # Initiate the statistics
        for i in range(job.num_boxes):
            stats[i].zero(1000, 0.0, False)
            type_stats[i].zero_types(spec.get_num_species(), False) 

        xyzStream = None

        if job.dumpArchive:
            out_stream.write(f"\n archive frequency {job.archiveFrequency}\n")
            xyzStream = open("mc_archive.xyz", "w")
            

        beta = 1.0 / (job.temperature * BOLTZMANN)
    
        out_stream.write(f"\n beta (1/KT) {beta:.8f}\n")
        
        for ibox in range(job.num_boxes):
            fld.setup()

            energy_new = fld.calculate_energy(basin[ibox])
           
            totalEnergy.append(energy_new)
            totalEnergy[ibox].print_energy(ibox, out_stream)

        numSteps = 1

        while numSteps <= job.mcSteps:

            ibox = int(job.num_boxes * np.random.random())

            choice = int(self.numMCMoves * np.random.random())

            selection = self.mcMoveList[choice]

            if selection == 1:
                self.move_atom(basin[ibox], fld, totalEnergy[ibox], beta, out_stream)

            elif selection == 2:
                self.swapAtoms(basin[ibox], fld, totalEnergy[ibox], job, beta, out_stream)

            elif selection == 3:
                self.transmutateAtoms(basin[ibox], fld, totalEnergy[ibox], spec, beta, out_stream)

            elif selection == 5:               
                 self.move_volume(basin[ibox], fld, totalEnergy[ibox], spec, job, beta, out_stream)

            for ib in range(job.num_boxes):
                stats[ib].sample(job.equilSteps, numSteps, totalEnergy[ib], basin[ib].get_volume(), basin[ib].get_cell().flatten(), out_stream)
                type_stats[ib].sample_types(numSteps, job.equilSteps, basin[ib], spec)

                #if self.numSemiWidom > 0:
                #    chem_stats.sample_semi_widom(chemPotAtom1, chemPotAtom2, numSemiWidom, numSteps, job.equilSteps)

            

                if numSteps % job.printFreq == 0:
                    stats[ib].check_point(numSteps, job.equilSteps, 0.0, out_stream)
                    type_stats[ib].check_point_types(spec, out_stream)
                    #if self.numSemiWidom > 0:
                    #    chem_stats[ib].check_point_semi_widom(self.semiWidomType1, self.semiWidomType2, self.numSemiWidom, spec, job.temperature, out_stream)

                if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                    basin[ib].writeBasisXYZ(spec, xyzStream)
                    xyzStream.flush()

                    with open(f"restart", "w") as restartStream:
                        basin[ib].dump_basis(spec, 0.0, 0.0, job.mc_steps, restartStream)

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
                restart_stream = open("restart", "w") 
            
                
                for ib in range(job.num_boxes):
                    checkEnergy = fld.calculate_energy(basin[ibox])

                    eDiff = checkEnergy.get_total_energy() - totalEnergy[ib].get_total_energy()

                    if abs(eDiff) > 1.0e-6:
                        out_stream.write(f"\n sanity check failed on iteration {numSteps} !!!!!!!\n")
                        out_stream.write(f" total diff {eDiff.totalEnergy:.10e}\n")
                        
                    totalEnergy[ib] = checkEnergy

                    #basin[ib].dump_basis(spec, 0.0, 0.0, numSteps, restart_stream)

                restart_stream.flush()
                restart_stream.close()

            numSteps += 1

        out_stream.write("\n\n *****************************************************************************************************\n")
        out_stream.write(f" final energy of system containing {basin[cycle].get_number_of_atoms()} atoms\n")
        out_stream.write(" *****************************************************************************************************\n")

        
        for ibox in range(job.num_boxes):
            final_energy = Energy()
            final_energy = fld.calculate_energy(basin[ibox])

            final_energy.print_energy(ibox, out_stream)

            out_stream.write("\n final sanity check")
            checkEnergy = final_energy - totalEnergy[ibox]
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

        #if self.attemptedGCFragRemoves > 0:
        #    out_stream.write(f"\n the number of fragment removes {self.successFulGCFragRemoves} attempted {self.attemptedGCFragRemoves}\n")
        #if self.attemptedGCFragInsert > 0:
        #    out_stream.write(f"\n the number of fragment inserts {self.successFulGCFragInsert} attempted {self.attemptedGCFragInsert}\n")

        out_stream.flush()

    

    def move_atom(self, basin: Atoms, fld: Field, total_energy: Energy, beta: float, out_stream):
        
        old_pos = np.zeros(3, dtype=np.float64)
        
        
        choice = int(np.random.random() * self.noAtomMovers)
        typ = self.atomMoveTypes[choice]
        
        energy_old = Energy()
        energy_new = Energy()
        
        atm = self.select_atom(basin)
        
        if atm < 0:
            return
        
        self.totalAtomMoves += 1
        self.attemptedAtomMoves[typ] += 1

        energy_old = fld.calculate_energy(basin)
        
        old_pos[:] = basin.positions[atm][:]
        delta_pos = self.atom_displacement(self.distance_atom_max[typ])
        self.make_atom_move(basin, atm, delta_pos)
        
        energy_new = fld.calculate_energy(basin)
        
        deltaV = energy_new.get_total_energy() - energy_old.get_total_energy()
        deltaVB = deltaV * beta
        
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
            total_energy.totalEnergy = energy_new.totalEnergy
            
            self.noAtomMoves[typ] += 1
            self.successfulAtomMoves += 1
            
        else:
            self.reject_atom_move(basin, atm, old_pos)


    def move_volume(self, basin: Atoms, fld: Field, totalEnergy: Energy, spec: Species, job: JobControl, beta: float, out_stream):
        indx = 0
        vol_new = 1.0
        natoms = len(basin)

        betaInv = 1.0 / beta

        oldEnergy = totalEnergy
        vol_old = self.cell_size(basin)

        bulks = np.ones(3, dtype=np.float64)
        
        if job.volMoveSymmetry == 0:
            indx = 0
            maxVol = self.maxVolChange[0]
            vol_new = self.expand_cell_cubic(basin, bulks, maxVol)
        elif job.volMoveSymmetry == 1:
            indx = int(2.0 * np.random.random())
            vol_new = self.expand_cell_tetragonal(basin, indx, bulks, self.maxVolChange)
        elif job.volMoveSymmetry == 2:
            indx = int(3.0 * np.random.random())
            vol_new = self.expand_cell_orthorhombic(basin, indx, bulks, self.maxVolChange)
        else:
            pass  # Implement other volume changes if needed
        
        self.attemptedVolChange[indx] += 1
        self.totalVolChanges[indx] += 1

        new_energy = Energy()
        new_energy = fld.calculate_energy(basin)

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
            self.restore_cell(basin, bulks, indx)
            
    def getRandomNumber(self):
        return np.random.random()

    def swapAtoms(self, basin: Atoms, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        self.attemptedSwaps += 1

        j = int(self.getRandomNumber() * self.numSwaps)
        
        atm1 = self.select_atom_of_type(basin, self.swapType1[j])
        atm2 = self.select_atom_of_type(basin, self.swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy

        #print("swapping", atm1, atm2, basin.chem_symbols[atm1], basin.chem_symbols[atm2])
        self.swap_atom_positions(basin, atm1, atm2)

        new_energy = Energy()
        #calculating energy does not take into account pair potential methods
        #new_energy = fld.calculate_swap_energy(basin.pos_r, basin.lat_vector, basin.rcp_vector, basin.charge, 
        #                    basin.atm_label, basin.chem_symbols, basin.frozen, basin.number_of_atoms, old_energy)
        new_energy = fld.calculate_energy(basin)

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

        else:
            self.swap_atom_positions(basin, atm1, atm2)

    
    def transmutateAtoms(self, basin: Atoms, fld: Field, totalEnergy: Energy, spec: Species, beta: np.float64, out_stream):

        old_energy = totalEnergy

        choice = self.getRandomNumber()

        if choice < 0.5:
            self.attemptForwardMutations += 1
            k = int(self.getRandomNumber() * self.numTrans)
            deltaMu = self.transmuteChemPot[k]

            numTypes1 = self.find_num_types(basin, self.transType1[k])
            if numTypes1 == 0:
                return

            numTypes2 = self.find_num_types(basin, self.transType2[k])
            weight = float(numTypes1) / float(numTypes2 + 1)

            atm = self.select_atom_of_type(basin, self.transType1[k])
            if atm == -1:
                return

            self.mutate_atom(basin, atm, self.transType2[k])

            new_energy = Energy()
            new_energy = fld.calculate_energy(basin)

            deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
            deltaVB = beta * (deltaV - deltaMu)
            prob = weight * math.exp(-deltaVB)
            print("forward",atm,self.transType1[k],self.transType1[k],deltaV,(deltaV - deltaMu),deltaVB,prob )
            if self.getRandomNumber() < prob:
                totalEnergy.totalEnergy = new_energy.totalEnergy
                self.forwardMutations += 1
            else:
                self.mutate_atom(basin, atm, self.transType1[k])
        else:
            self.attemptBackwardMutations += 1
            k = int(self.getRandomNumber() * self.numTrans)
            deltaMu = self.transmuteChemPot[k]

            numTypes2 = self.find_num_types(basin, self.transType2[k])
            if numTypes2 == 0:
                return

            numTypes1 = self.find_num_types(basin, self.transType1[k])
            weight = float(numTypes2) / float(numTypes1 + 1)

            atm = self.select_atom(basin, self.transType2[k])
            if atm == -1:
                return

            #oldEnergy = Energy()
            #fld.calculateAtomEnergy(atm, basin.pos_x, basin.pos_y, basin.pos_z, basin.lat_vector, basin.rcp_vector, basin.charge, 
            #                        basin.atm_label, basin.frozen, basin.number_of_atoms, oldEnergy)

            self.mutate_atom(basin, atm, self.transType1[k])

            new_energy = Energy()
            new_energy = fld.calculate_energy(basin)

            deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
            deltaVB = beta * (deltaV + deltaMu)
            prob = weight * math.exp(-deltaVB)
            print("back",atm,self.transType2[k],self.transType1[k],deltaV,(deltaV + deltaMu),deltaVB,prob )
            if self.getRandomNumber() < prob:
                totalEnergy.totalEnergy = new_energy.totalEnergy
                self.backwardMutations += 1
            else:
                self.mutate_atom(basin, atm, self.transType2[k])

    def select_atom(self, basin: Atoms) -> int:
        
        choice = -1
            
        choice = int(len(basin) * np.random.random())

        
        return choice     
        
    def select_atom_of_type(self, basin: Atoms, typ: str) -> int:
        
        atm = -1 
        choice = -1
            
        found = False

        atm_list = []
    
        for i in range(len(basin)):
            if basin.symbols[i] == typ:
                found = True
                atm_list.append(i)
        
        
        if found == False:
            return atm 

        choice = int(len(atm_list) * np.random.random())

        atm = atm_list[choice]

        return atm

    def find_num_types(self, bas:Atoms, typ) -> int:
        num_typ = 0
        for i in range(len(bas)):
            if typ == bas.symbols[i]:
                num_typ += 1

        return num_typ
    def atom_displacement(self, dist_max: float) -> np.ndarray:
        r = [self.getRandomNumber() for _ in range(3)]

        delta_pos = np.zeros(3, dtype=np.float64)
        delta_pos[0] = (r[0] - 0.5) * dist_max
        delta_pos[1] = (r[1] - 0.5) * dist_max
        delta_pos[2] = (r[2] - 0.5) * dist_max

        #print ("old", self.pos_r[atom,:], dist_max)
        #print("new ", new_pos[:])

        return delta_pos

    def make_atom_move(self, basin: Atoms, atm: int, delta_pos: np.ndarray):
        basin.positions[atm][:] = basin.positions[atm][:] + delta_pos[:]

    def reject_atom_move(self, basin: Atoms, atm: int, old_pos: np.ndarray):
        basin.positions[atm][:] = old_pos[:]

    def swap_atom_positions(self, basin: Atoms, atm1: int, atm2: int):
        tmp = np.zeros(3, dtype=np.float64)
        tmp[:] =basin.positions[atm1][:]
        basin.positions[atm1][:] = basin.positions[atm2][:]
        basin.positions[atm2][:] = tmp[:]

    def mutate_atom(self, atm: int, typ: int, spec: Species):
        ele = spec.get_species(typ)
        self.chem_symbols[atm] = ele.name
        self.mass[atm] = ele.mass
        self.charge[atm] = ele.charge
        self.atm_label[atm] = typ

    def expand_cell_cubic(self, basin, bulks, max_vol_change):
        r = np.random.random()
        
        scale = 1.0 + (r - 0.5) * max_vol_change

        cell = basin.get_cell()
        basin.cell[0][0] *= scale
        basin.cell[1][1] *= scale
        basin.cell[2][2] *= scale

        volume = self.cell_size(basin)
        #print("new volume ", scale, volume, max_vol_change)
 
        bulks[0] = scale
        bulks[1] = scale
        bulks[2] = scale

        self.scale_positions(basin, bulks)

        return volume
   
    def expand_cell_tetragonal(self, basin, indx, bulks, max_vol_change):
        r = np.random.random()
        #bulks = np.ones(3, dtype=np.float64)
        cell = basin.get_cell()
        
        scale = 1.0 + (r - 0.5) * max_vol_change

        if indx == 0:
            cell[0][0] *= scale
            cell[1][1] *= scale
            bulks[0] = scale
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        volume = self.cell_size()
        
        self.scale_positions(basin, bulks)

        return volume
       
    def expand_cell_orthorhombic(self, basin, indx, bulks, max_vol_change):
        r = np.random.random()
        cell = basin.get_cell()

        scale = 1.0 + (r - 0.5) * max_vol_change

        if indx == 0:
            cell[0][0] *= scale
            bulks[0] = scale
            bulks[1] = 1.0
            bulks[2] = 1.0
        elif indx == 1:
            cell[1][1] *= scale
            bulks[0] = 1.0
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        volume = self.cell_size()
        
        self.scale_positions(basin, bulks)

        return volume

    def scale_positions(self, basin:Atoms, bulks):

        for i in range(len(basin)):
           basin.positions[i,:] *= bulks[:]

    def cell_size(self, basin) -> np.float64:
        
        return basin.get_volume()
    
    def restore_cell(self, basin, bulks, indx):
        cell = basin.get_cell()

        if indx == 0:
            scale = 1.0 / bulks[0]
            cell[0][0] *= scale
            bulks[0] = scale
            bulks[1] = 1.0
            bulks[2] = 1.0
        elif indx == 1:
            scale = 1.0 / bulks[1]
            cell[1][1] *= scale
            bulks[0] = 1.0
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            scale = 1.0 / bulks[2]
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        self.scale_positions(basin, bulks)