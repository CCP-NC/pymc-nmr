import numpy as np
import math
from typing import List

from ase import Atoms
from ase.io import write

from Energy import Energy
from Field import Field
from Species import Species
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
        self.successfulDownSwaps = 0
        self.attemptedSwaps = 0
        self.successfulUpSwaps = 0
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
    
    def run(self, spec: Species, fld: Field, job: JobControl, stats: Statistics, type_stats: TypeStatistics, basin, numSteps, cycle, initialise, out_stream):

        totalEnergy = []
        checkEnergy = Energy()

        # Initiate the statistics
        for i in range(job.num_boxes):
            stats[i].zero(1000, 0.0, False)
            type_stats[i].zero_types(spec.get_num_species(), False) 

        beta = 1.0 / (job.temperature * BOLTZMANN)
    
        out_stream.write(f"\n beta (1/KT) {beta:.8f}\n")
        
        for ibox in range(job.num_boxes):
            fld.setup()

            energy_new = fld.calculate_energy_relax(basin[ibox], job.relmethod, job.relsteps, job.reltol)
           
            totalEnergy.append(energy_new)
            totalEnergy[ibox].print_energy(ibox, out_stream)

            if job.restart == False:
                write(filename="archive.xyz", images=basin, format="extxyz")

        numSteps = 1

        while numSteps <= job.mcSteps:

            ibox = int(job.num_boxes * np.random.random())

            choice = int(self.numMCMoves * np.random.random())

            selection = self.mcMoveList[choice]

            if selection == 1:
                basin[ibox] = self.run_md(basin[ibox], fld, totalEnergy[ibox], job, beta, out_stream)

            elif selection == 2:
                basin[ibox] = self.swapAtoms_relax(basin[ibox], fld, totalEnergy[ibox], job, beta, out_stream)

            energy_new = fld.calculate_energy(basin[ibox])
            print("energy in main routine ", energy_new.totalEnergy)

            for ib in range(job.num_boxes):
                stats[ib].sample(job.equilSteps, numSteps, totalEnergy[ib], basin[ib].get_volume(), basin[ib].get_cell().flatten(), out_stream)
                type_stats[ib].sample_types(numSteps, job.equilSteps, basin[ib], spec)

                if numSteps % job.printFreq == 0:
                    stats[ib].check_point(numSteps, job.equilSteps, 0.0, out_stream)
                    type_stats[ib].check_point_types(spec, out_stream)
                    
            if job.dumpArchive and numSteps % job.archiveFrequency == 0:
                write(filename="archive.xyz", images=basin, format="extxyz", append=True)
            

                   
            #if numSteps > job.equilSteps and job.sampleBasin and numSteps % job.sampleBasinFreq == 0:
            #    basin[cycle].samplePositions()

            if numSteps % job.sanityCheckFreq == 0: 
                write(filename="restart.xyz", images=basin, format="extxyz")
              
                for ib in range(job.num_boxes):
                    checkEnergy = fld.calculate_energy(basin[ibox])

                    eDiff = checkEnergy.get_total_energy() - totalEnergy[ib].get_total_energy()

                    if abs(eDiff) > 1.0e-6:
                        out_stream.write(f"\n sanity check failed on iteration {numSteps} !!!!!!!\n")
                        out_stream.write(f" total diff {eDiff.totalEnergy:.10e}\n")
                        
                    totalEnergy[ib] = checkEnergy

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
       
        if self.numSwaps > 0:
            swapRatio = (self.successfulUpSwaps + self.successfulDownSwaps) / self.attemptedSwaps
            successfulSwaps = self.successfulUpSwaps + self.successfulDownSwaps
            out_stream.write(f"\n swaps : attempted, successful and ratio {self.attemptedSwaps} {successfulSwaps} {swapRatio:.10e}\n")
            out_stream.write(f"\n swaps : Downhill and Uphill {self.successfulDownSwaps} {self.successfulUpSwaps} \n")

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

        write(filename="restart.xyz", images=basin, format="extxyz")
        
        out_stream.flush()


    def swapAtoms_relax(self, basin: Atoms, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        self.attemptedSwaps += 1

        j = int(np.random.random() * self.numSwaps)
        
        atm1 = self.select_atom_of_type(basin, self.swapType1[j])
        atm2 = self.select_atom_of_type(basin, self.swapType2[j])

        if atm1 == -1 or atm2 == -1:
            return

        old_energy = Energy()
        old_energy.totalEnergy = totalEnergy.totalEnergy
        old_pos = np.zeros((len(basin), 3), dtype=np.float64)
        symbols = []
        cell = np.zeros((3,3), dtype=np.float64)
        
        #np.copyto(old_pos, basin.positions)
        for i in range(len(basin)):
            symbols.append(basin.symbols[i])
            for j in range(3):
                old_pos[i,j] = basin.positions[i,j]
        for i in range(3):
            for j in range(3):
                cell[i,j] = basin.cell[i,j]
        new_basin = Atoms(symbols = symbols, positions=old_pos, cell= cell, pbc=True)
        new_basin.wrap()

        #print("swapping", atm1, atm2, basin.chem_symbols[atm1], basin.chem_symbols[atm2])
        print("swapping ", atm1, basin.symbols[atm1], atm2, basin.symbols[atm2])
        print("pos atm1 ", basin.symbols[atm1], basin.positions[atm1,:])
        print("pos atm2 ", basin.symbols[atm2], basin.positions[atm2,:])
        self.swap_atom_types(new_basin, atm1, atm2)

        new_energy = Energy()
       
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol)
        
        deltaV = new_energy.get_total_energy() - old_energy.get_total_energy()
        deltaVB = beta * deltaV
        print("swap ", old_energy.get_total_energy(), new_energy.get_total_energy(), deltaV, deltaVB)
        accept = False
        arg = np.random.random()
        if arg < np.exp(-deltaVB):
            accept = True
        
        if accept:
            totalEnergy.totalEnergy = new_energy.totalEnergy
            if deltaV < 0.0:
                self.successfulDownSwaps += 1
                if job.save_downhill:
                    write(filename="downhill.xyz", images=basin, format="extxyz", append=True)
            else:
                self.successfulUpSwaps += 1
            print("swap accepted")
            return new_basin
        else:
            #np.copyto(basin.positions, old_pos)
            #self.swap_atom_types(basin, atm1, atm2)
            #for i in range(len(basin)):
            #    for j in range(3):
            #        basin.positions[i,j] = old_pos[i,j]
            print("swap failed", new_basin.get_potential_energy())
            print("pos atm1 ", basin.symbols[atm1], basin.positions[atm1,:])
            print("pos atm2 ", basin.symbols[atm2], basin.positions[atm2,:])
            return basin
        
    def run_md(self, basin: Atoms, fld: Field, totalEnergy: Energy, job: JobControl, beta: np.float64, out_stream):
        
        self.md_runs += 1

        j = int(np.random.random() * self.numSwaps)
        
        
        old_pos = np.zeros((len(basin), 3), dtype=np.float64)
        symbols = []
        cell = np.zeros((3,3), dtype=np.float64)
        
        #np.copyto(old_pos, basin.positions)
        for i in range(len(basin)):
            symbols.append(basin.symbols[i])
            for j in range(3):
                old_pos[i,j] = basin.positions[i,j]
        for i in range(3):
            for j in range(3):
                cell[i,j] = basin.cell[i,j]
        new_basin = Atoms(symbols = symbols, positions=old_pos, cell= cell, pbc=True)
        new_basin.wrap()

        #run an md simulation
        fld.run_md(new_basin, job.timestep, job.mdtemperature_K, job.mdfriction, job.mdsteps)

        new_energy = Energy()
        #the energy needs to be relaxed to get the "new" energy
        new_energy = fld.calculate_energy_relax(new_basin, job.relmethod, job.relsteps, job.reltol)
        
        totalEnergy.totalEnergy = new_energy.totalEnergy
        
        return new_basin

       
    def transmutateAtoms(self, basin: Atoms, fld: Field, totalEnergy: Energy, spec: Species, beta: np.float64, out_stream):

        old_energy = totalEnergy

        choice = np.random.random()

        if choice < 0.5:
            self.attemptForwardMutations += 1
            k = int(np.random.random() * self.numTrans)
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
            if np.random.random() < prob:
                totalEnergy.totalEnergy = new_energy.totalEnergy
                self.forwardMutations += 1
            else:
                self.mutate_atom(basin, atm, self.transType1[k])
        else:
            self.attemptBackwardMutations += 1
            k = int(np.random.random() * self.numTrans)
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
            if np.random.random() < prob:
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
    

    def swap_atom_positions(self, basin: Atoms, atm1: int, atm2: int):
        #tmp = np.zeros(3, dtype=np.float64)
        #tmp[:] =basin.positions[atm1,:]
        #    basin.positions[atm1][:] = basin.positions[atm2][:]
        #    basin.positions[atm2][:] = tmp[:]
        for j in range(3):
            tmp = basin.positions[atm1,j]
            basin.positions[atm1,j] = basin.positions[atm2,j]
            basin.positions[atm2,j] = tmp

    
    def swap_atom_types(self, basin: Atoms, atm1: int, atm2: int):
        tmp = basin.symbols[atm1]
        basin.symbols[atm1] = basin.symbols[atm2]
        basin.symbols[atm2] = tmp

    def mutate_atom(self, atm: int, typ: int, spec: Species):
        ele = spec.get_species(typ)
        self.chem_symbols[atm] = ele.name
        self.mass[atm] = ele.mass
        self.charge[atm] = ele.charge
        self.atm_label[atm] = typ

    
    def cell_size(self, basin) -> np.float64:
        
        return basin.get_volume()
    
