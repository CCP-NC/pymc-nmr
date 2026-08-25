
import numpy as np

from ase.units import fs

class JobControl:
    def __init__(self):
        self.num_cycles = 1
        self.mc_steps = 0
        self.num_boxes = 1
        self.wrap = False
        self.extPressure = 0.0
        self.temperature = 0.0
        self.maxDistance: np.float64 = 0.001
        self.accAtomMoveUpdate = 1000
        self.acceptAtomMoveRatio = 0.37
        self.acceptVolUpdate = 1000
        self.acceptVolMoveRatio = 0.37
        self.maxVolDisplacement: np.float64 = 0.1
        self.frozen_types = []
        self.sanityCheckFreq = 1000
        self.printFreq = 1
        self.dumpArchive = False
        self.archiveFrequency = 1000
        self.save_downhill = False
        self.equilSteps = 0
        self.atomMoveFreq = 0
        self.moveTypes = []
        self.volMoveFreq = 0
        self.num_swap_atoms = 0
        self.swapFrequency = 0
        self.swapType1 = []
        self.swapType2 = []

        self.num_combi_swap_atoms = 0
        self.combi_swapFrequency = 0
        self.combi_swapType1 = []
        self.combi_swapType2 = []
        self.combi_swapType3 = []
        self.combi_swap_dist = 3.0

        self.num_transmutate_atoms = 0
        self.transmutateFrequency = 0
        self.mutateType1 = []
        self.mutateType2 = []
        self.transmuteChemPot = []
        self.volMoveSymmetry = 0
        self.structure_method = "basinhop"
        self.restart = False
        self.max_force = 1.0e-3

        self.writestats = False
        self.writestats_freq = 10

        #md parameters
        self.mdMoveFreq = 0
        self.mdtimestep = 2.0 * fs
        self.mdtemperature_K = 1000
        self.mdfriction = 0.01 / fs
        self.mdsteps = 1000

        #relaxation parameters
        self.relmethod = "lbfgs"
        self.relsteps = 1000
        self.reltol = 1e-3
        self.relstyle = "conp"

        #grid parameters
        self.gridx = 10
        self.gridy = 10
        self.gridz = 10
        self.grid_cut = 2.0

        # MPI/parallel parallelisation suffix
        self.rank_suffix = ""
        
    def write_mc_control(self, out_io):
        """
        writes out to a file the control parameters that are pertinent to the MC calculation

        Parameters
        ----------

        out_io : io
            the file stream for the writing of data
        """

        out_io.write(f"\n the number of steps                      {self.mcSteps} \n")
        out_io.write(f"\n pressure                                 {self.extPressure} \n")
        out_io.write(f"\n temperature                              {self.temperature} \n")
        out_io.write(f"\n max atom displacement                    {self.maxDistance} \n")
        out_io.write(f"\n max volume displacement                  {self.maxVolDisplacement} \n")
        out_io.write(f"\n sanity check frequency                   {self.sanityCheckFreq} \n")
        out_io.write(f"\n print frequency                          {self.printFreq} \n")
        if self.dumpArchive == True:
            out_io.write(f"\n archive configuration frequency          {self.archiveFrequency} \n")
        out_io.write(f"\n equilibration steps                      {self.equilSteps} \n")
        if len(self.moveTypes) > 0:
            out_io.write(f"\n atom types to be moved \n")
            for i in range(len(self.moveTypes)):
                out_io.write(f"\n         {self.moveTypes[i]}")

            out_io.write(f"\n atom move frequency                      {self.atomMoveFreq} \n")

            out_io.write(f"\n maximum atom displacement {self.maxDistance} \n")
            out_io.write(f"\n max displacement update every {self.accAtomMoveUpdate} steps \n")
            out_io.write(f"\n acceptance move ratio {self.acceptAtomMoveRatio} \n")
        
        if self.volMoveFreq > 0:
            out_io.write(f"\n volume move frequency                    {self.volMoveFreq} \n")

            if self.volMoveSymmetry == 0:
                out_io.write(f"\n volume move symmetry                     cubic \n")
            elif self.volMoveSymmetry == 1:
                out_io.write(f"\n volume move symmetry                     tetragonal \n")
            elif self.volMoveSymmetry == 2:
                out_io.write(f"\n volume move symmetry                     orthorhombic \n")
            elif self.volMoveSymmetry == 3:
                out_io.write(f"\n volume move symmetry                     vectors \n")
            else:
                out_io.write(f"\n volume move symmetry not recognised \n")
                out_io.flush()
                exit()

            
            out_io.write(f"\n maximum volume expansion {self.acceptVolUpdate} \n")
            out_io.write(f"\n volume acceptance ratio {self.acceptVolMoveRatio} \n")
            out_io.write(f"\n max volume displacement updated every {self.maxVolDisplacement} steps \n")

        if self.num_swap_atoms > 0:
            out_io.write(f"\n atom types to be swapped \n")
            for i in range(len(self.swapType1)):
                out_io.write(f"\n         {self.swapType1[i]}  {self.swapType2[i]}")
                    
            out_io.write(f"\n atom swap frequency                      {self.swapFrequency} \n")
 
        if self.writestats:
            out_io.write(f"\n energy data will be written to a file every {self.writestats_freq} steps \n")
            
    def write_bh_control(self, out_io):
        """
        writes out to a file the control parameters that are pertinent to the basin hopping calculation

        Parameters
        ----------

        out_io : io
            the file stream for the writing of data
        """

        out_io.write(f"\n the number of steps                      {self.mcSteps} \n")
        out_io.write(f"\n temperature                              {self.temperature} \n")
        out_io.write(f"\n print frequency                          {self.printFreq} \n")
        if self.dumpArchive == True:
            out_io.write(f"\n archive configuration frequency          {self.archiveFrequency} \n")
        out_io.write(f"\n equilibration steps                      {self.equilSteps} \n")
        if self.save_downhill == True:
            out_io.write(f"\n energy moves to a lower energy will be saved \n")
            
        if self.num_swap_atoms > 0:
            out_io.write(f"\n atom types to be swapped \n")
            for i in range(len(self.swapType1)):
                out_io.write(f"\n         {self.swapType1[i]}  {self.swapType2[i]}")
                    
            out_io.write(f"\n atom swap frequency                      {self.swapFrequency} \n")

        if self.num_combi_swap_atoms > 0:
            out_io.write(f"\n atom types to be swapped \n")
            for i in range(len(self.combi_swapType1)):
                out_io.write(f"\n         {self.combi_swapType1[i]}  {self.combi_swapType2[i]} combined with type {self.combi_swapType3[i]}")
                    
            out_io.write(f"\n combined atom swap frequency             {self.combi_swapFrequency} \n")
            out_io.write(f"\n combined atom swap distance              {self.combi_swap_dist} \n")

            out_io.write(f"\n grid parameters {self.gridx} {self.gridy} {self.gridz} \n")
            out_io.write(f"\n grid exclusion size {self.grid_cut} \n")

        if self.mdMoveFreq > 0:
            out_io.write(f"\n an MD will be undertaken with frequency  {self.mdMoveFreq} \n")
            out_io.write(f"\n  MD timestep                             {self.mdtimestep} \n")
            out_io.write(f"\n  MD temperature                          {self.mdtemperature_K} \n")
            out_io.write(f"\n  MD friction                             {self.mdfriction} \n")
            out_io.write(f"\n  MD steps                                {self.mdsteps} \n")

        
        out_io.write(f"\n  relaxation method                       {self.relmethod} \n")
        out_io.write(f"\n  relaxation steps                        {self.relsteps} \n")
        out_io.write(f"\n  relaxation tolerance                    {self.reltol} \n")

        if self.relstyle == "conp":
            out_io.write(f"\n cell lengths and angles will be allowed to change during relaxation \n")
        elif self.relstyle == "cona":
            out_io.write(f"\n cell angles will be fixed during relaxation \n")
        elif self.relstyle == "conv":
            out_io.write(f"\n cell lengths and angles will be fixed during relaxation \n")
        else:
            out_io.write(f"\n wrong relaxation style entered. Valid options are : conp, cona, conv \n")

        if self.writestats:
            out_io.write(f"\n energy data will be written to a file every {self.writestats_freq} steps \n")

    def write_airss_control(self, out_io):
        """
        writes out to a file the control parameters that are pertinent to the AIRSS-style calculation

        Parameters
        ----------

        out_io : io
            the file stream for the writing of data
        """
        out_io.write(f"\n the number of steps                      {self.mcSteps} \n")
        
        out_io.write(f"\n print frequency                          {self.printFreq} \n")
        
            
        if self.num_swap_atoms > 0:
            out_io.write(f"\n atom types to be swapped \n")
            for i in range(len(self.swapType1)):
                out_io.write(f"\n         {self.swapType1[i]}  {self.swapType2[i]}")
                    
            out_io.write(f"\n atom swap frequency                      {self.swapFrequency} \n")

        out_io.write(f"\n  relaxation method                       {self.relmethod} \n")
        out_io.write(f"\n  relaxation steps                        {self.relsteps} \n")
        out_io.write(f"\n  relaxation tolerance                    {self.reltol} \n")

        if self.relstyle == "conp":
            out_io.write(f"\n cell lengths and angles will be allowed to change during relaxation \n")
        elif self.relstyle == "cona":
            out_io.write(f"\n cell angles will be fixed during relaxation \n")
        elif self.relstyle == "conv":
            out_io.write(f"\n cell lengths and angles will be fixed during relaxation \n")   
        else:
            out_io.write(f"\n wrong relaxation style entered. Valid options are : conp, cona, conv \n")
        
        if self.writestats:
            out_io.write(f"\n energy data will be written to a file every {self.writestats_freq} steps \n")

    def load_from_schema(self, cfg):
        self.num_cycles = cfg.job.num_cycles
        self.mc_steps = cfg.monte.steps
        self.mcSteps = cfg.monte.steps
        self.wrap = cfg.job.wrap
        self.extPressure = cfg.job.pressure
        self.temperature = cfg.job.temperature
        self.maxDistance = cfg.monte.max_distance
        self.accAtomMoveUpdate = cfg.monte.distance_update_freq
        self.acceptAtomMoveRatio = cfg.monte.distance_ratio
        self.acceptVolUpdate = cfg.monte.vol_update_freq
        self.acceptVolMoveRatio = cfg.monte.vol_ratio
        self.maxVolDisplacement = cfg.monte.max_vol_displacement
        self.frozen_types = cfg.job.frozen_types
        self.sanityCheckFreq = cfg.job.sanity_check_freq
        self.printFreq = cfg.job.print_freq
        self.dumpArchive = cfg.job.dump_archive
        self.archiveFrequency = cfg.job.archive_frequency
        self.save_downhill = cfg.job.save_downhill
        self.equilSteps = cfg.job.equil_steps
        
        # atom move config
        self.atomMoveFreq = cfg.monte.atom_move_freq
        self.moveTypes = cfg.monte.move_types
        
        # vol move config
        self.volMoveFreq = cfg.monte.vol_move_freq
        # map volMoveSymmetry string to int if needed
        sym_map = {"cubic": 0, "tetragonal": 1, "orthorhombic": 2, "vectors": 3}
        self.volMoveSymmetry = sym_map.get(cfg.monte.vol_move_symmetry.lower(), 0)
        
        # swap config
        self.swapType1 = []
        self.swapType2 = []
        self.swapFrequency = 0
        
        self.combi_swapType1 = []
        self.combi_swapType2 = []
        self.combi_swapType3 = []
        self.combi_swapFrequency = 0
        self.combi_swap_dist = 3.0
        
        for move in cfg.moves.swaps:
            if move.counterion is not None:
                self.combi_swapFrequency += move.frequency
                self.combi_swap_dist = move.max_distance
                for t1, t2 in move.pairs:
                    self.combi_swapType1.append(t1)
                    self.combi_swapType2.append(t2)
                    self.combi_swapType3.append(move.counterion)
            else:
                self.swapFrequency += move.frequency
                for t1, t2 in move.pairs:
                    self.swapType1.append(t1)
                    self.swapType2.append(t2)
                    
        self.num_swap_atoms = len(self.swapType1)
        self.num_combi_swap_atoms = len(self.combi_swapType1)

        # transmutate
        self.transmutateFrequency = cfg.monte.transmutate_frequency
        self.mutateType1 = cfg.monte.mutate_types_1
        self.mutateType2 = cfg.monte.mutate_types_2
        self.num_transmutate_atoms = len(self.mutateType1)
        self.transmuteChemPot = cfg.monte.transmute_chem_pot

        self.structure_method = cfg.job.structure_method
        self.restart = cfg.job.restart
        self.max_force = cfg.job.max_force
        self.writestats = cfg.job.writestats
        self.writestats_freq = cfg.job.writestats_freq

        # md parameters
        self.mdMoveFreq = cfg.md.move_freq
        self.mdtimestep = cfg.md.timestep_fs * fs
        self.mdtemperature_K = cfg.md.temperature_k
        self.mdfriction = cfg.md.friction_fs / fs
        self.mdsteps = cfg.md.steps

        # relaxation parameters
        self.relmethod = cfg.relax.method
        self.relsteps = cfg.relax.steps
        self.reltol = cfg.relax.tol
        self.relstyle = cfg.relax.style

        # grid parameters
        self.gridx = cfg.grid.x
        self.gridy = cfg.grid.y
        self.gridz = cfg.grid.z
        self.grid_cut = cfg.grid.cut


