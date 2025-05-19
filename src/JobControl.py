
import numpy as np

from ase.units import fs

class JobControl:
    def __init__(self):
        self.num_cycles = 1
        self.mc_steps = 0
        self.num_boxes = 1
        self.extPressure = 0.0
        self.temperature = 0.0
        self.maxDistance: np.float64 = 0.0
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
        self.num_transmutate_atoms = 0
        self.transmutateFrequency = 0
        self.mutateType1 = []
        self.mutateType2 = []
        self.transmuteChemPot = []
        self.volMoveSymmetry = 0
        self.structure_method = "basinhop"
        self.restart = False
        self.max_force = 1.0e-3

        #md parameters
        self.mdMoveFreq = 0
        self.timestep = 2.0 * fs
        self.mdtemperature_K = 1000
        self.mdfriction = 0.01 / fs
        self.mdsteps = 1000

        #relaxation parameters
        self.relmethod = "lbfgs"
        self.relsteps = 1000
        self.reltol = 1e-3
        
    def check_job_control(self, out_stream):
        pass  # Implement parameter check logic if needed

    def read_job_control(self, in_stream, out_stream):
        while True:
            line = in_stream.readline()
            if not line:
                break
            out_stream.write(f" input line: {line.strip()}\n")
            if line[0] == '#':
                continue

            words = line.split()

            if not words:
                continue

            keyWord = words[0].lower()

            if keyWord == "steps":
                self.mcSteps = int(words[1])
            elif keyWord == "cycles":
                self.num_cycles = int(words[1])
            elif keyWord == "pressure":
                self.extPressure = float(words[1])
            elif keyWord == "temperature":
                self.temperature = float(words[1])
            elif keyWord == "maxdistance":
                self.maxDistance = np.float64(words[1])
            elif keyWord == "maxvolume":
                self.maxVolDisplacement = np.float64(words[1])
            elif keyWord == "freeze":
                num = int(words[1])
                for _ in range(num):
                    line = in_stream.readline()
                    words = self.split(line)
                    self.frozen_types.append(words[0])
            elif keyWord == "check":
                self.sanityCheckFreq = int(words[1])
            elif keyWord == "print":
                self.printFreq = int(words[1])
            elif keyWord == "archive":
                self.dumpArchive = True
            elif keyWord == "archivefrequency":
                self.archiveFrequency = int(words[1])
            elif keyWord == "equilsteps":
                self.equilSteps = int(words[1])
            elif keyWord == "method":
                self.structure_method = words[1]
            elif keyWord == "restart":
                self.restart = True
            elif keyWord == "savedownhill":
                self.save_downhill = True
            elif keyWord == "move":
                subWord = words[1]
                if subWord == "atoms":
                    num = int(words[2])
                    self.atomMoveFreq = int(words[3])
                    for _ in range(num):
                        line = in_stream.readline()
                        words = line.split()
                        self.moveTypes.append(words[0])
                elif subWord == "volume":
                    self.volMoveFreq = int(words[2])
                elif subWord == "swap":
                    self.num_swap_atoms = int(words[2])
                    self.swapFrequency = int(words[3])
                    for _ in range(self.num_swap_atoms):
                        line = in_stream.readline()
                        words = line.split()
                        self.swapType1.append(words[0])
                        self.swapType2.append(words[1])
                elif subWord == "transmutate":
                    self.num_transmutate_atoms = int(words[2])
                    self.transmutateFrequency = int(words[3])
                    for _ in range(self.num_transmutate_atoms):
                        line = in_stream.readline()
                        words = line.split()
                        self.mutateType1.append(words[0])
                        self.mutateType2.append(words[1])
                        self.transmuteChemPot.append(float(words[2]))
                elif subWord == "moldyn":
                    self.mdMoveFreq = int(words[2])
            elif keyWord == "symmetry":
                subWord = words[1]
                if subWord == "cubic":
                    self.volMoveSymmetry = 0
                elif subWord == "tetragonal":
                    self.volMoveSymmetry = 1
                elif subWord == "orthorhombic":
                    self.volMoveSymmetry = 2
            elif keyWord == "timestep":
                self.timestep = float(words[1])
            elif keyWord == "mdtemperature":
                self.mdtemperature_K = float(words[1])
            elif keyWord == "mdfriction":
                self.mdfriction = float(words[1])
            elif keyWord == "mdsteps":
                self.mdsteps = int(words[1])
            elif keyWord == "relmethod":
                self.relmethod = words[1]
            elif keyWord == "relsteps":
                self.relsteps = int(words[1])
            elif keyWord == "reltol":
                self.reltol = float(words[1]) 

