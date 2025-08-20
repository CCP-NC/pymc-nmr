import numpy as np

class Energy:
    def __init__(self, src=None):
        if src:
            self.totalEnergy = src.totalEnergy
        else:
            self.totalEnergy: np.float64 = 0.0
            

    def __sub__(self, src):
        result = Energy()
        result.totalEnergy = self.totalEnergy - src.totalEnergy
        
        return result

    def __add__(self, src):
        result = Energy()
        result.totalEnergy = self.totalEnergy + src.totalEnergy
        
        return result

    def __iadd__(self, src):
        self.totalEnergy += src.totalEnergy
        
        return self

    def zero(self):
        self.totalEnergy = 0.0
        

    def print_energy(self, box, outStream):
        
        outStream.write(f"\n\n energies of box {box}\n")
        outStream.write(f" total (internal) energy      {self.totalEnergy:25.15e}\n")
        

    def get_total_energy(self) -> np.float64:
        return self.totalEnergy