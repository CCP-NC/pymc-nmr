#calculates weights for importance sampling of nmr spectra. Follows equation 30 in 

import numpy as np

from ase import Atoms
from ase.io import read, write
from ase.build.tools import sort

def normalise(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def get_weights(expon):
    sum_exp = np.sum(expon)
    print("sum ", sum_exp)
    return expon / sum_exp

temperature = 300.0
BOLTZMANN = .00008617333262145 # in  eV
beta = 1.0 / (temperature * BOLTZMANN)

#read in all the frames in a list of atoms objects
input_name = "accepted.xyz"
dataset = read(input_name, format="extxyz", index=":")

nframes = len(dataset)

energies = np.zeros(nframes, dtype=np.float64)
expon = np.zeros(nframes, dtype=np.float64)
weights = np.zeros(nframes, dtype=np.float64)

for i in range(nframes):
    energies[i] = dataset[i].get_potential_energy() 

min_eng = np.min(energies)
print("the minimum energy ", min_eng)

for i in range(nframes):
    #expon[i] = np.exp(-(energies[i] * beta), dtype=np.float64)
    expon[i] = np.exp(-((energies[i]-min_eng) * beta), dtype=np.float64)

#calculate the weights for the spectra
weights = get_weights(expon)

out_io = open("weights", "w")
for i in range(nframes):
    if weights[i] > 1.0e-6:
        out_io.write(f" {i}  {energies[i]}  {energies[i]-min_eng}   {expon[i]}   {weights[i]} \n")

print("sum of weights ", np.sum(weights))