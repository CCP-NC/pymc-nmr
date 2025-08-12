# read in an ASE database of multiple images with calculated potential energies (the input file must be accepted.xyz).
# Follows equation 30 in Irea Mosquera-Lois, Sean R. Kavanagh, Johan Klarbring, Kasper Tolborg and Aron Walsh Chem Soc Rev 2023, vol 52, 5812.
# Once the energies are read in the minimum energy is found and the "defect" energy calculated. The latter energy is used to determine the "population" and weight
# for use within the calculation of the spectra.
# NB assumes energies ate in eV.  

import numpy as np
import argparse

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
    #print("sum ", sum_exp)
    return expon / sum_exp

def main(args=None):
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """add otcar relaxation to extended xyz for training."""
        )
    )
    parser.add_argument(
        "--input",
        help="configurations + energies. ",
        type=str,
        default='accepted.xyz',
    )
    parser.add_argument(
        "--format",
        help="input style. ",
        type=str,
        default='extxyz',
    )
    parser.add_argument(
        "--output",
        help="output xyz file name. ",
        type=str,
        default="weights",
    )
    parser.add_argument(
        "--temperature",
        help="the temperature used for weight calculation.",
        type=np.float64,
        default=300,
    )

    # Parse the args
    args = parser.parse_args(args=args)


    input_name = args.input
    output_name = args.output
    format = args.format
    temperature = args.temperature

    BOLTZMANN = .00008617333262145 # in  eV
    beta = 1.0 / (temperature * BOLTZMANN)

    #read in all the frames in a list of atoms objects
    dataset = read(input_name, format=format, index=":")

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

    out_io = open(output_name, "w")
    for i in range(nframes):
        if weights[i] > 1.0e-6:
            out_io.write(f" {i}  {energies[i]}  {energies[i]-min_eng}   {expon[i]}   {weights[i]} \n")

    print("sum of weights ", np.sum(weights))


if __name__ == "__main__": 
    main()
