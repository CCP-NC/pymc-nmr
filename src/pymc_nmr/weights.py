"""Boltzmann population weights for an ensemble of accepted structures.

Reads an extended-xyz trajectory of structures with calculated potential
energies (normally ``accepted.xyz``), finds the minimum energy, and converts
the relative energies into normalised populations

    w_i = exp(-(E_i - E_min) / kB T) / sum_j exp(-(E_j - E_min) / kB T)

following equation 30 of Mosquera-Lois, Kavanagh, Klarbring, Tolborg and Walsh,
*Chem. Soc. Rev.* **2023**, 52, 5812. These weights are intended for averaging
computed NMR parameters over the ensemble.

Energies are assumed to be in eV (the ASE convention used throughout pymc-nmr).

.. warning::

   **Known limitation: the weighting semantics of ``accepted.xyz`` are
   ambiguous, and this script does not resolve them.**

   ``accepted.xyz`` is appended to on every Metropolis-accepted move, so it is
   a Markov chain already distributed according to ``exp(-E/kB T_run)`` at the
   temperature of the run. Applying a Boltzmann factor to those same structures
   therefore double-counts the energy dependence, giving an effective
   ``exp(-2 E / kB T)`` if the weighting temperature equals the run
   temperature.

   Two self-consistent interpretations exist, and you must decide which one
   your analysis assumes:

   1. *Chain semantics.* Treat ``accepted.xyz`` as an MC trajectory. Ensemble
      averages should then use uniform weights ``1/N``, and this script should
      not be used.
   2. *Pool semantics.* Treat ``accepted.xyz`` as a pool of candidate
      orderings found by the search. Boltzmann reweighting is then appropriate,
      but the pool must first be reduced to symmetry-distinct structures and
      each weighted by its configurational multiplicity. This script performs
      neither the deduplication nor the multiplicity weighting.

   No deduplication is applied here, so repeated visits to the same basin are
   counted once per acceptance. Resolve this before quoting ensemble-averaged
   NMR parameters.
"""

import argparse

import numpy as np
from ase.io import read
from ase.units import kB


def get_weights(expon):
    """Normalise a vector of Boltzmann factors to sum to one."""
    return expon / np.sum(expon)


def main(args=None):
    parser = argparse.ArgumentParser(
        description=(
            "Calculate normalised Boltzmann population weights from a "
            "trajectory of structures with calculated energies (eV). "
            "See the module docstring for an important caveat on how "
            "accepted.xyz should be interpreted."
        )
    )
    parser.add_argument(
        "--input",
        help="input trajectory of configurations with energies.",
        type=str,
        default="accepted.xyz",
    )
    parser.add_argument(
        "--format",
        help="ASE format of the input file.",
        type=str,
        default="extxyz",
    )
    parser.add_argument(
        "--output",
        help="output file name for the table of weights.",
        type=str,
        default="weights",
    )
    parser.add_argument(
        "--temperature",
        help="temperature in K used for the weight calculation.",
        type=np.float64,
        default=300,
    )

    args = parser.parse_args(args=args)

    input_name = args.input
    output_name = args.output
    input_format = args.format
    temperature = args.temperature

    if temperature <= 0.0:
        parser.error("--temperature must be greater than zero.")

    beta = 1.0 / (temperature * kB)

    # read in all the frames as a list of atoms objects
    dataset = read(input_name, format=input_format, index=":")
    nframes = len(dataset)

    energies = np.zeros(nframes, dtype=np.float64)
    for i in range(nframes):
        energies[i] = dataset[i].get_potential_energy()

    min_eng = np.min(energies)
    print("the minimum energy ", min_eng)

    # Subtracting the minimum energy before exponentiating is numerically
    # essential: exp(-E*beta) for raw MLIP energies (thousands of eV) underflows.
    expon = np.exp(-(energies - min_eng) * beta)

    # calculate the weights for the spectra
    weights = get_weights(expon)

    # Rows with negligible weight are omitted, so the printed column does not
    # sum to exactly one; the reported total below is over all frames.
    with open(output_name, "w") as out_io:
        out_io.write("# frame  energy_eV  relative_energy_eV  boltzmann_factor  weight\n")
        for i in range(nframes):
            if weights[i] > 1.0e-6:
                out_io.write(
                    f" {i}  {energies[i]}  {energies[i] - min_eng}   "
                    f"{expon[i]}   {weights[i]} \n"
                )

    print("sum of weights ", np.sum(weights))


if __name__ == "__main__":
    main()
