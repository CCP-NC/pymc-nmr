import random
import time
import numpy as np
import os
import sys

from ase import Atoms
from ase.io import write, read, iread
from ase.optimize import BFGS, FIRE, LBFGS

from datetime import datetime

from JobControl import JobControl
from Species import Species
from Statistics import Statistics, TypeStatistics
from Field import Field
from basin_hop import BasinHop

from monte_carlo import MonteCarlo


def main():
    restart_iteration = 0
    md_time = 0.0

    job = JobControl()
    spec = Species()

    fld = Field()

    bh = BasinHop()
    mc = MonteCarlo()

    basisFileName = "basin.xyz"
    fieldFileName = "potentials"
    logFileName = "mc.log"
    jobFileName = "control"

    out_stream = open(logFileName, 'a')

    #job control parameters
    try:
        with open(jobFileName, 'r') as in_stream:
            job.read_job_control(in_stream, out_stream)
    except FileNotFoundError:
        out_stream.write("\n*** could not find control file\n")
        sys.exit(1)

    #field file
    try:
        with open(fieldFileName, 'r') as in_stream:
            fld.readPotential(in_stream, out_stream, spec)
            
    except FileNotFoundError:
        out_stream.write("\n*** could not find potentials file\n")
        sys.exit(1)

    # Read data from basis file in xyz format
    try:
        basins = read(basisFileName, index=":")
        print("the number of datsets read in", len(basins))
    except FileNotFoundError:
        out_stream.write("\n*** could not find configuration file: basin.xyz \n")
        sys.exit(1)

    if job.relax_structure:
        bh.initialise(spec, job, out_stream)
    else:
        mc.initialise(spec, job, out_stream)

    #for ib in range(job.num_boxes):
    #    basins[ib].freeze_atom_types(spec, job.frozen_types)

    #########################################################################################################
    # start the dimulation
    #########################################################################################################

    num_cycles = job.num_cycles

    stats = []
    type_stats = []
     
    for _ in range(num_cycles):
        stats.append(Statistics())
        type_stats.append(TypeStatistics())

    for cycle in range(num_cycles):

        
        
        num_steps = 0

        # Start the timer
        start_time = time.time()

        # Do the MC calculation
        initialise = True

        if job.relax_structure:
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Basin Hopping Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")

            bh.run(spec, fld, job, stats, type_stats, basins, num_steps, cycle, initialise, out_stream)
        else:
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Monte Carlo Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")

            mc.run(spec, fld, job, stats, type_stats, basins, num_steps, cycle, initialise, out_stream)

        finish_time = time.time()
        diff = finish_time - start_time

        stats[cycle].last_summary(out_stream, num_steps)
        type_stats[cycle].last_summary_types(spec, out_stream)

        out_stream.write("\n\n" + " *" * 53 + "\n")
        out_stream.write(f"\n time to Monte Carlo simulation : {diff:.3f} seconds\n")
        out_stream.write("\n *** Monte Carlo Finished. Writing restart\n")
        out_stream.flush()

        #with open("restart", "w") as restart_stream:
        #    for ib in range(job.num_boxes):
        #        basins[cycle].dump_basis(spec, md_time, 0.0, job.mc_steps, restart_stream)

    out_stream.flush()
    out_stream.close() #end of simulation

if __name__ == "__main__":
    main()
