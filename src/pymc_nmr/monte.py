#! /usr/bin/env python3
# 

import random
import time
import numpy as np
import os
import sys

from datetime import datetime

from config import Config
from job_control import JobControl
from species import Species
from statistics import Statistics, TypeStatistics
from field import Field
from basin_hop import BasinHop
from airss_style import AirssStyle
from monte_carlo import MonteCarlo


def main():
    """
    This is the principal function that reads in all the data required (initial config, force field, control) and starts to set the calculation up
    Either a MC, basin hopping and AIRSS-style calculation is initiated
    """
    restart_iteration = 0
    md_time = 0.0

    job = JobControl()
    spec = Species()

    fld = Field()

    bh = BasinHop()
    mc = MonteCarlo()
    ai = AirssStyle()

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
    basin = Config()
    restart_iteration = 0
    restart_time = 0.0 
    restart_energy = 0.0
    try:
        instream = open(basisFileName, "r")
        restart_iteration, restart_time, restart_energy = basin.read_config(instream)
        basin.setup_configuration(spec, out_stream)
    except FileNotFoundError:
        out_stream.write("\n*** could not find configuration file: basin.xyz \n")
        sys.exit(1)

    if job.structure_method == "basinhop":
        job.write_bh_control(out_stream)
        bh.initialise(spec, job, out_stream)
    elif job.structure_method == "airss":
        job.write_airss_control(out_stream)
        ai.initialise(spec, job, out_stream)
    elif job.structure_method == "monte":
        job.write_mc_control(out_stream)
        mc.initialise(spec, job, out_stream)
    else:
        out_stream.write("\n*** unrecognised structure search method \n")
        sys.exit(1)

    #########################################################################################################
    # start the dimulation
    #########################################################################################################

    num_cycles = job.num_cycles

    stats = Statistics()
    type_stats = TypeStatistics()
     
    for cycle in range(num_cycles):
        
        num_steps = 0

        # Start the timer
        start_time = time.time()

        # Do the MC calculation
        initialise = True

        if job.structure_method == "monte":
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Monte Carlo Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")

            mc.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, restart_iteration, out_stream)
        
        elif job.structure_method == "airss":
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" AIRSS Style Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")

            ai.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, initialise, restart_iteration, 
                   restart_energy, out_stream)
        else:
            
            out_stream.write("\n\n" + " *" * 53 + "\n")
            out_stream.write(" Basin Hopping Simulation. Cycle : " + str(cycle) + "\n")
            out_stream.write(" *" * 53 + "\n")

            bh.run(spec, fld, job, stats, type_stats, basin, num_steps, cycle, initialise, restart_iteration, restart_energy, 
                   out_stream)
           

        finish_time = time.time()
        diff = finish_time - start_time

        stats.last_summary(out_stream, num_steps)
        type_stats.last_summary_types(spec, out_stream)

        out_stream.write("\n\n" + " *" * 53 + "\n")
        out_stream.write(f"\n time to Monte Carlo simulation : {diff:.3f} seconds\n")
        out_stream.write("\n *** Monte Carlo Finished. Writing restart\n")
        out_stream.flush()

    out_stream.flush()
    out_stream.close() #end of simulation

if __name__ == "__main__":

    main()
