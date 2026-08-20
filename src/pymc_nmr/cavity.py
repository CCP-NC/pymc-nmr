#employs cavity bias style approach to finding suitable positions
import random
import time
import numpy as np
import os
import sys

from datetime import datetime

from pymc_nmr.config import Config
from pymc_nmr.species import Species
from pymc_nmr.field import Field
from pymc_nmr.grid import Grid

# routine that uses cavity bias style to put K close to AL
def main():

    fld = Field()
    spec = Species()

    basisFileName = "basin.xyz"
    fieldFileName = "potentials"
    logFileName = "insert.log"
    
    out_stream = open(logFileName, 'a')

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

    rcut = 3.5
    grd = Grid(20,15,10, 2.0)
    grd.build_grid(basin)
    #occ = grd.get_grid_occupancy()

    for i in range(basin.natoms):
        if basin.symbol[i] != "Al":
            continue

        #get gridpints within rcut of Al
        grd_list = grd.find_empty_grids(basin.pos[i,0], basin.pos[i,1], basin.pos[i,2], basin.vectors, rcut)

        if len(grd_list) == 0:
            print("the value of combination swap distance is too small")
            exit()

        choice = int(len(grd_list) * np.random.random())
        atm = grd_list[choice]
        #basin.pos[atm,:] = grd.grid_pos[atm,:]

        out_stream.write(f"\n K     {grd.grid_pos[atm,0]}   {grd.grid_pos[atm,1]}    {grd.grid_pos[atm,2]}")

    #########################################################################################################
    # start the dimulation
    #########################################################################################################


if __name__ == "__main__":
    main()
