#employs cavity bias style approach to finding suitable positions
import random
import time
import numpy as np
import os
import sys

from datetime import datetime

from config import Config
from Species import Species
from Field import Field
from grid import Grid


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

    
    grd = Grid(5,5,5, 2.0)
    grd.build_grid(basin)
    #########################################################################################################
    # start the dimulation
    #########################################################################################################


if __name__ == "__main__":
    main()
