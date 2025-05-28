===============
1. Introduction
===============

Some pre-amble about the methods in the 

===============
2. Installation
===============

External libraries required are NUMPY, MACE and ASE (and their dependecies). In addition libraries to satisfy the ASE calculator 
or MLIP potentials may be necessary.


=========================
3. How to run the program
=========================


The program requires three input files, *control*, *potentials* and *basin.xyz* that contain the
keywords for the functionality of the program, potential model and the atomic
positions respectively (in extended xyz).

-----------
3.1 control
-----------

The input is broken into sections depending on the functionality. The
main input key words are

+--------------+-------------+------------------------------------------------------+
| **Key word** | **Type**    | **Functionality**                                    |
+==============+=============+======================================================+
| #            |             | comments the line (must be first char on the line)   |                                    
+--------------+-------------+------------------------------------------------------+
| print        | int         | The frequency of printing energies to mc.log         |
|              |             | defaults to 1 i.e. every iteration                   |
+--------------+-------------+------------------------------------------------------+
| check        | int         | The frequency at which a sanity check is undertaken  |
|              |             | defaults to 1000                                     |
+--------------+-------------+------------------------------------------------------+
| restart      |             | Restart the calculation with old configuration       |
|              |             | make sure to cp restart.xyz to basin.xyz             |
+--------------+-------------+------------------------------------------------------+
| temperature  | float       | The temperature (K) of simulation used in Boltzman   |
|              |             | sampling                                             |
+--------------+-------------+------------------------------------------------------+
| method       | string      | The type of simulation required. Use either          |
|              |             | basinhop (default), airss or monte. The latter       |
|              |             | is depricated at the moment as it does not appear    |
|              |             | to work for swapping Si/Al                           |
+--------------+-------------+------------------------------------------------------+
| close        |             | Finish all input                                     |
+--------------+-------------+------------------------------------------------------+
| archive      |             | Turns on the creation of a archive/history file      |
+--------------+-------------+------------------------------------------------------+
| archivefreq  | int         | The frequency at which the configuration is saved    |
|              |             | defaults to 1000                                     |
+--------------+-------------+------------------------------------------------------+
| equilsteps   | int         | The number of equilibration steps                    |
+--------------+-------------+------------------------------------------------------+
| savedownhill |             | During basin hopping method any configuration that   |
|              |             | has an energy lower than the previous configuration  |
|              |             | will be saved (NB there may have been an up hill move|
|              |             | between downhill events).                            |
+--------------+-------------+------------------------------------------------------+
| move         | string      | The move keyword is used to provide an action on the |
|              |             | configuration by the means of a second key word.     |
|              |             | Each keyword is described below                      |
+--------------+-------------+------------------------------------------------------+

Description of keywords to control the molecular dynamics functionality. 

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| mdtimestep     | float       | The time step in femtoseconds (default 2.0)          |
+----------------+-------------+------------------------------------------------------+
| mdtemperature  | float       | The temperature for the MD in K (default 1000 K)     |
+----------------+-------------+------------------------------------------------------+
| mdfirction     | float       | NVT friction (default 0.01 / fs)                     |
+----------------+-------------+------------------------------------------------------+
| mdsteps        | int         | The number of MD steps (default 1000)                |
+----------------+-------------+------------------------------------------------------+

Description of keywords to control the relaxation in basin hopping or AIRSS style calculations. 

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| relmethod      | string      | the type of relaxation method. Options are lbfgs     |
|                |             | (default), bfgs, or fire                             |
+----------------+-------------+------------------------------------------------------+
| reltol         | float       | tolerence for convergence                            |
+----------------+-------------+------------------------------------------------------+
| relsteps       | int         | The number of relaxation steps (default 1000)        |
+----------------+-------------+------------------------------------------------------+


+----------------+-------------+------------------------------------------------------+
| **Key word**   |  **Functionality**                                    |
+================+=============+======================================================+
| atoms          | The *atoms* keyword must be followed the number of atom types to   |
|                | be moved and the probability of the move occuring. It is only      |
|                | active in Monte Carlo simulations. An example is                   |
|                | move atoms 2 100                                                   |
|                | Si                                                                 |
|                | O                                                                  |
+----------------+-------------+------------------------------------------------------+
| swap           | Swaps interchange the positions of two types. The number of atom   |
|                | pairs and the probability of the move occuring. It is only         |
|                | active in all types of calculation. An example is                  |
|                | move swap 2 100                                                    |
|                | Si  Al                                                             |
|                | K   ghost   (the fictitious particle must be called ghost)         |
+----------------+-------------+------------------------------------------------------+
| moldyn         | This activates an NVT molecular dynamics calculation (only within  |
|                | basin hopping). An integer gives the probability of this ooccuring.|
+----------------+-------------+------------------------------------------------------+
| volume         | A volume move within the monte Carlo method only An integer gives  | 
|                | the probability of this ooccuring.                                 | 
+----------------+-------------+------------------------------------------------------+


---------
3.2 basis
---------

The basis format follows that of a simplified extended xyz. The minimum format is::

   number_of_atoms
   Lattice="vectors * 9"
   name      x  y  z
   name      x  y  z

*NB* any fictitious or ghost atoms must be placed last in the basin.xyz file.

For example::

   2549
   Lattice="32.7232154959 0.0 0.0 0.0 32.7232154959 0.0 0.0 0.0 32.7232154959"
   O 1.4168603868  1.3899782202  1.3108695097
   O 4.1014301991  4.0478192159  1.4171530724
   O 3.7163970637  -1.3764087458  4.0580241564

--------------
3.3 potentials
--------------

Description of keywords used within the potentials file. This is consumed by the field.py 
class and controls the calculation style. 

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| species        | int         | the number of different types. Each species are on   |
|                |             | following lines e.g.                                 |
|                |             | species 1                                            |
|                |             | O 16.0 0.0 8.0 (the floats are required but not used |
|                |             | at present)                                          |
+----------------+-------------+------------------------------------------------------+
| close          |             | Finish all input                                     |
+----------------+-------------+------------------------------------------------------+
| #              |             | comments the line (must be first char on the line)   |                                    
+----------------+-------------+------------------------------------------------------+
| device         | string      | where the energy evaluation should be evaluated      |
|                |             | it is not checked by the program but sent directly   |
|                |             | the ASE/janus calculator. (default: "cuda")          |
+----------------+-------------+------------------------------------------------------+
| precision      | string      | The precision of the calculation.                    |
|                |             | It is not checked by the program but sent directly   |
|                |             | the ASE/janus calculator. (default: "float64")       |
+----------------+-------------+------------------------------------------------------+
| arch           | string      | The style of calculation used in the calculator.     |
|                |             | In theory any ASE calculator could be used but it    |
|                |             | could require the Field.py to be modified. Therefore |
|                |             | it is not checked by the program but sent directly   |
|                |             | the ASE/janus calculator. (default: "mace_mp")       |
+----------------+-------------+------------------------------------------------------+
| model          | string      | The full path to the MLIP model is required          |
+----------------+-------------+------------------------------------------------------+


=============
6. References
=============

