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
5.1 control
-----------

The input is broken into sections depending on the functionality. The
main input key words are

+--------------+-------------+------------------------------------------------------+
| **Key word** | **Type**    | **Functionality**                                    |
+==============+=============+======================================================+
| #            |             | comments the line (must be first char on the line)   |                                    |
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
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| freeze       | int         | Freeze the following types of atoms in MC only. E.g. |
|              |             |                                                      |
|              |             | freeze 1                                             |
|              |             | Au                                                   |
|              |             |                                                      |
|              |             | This prevents any movement of Au type atoms          |
| kmcsteps       | int         | The number of kmc steps/cycles                       |
+----------------+-------------+------------------------------------------------------+
| mincap         | double      | The minimum activation energy.                       |
+----------------+-------------+------------------------------------------------------+
| prefactor      | double      | The prefactor used to calculate rates.               |
+----------------+-------------+------------------------------------------------------+
| window         | double      | The window for accepting kmc energies.               |
+----------------+-------------+------------------------------------------------------+
| kmctemperature | double      | The temperature to be used by the kmc simulation     |
+----------------+-------------+------------------------------------------------------+
| kmcevents      | int         | The number of events collected per cycle.            |
+----------------+-------------+------------------------------------------------------+
| basin2delta    | double      | The magnitude of the displacement away from the      |
|                |             | saddle point prior to the attempted relaxation to    |
|                |             | the second basin.                                    |
+----------------+-------------+------------------------------------------------------+
| kmcmethod      | string      | The method used for the saddle search. Only ART or   |
|                |             | dimer are available at the moment                    |
+----------------+-------------+------------------------------------------------------+
| kmcbasinradius | double      | The minimum displacement of an atom before it is     |
|                |             | considered to have entered a new basin               |
+----------------+-------------+------------------------------------------------------+
| recycle        |             | recycles saddle points. Read the code as this is     |
|                |             | experimental and may have significant impact on the  |
|                |             | results                                              |
+----------------+-------------+------------------------------------------------------+
| usegaussian    | double      | The atoms are given a random displacement weighted   |
|                |             | by a Gaussian of the specified width                 |
+----------------+-------------+------------------------------------------------------+

Description of keywords to control either the Activation Relaxation Technique. Again units are specified
in the main directives above.

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| numvectors          | int         | The number of Lnanczos vectors used to obtain        |
|                     |             | eigenvalues. Default 20                              |
+---------------------+-------------+------------------------------------------------------+
| maxeigenvalue       | double      | The max eigenvalue. Once an eigenvalue falls below   |
|                     |             | this value the forces parallel to the eigenvalue are |
|                     |             | used.                                                |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| eigentolerence      | double      | The tolerence to converge the eigenvalues.           |
+---------------------+-------------+------------------------------------------------------+
| lanczosdisplacement | double      | The displacement of atoms used to calculate the      |
|                     |             | eigenvalues from the tri-diagonal matrix             |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | int         | The number of iterations to calculate the saddle     |
|                     |             | point.                                               |
+---------------------+-------------+------------------------------------------------------+
| minmethod           | string      | The mminimisation technique to find the saddle point |
|                     |             | Only FIRE(2) is available at the moment.             |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| damp                | double      | damping factor for the parallel forces               |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

Input keywords for the relaxation

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| debug               |             | increases the amount of information output to files. |
+---------------------+-------------+------------------------------------------------------+
| relaxsteps          | int         | The number of iterations of the minimisaer.          |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | double      | The maximum size of the displacement in FIRE.        |
+---------------------+-------------+------------------------------------------------------+
| method              | string      | The mminimisation technique. There is a choice       |
|                     |             | between FIRE, FIRE2 and ASE the moment.              |
|                     |             | (The latter uses the minimisation technique from the |
|                     |             | librray program ASE.)                        |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

--------------
5.1 potentials
--------------

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| device              | string      | The device wheer the potential energy is calculated  |
|                     |             | default is "cuda". This variable is passed to the    |
|                     |             | ASE/Janus-core calculator                            |
+---------------------+-------------+------------------------------------------------------+
| precision           | string      | either float64 (default) or float32                  |
|                     |             | this string is passed directlt to the calculator     |
|                     |             | i.e. no testing is done on it for max flexibility.   |
+---------------------+-------------+------------------------------------------------------+
| arch                |             | Architecture of the ASE calculator                   |
|                     |             | Default (mace_mp).                                   |
+---------------------+-------------+------------------------------------------------------+
| model               |  string     | The full path and name of the model                  |
+---------------------+-------------+------------------------------------------------------+
| close               |             | stops any further input.                             |
+---------------------+-------------+------------------------------------------------------+
| #                   |             | indicates comment. Must be first character           |
+---------------------+-------------+------------------------------------------------------+
| species             | int         | species keyword followed by the number of different  |
|                     |             | speccies. Each element type should be input as       |
|                     |             | follows:                                             |
|                     |             | name  mass charge atomic_number                      |
+---------------------+-------------+------------------------------------------------------+

The is also the possibility of using machine learned interatomic potentials with ASE calculators (mace is the default). In principal any ASE calculator can be used,
but in practice the Field.py file may need to be modified. A potentials file is always needed to setup the calculation:

   species 2
   Si 28.0 0.0 14
   O 16.0 0.0 8
   model  ~/some/path/mace.model
   close

By default the calculation will run on a single GPU


---------
5.2 basis
---------

The basis format follows that of a simplified extended xyz. The minimum format is::

   number_of_atoms
   Lattice="vectors * 9"
   name      x  y  z
   name      x  y  z

For example::

   2549
   Lattice="32.7232154959 0.0 0.0 0.0 32.7232154959 0.0 0.0 0.0 32.7232154959"
   O 1.4168603868  1.3899782202  1.3108695097
   O 4.1014301991  4.0478192159  1.4171530724
   O 3.7163970637  -1.3764087458  4.0580241564

=============
6. References
=============

