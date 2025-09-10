===============
1. Introduction
===============

Pymc-nmr has three modes of operation that are based on swap moves. A swap move is where the positions of two atoms, chosen at random, are interchanged.
For example, the two atoms outlined in red in the figure below.

.. image:: swap-move.png
  :width: 400

The difference in energy of the old (*o*) and new (*n*) cells is calculated and employed within the Metropolis-Hastings sampling. Thus the move is accepted with
the probability

:math:`acc(o \rightarrow n) = min(1, exp(- \beta [E_{n} - E_{o}]))`

---------------
1.1 Monte Carlo
---------------

This is the standard Monte Carlo method described in the books by (\ :ref:`Allen and Tildesley 1989 <ref-Allen_Tildesley_87>`\ ; \ :ref:`Frenkel and Smit 2002 <ref-Frenkel_Smit>`\ ).
DL_MONTE. Simulations can be performed in the canonical (:math:`NVT`) ensemble and isobaric-isothermal (:math:`NPT`) (where MC moves attempting variations in 
the volume of the system are applied). In addition, the swap moves described above can be employed.

Note as a *single* particle is moved
at a time it tends to be slow and is only useful where the difference in energy is small. However, using the small MACE_mp model the length of the simulations are more acceptable
and the method is capable of getting you out of a jam if the forces in your structure are very high.

-----------------
1.2 Basin Hopping
-----------------

The basin hopping approach is useful when the difference in energy between swapped particles is large. Each swap of particles is followed by an energy minimisation. The 
difference in energy is taken after the minimisation and the Metropolis selection rule applied. Although the calculation of the energy consumes more resources
the probability of the swap being accepted is significantly greater. 
In some cases the position of a particle is not known and we have created the possibility of using a NVT molecular dynamics calculation followed by energy minimisation
to further reduce the total energy of the cell.

----------------------------
1.2 Randomised Configuration
----------------------------

The last method is an ab initio random structure search method. In this approach all the particles of two specific type are swapped and an energy minimisation performed.
In out limited experience, this the basin hooping method is more efficient than this technique.

===============
2. Installation
===============

The easiest way to install the pymc-nmr package is to use the uv package manager. the installation instructions and guides can be found `here`_ .

.. _here: https://docs.astral.sh/uv/getting-started/

The program can be installed with ::
   uv pip install -r pyproject.toml 

It is also possible to set the program up manually using pip. However it is difficult to define how to setup the program as many computer systems are slightly different, especially clusters. However, the program is written in python and some guidelines are as follows.
It is likely that you need to create a python environment following the python `documentation`_ . 

.. _documentation: https://docs.python.org/3/library/venv.html

External libraries required are `NUMPY`_, `MACE`_, `ASE`_, `JANUS`_,  `DFTD3`_ and  `PYTEST`_ (and their dependecies). In addition libraries to satisfy the ASE calculator 
or MLIP potentials may be necessary. (if you want to run the tests then pytest also needs to be installed).

.. _NUMPY: https://numpy.org/
.. _MACE: https://github.com/ACEsuit/mace
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _JANUS: https://github.com/stfc/janus-core
.. _DFTD3: https://github.com/dftd3
.. _PYTEST: https://docs.pytest.org/en/stable/

Once this is complete the program can be run using the command *python3 /path/to/prog/monte.py*.

A user guide exists in the docs folder. To creat a set of html pages run ::

   sphinx-build -M html ./source ./build

Tests are in the tests folder abd are in pytest format.

=========================
3. How to run the program
=========================


The program requires three input files, *control*, *potentials* and *basin.xyz* (extended xyz) that contain the
keywords for the functionality of the program, potential model and the atomic
positions respectively.

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
| writestats   |             | Activates writing of energy/cell data to stats file  |
+--------------+-------------+------------------------------------------------------+
| statsfreq    | int         | frequency of writing data to stats file              |
+--------------+-------------+------------------------------------------------------+
| grid         | 3*int,float | the spacing of the combiswap grid and the cutoff     |
|              |             | employed for the exclusion of regions                |
|              |             | Default values are 10 10 10 and 2.0                  |
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
| conp           | string      | cell anngles and lengths can change during the       |
|                |             | relaxation (default)                                 |
+----------------+-------------+------------------------------------------------------+
| cona           | string      | only the cell lengths can change during the          |
|                |             | relaxation. The angles are fixed.                    |
+----------------+-------------+------------------------------------------------------+
| conv           | string      | the cell lengths and angles are fixed.               |
+----------------+-------------+------------------------------------------------------+

The keywords that are specific to MC simulations only:

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| maxdistance    | float       | the maximum atomic displacement. The default is      |
|                |             | 0.001                                                |
+----------------+-------------+------------------------------------------------------+
| distanceupdate | int         | The number of MC steps between the updates of the    |
|                |             | max displacement to acieve the target ratio below    |
+----------------+-------------+------------------------------------------------------+
| distanceratio  | float       | The target ratio for the acceptance / rejection of   |
|                |             | particle displacements                               |
+----------------+-------------+------------------------------------------------------+
| maxvolume      | float       | the maximumvolume displacement. The default is       |
|                |             | 0.001                                                |
+----------------+-------------+------------------------------------------------------+
| volupdate      | int         | The number of MC steps between the updates of the    |
|                |             | max volume displacement to acieve the target ratio   |
|                |             | below                                                |
+----------------+-------------+------------------------------------------------------+
| volratio       | float       | The target ratio for the acceptance / rejection of   |
|                |             | volume expansion / contraction                       |
+----------------+-------------+------------------------------------------------------+

These keywords follow the *move* command in the control file.

+----------------+--------------------------------------------------------------------+
| **Key word**   |  **Functionality**                                                 |
+================+====================================================================+
| atoms          | The *atoms* keyword must be followed the number of atom types to   |
|                | be moved and the probability of the move occuring. It is only      |
|                | active in Monte Carlo simulations. An example is                   |
+----------------+--------------------------------------------------------------------+
| swap           | Swaps interchange the positions of two types. The number of atom   |
|                | pairs and the probability of the move occuring. It is              |
|                | active in all types of calculation.                                |
+----------------+--------------------------------------------------------------------+
| combiswap      | This is the same as the standard swap except that an idditional    |
|                | atom is moved in conjunction with the second atom and a grid is    |
|                | employed to facilitate the location of this atom. It is only       |
|                | active in basin hopping.  See tutorial 1 for an example.           |
+----------------+--------------------------------------------------------------------+
| moldyn         | This activates an NVT molecular dynamics calculation (only within  |
|                | basin hopping). An integer gives the probability of this ooccuring.|
+----------------+--------------------------------------------------------------------+
| volume         | A volume move within the monte Carlo method only An integer gives  | 
|                | the probability of this ooccuring.                                 | 
+----------------+--------------------------------------------------------------------+

An example ov *move atom* ::

   move atoms 2 100                                                  
   Si                                                                 
   O        
   
An example of *swap* ::

   move swap 2 100                                                    
   Si  Al                                                             
   K   ghost   (the fictitious particle must be called ghost) 

How to use molecular dynamics in basin hopping ::

   move moldyn 20

Using the volume move in Monte Carlo ::

   move volume 20


Example number 1: Monte Carlo control file

--------------

::

   # Commented out lines start with '#
   print 1                           # the frequency for printing energies etc
   maxdistance 0.35                  # the maximum displacement of particles
   maxvolume 0.035                   # the maximum 
   steps 30000                       # the number of steps
   equilsteps 15000                  # the number of steps before average properties are calculated
   archive                           # flag to create an archive/history file (archive.xyz) of the configurations
   archivefrequency 1000             # frequency for writing configurations to the archive
   method monte                      # use standard MC method
   check 5                           # check that the running total energy and a freshly calculated value agree also restart.xyz is written at this point 

   temperature 1300                  # temperature / Kelvin
   # displacement MC moves for atoms
   move atoms 4 50
   Na
   O
   Si
   Al
   #particle swap moves for zeolite containing Si and Al
   move swap 1 30
   Al Si
   # volume moves in NpT
   move volume 20  
   close                             # stop input 

--------------

Example number 2: Basin hopping control file

--------------

::

   print 1                           # the frequency for printing energies etc
   steps 10                          # the number of steps
   equilsteps 5                      # the number of equilibration steps
   archive                           # flag to create an archive/history file (archive.xyz) of the configurations
   archivefrequency 10               # frequency for writing configurations to the archive
   savedownhill                      # downhill (energy less than previous configuration) moves are saved to downhill.xyz
   temperature 1200                  # tempretaure (Kelvin) of Boltzmann sampling
   method basinhop                   # Basin hopping method selected
   relmethod bfgs                    # relaxation/minimisation method
   relsteps 1500                     # maximum number of steps allowed to complete the relaxation
   reltol 1.0e-2                     # tolerance for the relaxation
   mdsteps 500                       # the number of steps in the MD
   mdtemperature 200                 # the temperature of the MD mini-simulation
   move swap 1 80                    # the atoms to be swapped
   Al Si
   move moldyn 20                    # also usee MD to shake the system down
   close                             # stop input

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
| dispersion     |             | switches on dispersion (default: False)              |
+----------------+-------------+------------------------------------------------------+

An example of a potentials file for STA30 containing both H and K ::

   species 5                                                 #the number of different types of atoms
   H 0.0 0.0 1                                               #Chemical symbol followed by 2 floats and an int
   K 0.0 0.0 1                                               #they are redundant at the moment
   Al 0.0 0.0 13
   Si 0.0 0.0 14
   O 0.0 0.0 8
   model ../../../models/MACE-matpes-r2scan-omat-ft.model    #the path to the MLIP model
   #arch mace_mp                                             # architecture
   #device gpu                                               # will run on the gpu
   #precision float 64                                       # use double precision
   close

===============
4. Output Files
===============

All styles of calculations create a file called *mc.log* that provides the details of the simulation. Initially the input is written to this file followed by the
details of the simulation. Subsequently the total energy and the cell properties (*a, b, c, alpha, beta and gamma*) are printed. Finally the number of successful
swaps (and other moves if selected) are given and lastly the average energy and cell properties along side the fluctuations.

The configurations are written to the following files (all use extended xyz)
- archive.xyz (MC and basin hopping) if activated and periodically dump the structure to the file. It is overwritten unless the *restart* keyword is used.
- restart.xyz (all modes) the final configuration and should be copied to basin.xyz to continue the calculation
- accepted.xyz (basin hopping). If a configuration is accepted, even by Metropolis sampling, then the configuration is written (**NB** this file is lways appended
to so if you dont want the data from the previous run then you should delete it!)
- downhill.xyz (basin hopping). if a configuration has a lower energy than the previous configuration then it will be saved. Be aware that the previous accepted structure may not be the lowest energy configuration due to Metropolis-Hastings sampling. The file is also appended to.

==========
5. Weights
==========

Included within the python source folder is a script called weights.py which can be used to read in a set of ASE images (Atoms objects) that have precalculated potential energies (in eV).
The single file of configurations can be specified with the --input option (default: accepted.xyz), the format of the configurations with --format (default: extxyz), output with --output
(default: weights) and temperature with --temperature (default: 300 K). The weights or population of "defects" for the calculation of NMR spectra are given by the equation

.. math::
   w_{i} = \frac{exp(\Delta E) / k_{B} T} {\sum exp(\Delta E) / k_{B} T}

where :math:`\Delta E` is the difference in energy between the convex hull minimum and the energy of the configuration. The convex hull minimum is taken as the lowest energy in 
configuration in the input file.

The output file has five columns. These are the configuration index, potential energy, the energy difference (potential energy minus the potential energy of the minimum), the exponential of this
and the weight (population) respectively.

=============
6. References
=============

.. container:: references csl-bib-body hanging-indent   :name: refs

   .. container:: csl-entry
      :name: ref-Allen_Tildesley_87

      Allen MP, and Tildesley DJ. 1989. *Computer Simulation of Liquids*. Oxford University Press, USA.

   .. container:: csl-entry 
      :name: ref-Frenkel_Smit

      Daan Frenkel, and Berend Smit. 2002. *Understanding Molecular Simulation: From Algorithms to Applications*. Second. San Diego: Academic Press.

   

