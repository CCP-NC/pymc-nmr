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
| cutoff              | double      | The short range cutoff for the potential and Ewald . |
+---------------------+-------------+------------------------------------------------------+
| noimage             |             | The nearest image conventin is used by default. This |
|                     |             | keyword uses a slower multiple image method that     |
|                     |             | is better suited to small simulation cells.          |
+---------------------+-------------+------------------------------------------------------+
| noewald             |             | The Ewald sum is used by default for the two-body    |
|                     |             | potential model. Thus this keyword switches it off.  |
+---------------------+-------------+------------------------------------------------------+
| ewald precision     |             | Controls the accuracy of the Ewald sum. Default      |
|                     |             | value: 1.0e-6                                        |
+---------------------+-------------+------------------------------------------------------+
| species             | int         | species keyword followed by the number of different  |
|                     |             | speccies. Each element type should be input as       |
|                     |             | follows:                                             |
|                     |             | name  mass charge atomic_number                      |
+---------------------+-------------+------------------------------------------------------+

As described in the installation section the program can be compiled with either two-body (including Ewald sum),
many-body (metal potentials) or ASE calculator. Note all parameters are in electron volts! 
Parameters compatible with the two-body are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| buck                | *A* , :math:`{\alpha}` , *C*                                       |
+---------------------+--------------------------------------------------------------------+
| morse               | *D* , r\ :sub:`eq` , *k*                                           |
+---------------------+--------------------------------------------------------------------+
| ljones              | :math:`{\eta}` , :math:`{\sigma}`                                  |
+---------------------+--------------------------------------------------------------------+
| bhm                 | *A* , :math:`{\alpha}` , *C* , *D*                                 |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 8.0
   noimage
   species 2
   Mg 24.0 2.0 12
   O 16.0 -2.0  8
   twobody 2
   buck
   Mg O  1428.5 0.2945  0.00
   buck
   O  O 22764.3 0.1490 27.879
   close

Parameters compatible with the metal potentials are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| stch                | :math:`{\eta}` , *a* , *n* , *m* , *c*                             |
+---------------------+--------------------------------------------------------------------+
| gupta               | *A* , r\ :sub:`eq` , *p*, *B*, *q*                                 |
+---------------------+--------------------------------------------------------------------+
| fnsc                | *c0* , *c1* , *c2* , *c* , *A* , *d* , :math:`{\Beta}`             |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 6.5
   species 1
   Al  25.0  0.0 13
   manybody 1 ev
   suttonchen
   Al  Al   0.033147    4.05       7.0        6.0         16.399
   close

The is also the possibility of using machine learned interatomic potentials with ASE calculators and the ASE dimer method. A potentials file is still needed
to setup the calculation:

   species 2
   B 16.0 0.0
   C 16.0 0.0
   model  CoB_v3.model
   close

By default thrdr will run on a single GPU (multiple GPU's has not been tested).

In all calculations a file search.env is required and is employed to switch on or off the ASE functionality. For rigid ion or metal calculations this should be:

   USE_ASE_RELAX = "False"
   USE_ASE_SEARCH = "False"

whilst for MLIP's and ASE dimer the False values should be changed to True:

   USE_ASE_RELAX = "True"
   USE_ASE_SEARCH = "True"

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

.. [1]
   D.T. Gillespie, *J. Phys. Chem.*, 1997, **81**, 2340-2361.

.. [2]
   A.F. Voter, *Phys. Rev. B*, 1986, **34**, 6819-6829.

.. [3]
   C.C. Battaile, *Comput. Methods Appl. Mech. Engrg.*, 2008, **197**,
   3386-3398.

.. [4]
   D. Frenkel and B. Smit, *Understanding Molecular Simulation: From
   Algorithms to Applications*, 2002, Academic Press.

.. [5]
   G. Henkelman and H. Jónsson, *J. Chem. Phys.*, 1999, **111**,
   7010-7020.

.. [6]
   G.T. Barkema and N. Mousseau, *Computational Materials Science*,
   2001, **20**, 285-292.

.. [7]
   R.A. Olsen, G.J. Kroes, and G. Henkelman, *The Journal of Chemical Physics*, 2004,  **121(20)**, 9776–9792.

.. [8]
   E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, P. Gumbsch, *Phys. Rev. Lett.*, 2006, **97** , Article 170201.
