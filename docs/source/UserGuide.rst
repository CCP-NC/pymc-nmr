=============================
PyMC-NMR: Modern User Guide
=============================

.. note::

   New to PyMC-NMR? Start with the :doc:`quickstart`, which walks through
   installation and a first simulation. This page is the complete reference.

PyMC-NMR has three modes of operation that are based on swap moves. A swap move is where the positions of two atoms, chosen at random, are interchanged.
The difference in energy of the old (*o*) and new (*n*) cells is calculated and employed within the Metropolis-Hastings sampling. Thus the move is accepted with the probability:

:math:`acc(o \rightarrow n) = min(1, exp(- \beta [E_{n} - E_{o}]))`

where :math:`\beta = 1 / k_{B}T` and :math:`T` is ``job.temperature``.

.. figure:: swap-move.png
   :width: 60%
   :align: center

   A swap move: two atoms of different types exchange positions.

--------------------
Hardware: use a GPU
--------------------

Runtime is dominated by evaluations of the machine-learned potential. Each
basin-hopping or AIRSS step performs a geometry relaxation, and each relaxation
step is a further evaluation, so one move costs hundreds of potential calls.

**A GPU is typically one to two orders of magnitude faster than a CPU here**,
and for anything beyond a small test you should use one. Set the device in the
``potential`` block:

.. code-block:: yaml

   potential:
     device: cuda      # NVIDIA GPU; "mps" on Apple Silicon; "cpu" otherwise

Verify that PyTorch can see the device before launching a long job, since a
CPU-only PyTorch build will fail on ``device: cuda``:

.. code-block:: bash

   python -c "import torch; print(torch.cuda.is_available())"

When running on CPU, set ``potential.threads_per_rank`` to the number of
physical cores actually available. The default of ``1`` is deliberately
conservative, and setting it higher than your core count (or higher than
``cores / MPI ranks``) causes oversubscription and runs *slower*.

-----
Units
-----

PyMC-NMR follows ASE conventions throughout. Energies are in **eV**, lengths in
**Angstrom**, temperatures in **Kelvin**, and consequently pressure is in
**eV/Angstrom³** — *not* GPa. Times given to the MD shake move are in
femtoseconds and converted internally.

------------------
Modes of Operation
------------------

1. **Monte Carlo**: standard canonical (:math:`NVT`) or isobaric-isothermal (:math:`NPT`) simulations using simple particle translation or volume moves.
2. **Basin Hopping**: a powerful approach where each swap is followed by geometry relaxation/minimisation, resulting in higher swap acceptance rates.
3. **AIRSS Style**: ab initio random structure search where all particles of specific types are randomized, followed by a geometry relaxation.

=============
Installation
=============

PyMC-NMR is designed for modern Python virtual environments.

.. code-block:: bash

   git clone https://github.com/CCP-NC/pymc-nmr.git
   cd pymc-nmr
   uv venv
   source .venv/bin/activate
   uv pip install -e .

Optional extras are available for the test suite, the documentation build, the
example notebooks, and MPI:

.. code-block:: bash

   uv pip install -e ".[test]"        # pytest
   uv pip install -e ".[docs]"        # sphinx
   uv pip install -e ".[notebooks]"   # jupyter
   uv pip install -e ".[mpi]"         # mpi4py, for parallel replicas

The core dependencies are ``numpy``, ``ase``, ``janus-core``, ``mace-torch``,
``torch``, ``pydantic``, ``pyyaml``, ``typer`` and ``rich``. See the
:doc:`quickstart` for a step-by-step version of this.

========================
How to Run the Program
========================

The program is run using a Command Line Interface. It requires only two files:

1. **Initial Structure File** (e.g., ``basin.xyz``): An ASE-compatible file with positions and unit cell.
2. **Configuration File** (e.g., ``config.yaml``): A unified YAML file describing simulation parameters and model configurations.

To launch a simulation:

.. code-block:: bash

   pymc-nmr --config config.yaml --structure basin.xyz

All arguments have default fallbacks to ``config.yaml`` and ``basin.xyz`` in the current directory, so you can simply run:

.. code-block:: bash

   pymc-nmr

To display help and available CLI arguments:

.. code-block:: bash

   pymc-nmr --help

========================
Configuration Reference
========================

The entire simulation is controlled via a single unified ``config.yaml`` file. Pydantic automatically validates its contents prior to execution.

Validation is strict, so that configuration mistakes surface immediately rather
than as a silently wrong simulation:

* **Unknown keys are rejected.** A misspelled key raises ``Extra inputs are not
  permitted`` instead of being ignored and leaving the default in force.
* **Enumerated options are checked.** ``structure_method``, ``relax.method``,
  ``relax.style``, ``vol_move_symmetry``, ``device`` and ``precision`` only
  accept the values listed in the table below.
* **Only two settings are mandatory**: ``potential.model`` and
  ``job.temperature``.

You do not list the chemical species anywhere: they are derived automatically
from the structure file.

An example of a complete ``config.yaml`` file:

.. code-block:: yaml

   job:
     structure_path: basin.xyz
     structure_method: basinhop      # Options: basinhop, airss, monte
     temperature: 1200.0             # Simulation temperature (Kelvin)
     dump_archive: true              # Periodic archiving
     archive_frequency: 20           # Archive interval
     save_downhill: true             # Save downhill configurations

   potential:
     arch: mace_mp                   # ASE calculator type
     model: MACE-matpes-r2scan-omat-ft.model
     device: cpu                     # cpu, cuda, or mps
     precision: float64              # float64 or float32

   monte:
     steps: 200                      # Number of steps

   moves:
     swaps:
       - pairs:
           - [Al, Si]
         counterion: H               # Optional counterion to swap/relocate
         max_distance: 4.0           # Optional radius for relocating H
         frequency: 100

   relax:
     method: bfgs                    # Options: lbfgs, bfgs, fire
     steps: 500                      # Max minimisation steps
     tol: 0.01                       # Minimisation tolerance

---------------------------
Configuration Keyword Table
---------------------------

+--------------------------------+----------------+-----------------------------------------------------------+
| Keyword                        | Default        | Description                                               |
+================================+================+===========================================================+
| **job**                                                                                                     |
+--------------------------------+----------------+-----------------------------------------------------------+
| structure_path                 | basin.xyz      | Path to the initial structure file (ASE-compatible XYZ)   |
+--------------------------------+----------------+-----------------------------------------------------------+
| structure_method               | basinhop       | Simulation mode: ``basinhop``, ``airss``, or ``monte``    |
+--------------------------------+----------------+-----------------------------------------------------------+
| num_cycles                     | 1              | Number of independent simulation cycles to run            |
+--------------------------------+----------------+-----------------------------------------------------------+
| pressure                       | 0.0            | External pressure for NpT volume moves (eV/Ang^3, not GPa)|
+--------------------------------+----------------+-----------------------------------------------------------+
| temperature                    | *required*     | Sampling temperature (K); sets beta = 1/kT. No default    |
+--------------------------------+----------------+-----------------------------------------------------------+
| wrap                           | false          | Wrap atoms into the unit cell after each move             |
+--------------------------------+----------------+-----------------------------------------------------------+
| restart                        | false          | Restart from ``restart.xyz`` copied to ``basin.xyz``      |
+--------------------------------+----------------+-----------------------------------------------------------+
| max_force                      | 1.0e-3         | Force convergence criterion for geometry relaxations      |
+--------------------------------+----------------+-----------------------------------------------------------+
| frozen_types                   | []             | List of element symbols that should not be moved/swapped  |
+--------------------------------+----------------+-----------------------------------------------------------+
| sanity_check_freq              | 1000           | Steps between sanity-check energy recomputations          |
+--------------------------------+----------------+-----------------------------------------------------------+
| print_freq                     | 1              | Steps between printing energies to ``mc.log``             |
+--------------------------------+----------------+-----------------------------------------------------------+
| dump_archive                   | false          | Write periodic snapshots to ``archive.xyz``               |
+--------------------------------+----------------+-----------------------------------------------------------+
| archive_frequency              | 1000           | Interval (steps) between archived snapshots               |
+--------------------------------+----------------+-----------------------------------------------------------+
| save_downhill                  | false          | Save every lower-energy structure to ``downhill.xyz``     |
+--------------------------------+----------------+-----------------------------------------------------------+
| equil_steps                    | 0              | Steps before accumulating averages and statistics         |
+--------------------------------+----------------+-----------------------------------------------------------+
| writestats                     | false          | Write machine-readable energy/cell CSV to ``stats``       |
+--------------------------------+----------------+-----------------------------------------------------------+
| writestats_freq                | 10             | Interval (steps) between ``stats`` file writes            |
+--------------------------------+----------------+-----------------------------------------------------------+
| **potential**                                                                                               |
+--------------------------------+----------------+-----------------------------------------------------------+
| arch                           | mace_mp        | ASE calculator / MLIP architecture                        |
+--------------------------------+----------------+-----------------------------------------------------------+
| model                          | *required*     | Path to the MLIP model file                               |
+--------------------------------+----------------+-----------------------------------------------------------+
| device                         | cpu            | Device for model evaluation: ``cpu``, ``cuda``, ``mps``   |
+--------------------------------+----------------+-----------------------------------------------------------+
| precision                      | float64        | Numerical precision: ``float64`` or ``float32``           |
+--------------------------------+----------------+-----------------------------------------------------------+
| dispersion                     | false          | Switch on dispersion correction                           |
+--------------------------------+----------------+-----------------------------------------------------------+
| threads_per_rank               | 1              | Threads used per MPI rank for model evaluation            |
+--------------------------------+----------------+-----------------------------------------------------------+
| **monte**                                                                                                   |
+--------------------------------+----------------+-----------------------------------------------------------+
| steps                          | 0              | Total number of Monte Carlo / basin-hopping steps         |
+--------------------------------+----------------+-----------------------------------------------------------+
| max_distance                   | 0.001          | Maximum single-atom displacement (Å)                      |
+--------------------------------+----------------+-----------------------------------------------------------+
| distance_update_freq           | 1000           | Steps between displacement-size updates                   |
+--------------------------------+----------------+-----------------------------------------------------------+
| distance_ratio                 | 0.37           | Target acceptance ratio for atom moves                    |
+--------------------------------+----------------+-----------------------------------------------------------+
| vol_move_freq                  | 0              | Frequency weight of volume moves (MC only)                |
+--------------------------------+----------------+-----------------------------------------------------------+
| vol_move_symmetry              | cubic          | Volume-move symmetry: cubic, tetragonal, orthorhombic,    |
|                                |                | vectors                                                   |
+--------------------------------+----------------+-----------------------------------------------------------+
| max_vol_displacement           | 0.1            | Maximum fractional volume displacement                    |
+--------------------------------+----------------+-----------------------------------------------------------+
| vol_update_freq                | 1000           | Steps between volume-displacement updates                 |
+--------------------------------+----------------+-----------------------------------------------------------+
| vol_ratio                      | 0.37           | Target acceptance ratio for volume moves                  |
+--------------------------------+----------------+-----------------------------------------------------------+
| atom_move_freq                 | 0              | Frequency weight of single-atom translation moves (MC)    |
+--------------------------------+----------------+-----------------------------------------------------------+
| move_types                     | []             | Element symbols allowed to undergo translation moves      |
+--------------------------------+----------------+-----------------------------------------------------------+
| transmutate_frequency          | 0              | Frequency weight of transmutation moves                   |
+--------------------------------+----------------+-----------------------------------------------------------+
| mutate_types_1                 | []             | First element in each transmutation pair                  |
+--------------------------------+----------------+-----------------------------------------------------------+
| mutate_types_2                 | []             | Second element in each transmutation pair                 |
+--------------------------------+----------------+-----------------------------------------------------------+
| transmute_chem_pot             | []             | Chemical potential terms for semi-grand transmutation     |
+--------------------------------+----------------+-----------------------------------------------------------+
| **moves.swaps**                                                                                             |
+--------------------------------+----------------+-----------------------------------------------------------+
| pairs                          | []             | List of [A, B] pairs whose positions are swapped          |
+--------------------------------+----------------+-----------------------------------------------------------+
| counterion                     | None           | Optional third atom relocated during a combined swap      |
+--------------------------------+----------------+-----------------------------------------------------------+
| max_distance                   | 3.0            | Search radius (Å) for relocating the counterion           |
+--------------------------------+----------------+-----------------------------------------------------------+
| frequency                      | 100            | Frequency weight of this swap move                        |
+--------------------------------+----------------+-----------------------------------------------------------+
| **md**                                                                                                      |
+--------------------------------+----------------+-----------------------------------------------------------+
| move_freq                      | 0              | Frequency weight of MD shake moves (basin hopping)        |
+--------------------------------+----------------+-----------------------------------------------------------+
| timestep_fs                    | 2.0            | MD time step (femtoseconds)                               |
+--------------------------------+----------------+-----------------------------------------------------------+
| temperature_k                  | 1000.0         | Target temperature for the MD thermostat (K)              |
+--------------------------------+----------------+-----------------------------------------------------------+
| friction_fs                    | 0.01           | Langevin friction coefficient (1/fs)                      |
+--------------------------------+----------------+-----------------------------------------------------------+
| steps                          | 1000           | Number of MD steps in each shake move                     |
+--------------------------------+----------------+-----------------------------------------------------------+
| **relax**                                                                                                   |
+--------------------------------+----------------+-----------------------------------------------------------+
| method                         | lbfgs          | Optimiser: ``lbfgs``, ``bfgs``, or ``fire``               |
+--------------------------------+----------------+-----------------------------------------------------------+
| steps                          | 1000           | Maximum relaxation steps                                  |
+--------------------------------+----------------+-----------------------------------------------------------+
| tol                            | 1.0e-3         | Force/energy convergence tolerance                        |
+--------------------------------+----------------+-----------------------------------------------------------+
| style                          | conp           | Cell degrees of freedom: ``conp``, ``cona``, ``conv``     |
+--------------------------------+----------------+-----------------------------------------------------------+
| **grid**                                                                                                    |
+--------------------------------+----------------+-----------------------------------------------------------+
| x, y, z                        | 10, 10, 10     | Grid resolution for combined-swap counterion placement    |
+--------------------------------+----------------+-----------------------------------------------------------+
| cut                            | 2.0            | Exclusion distance (Å) around existing atoms on the grid  |
+--------------------------------+----------------+-----------------------------------------------------------+


=============
Output Files
=============

The simulation creates the following files. In multi-rank (MPI) runs each filename is suffixed with ``_rank<N>`` before the extension (e.g. ``mc_rank0.log``), so multiple replicas can run in the same directory without overwriting each other.

* **mc.log**: Comprehensive log for the rank, detailing calculation steps, accepted moves, energies, and fluctuations.
* **archive.xyz**: History file of periodically dumped structural images.
* **restart.xyz**: Final structural state of the rank, which can be copied directly to ``basin.xyz`` to continue a run.
* **accepted.xyz**: Captured structural coordinates of every accepted state (basin hopping / AIRSS).
* **downhill.xyz**: Every configuration that results in an energy drop (basin hopping / AIRSS).
* **stats**: Machine-readable CSV of step, total energy, and cell parameters (only written when ``writestats`` is true).
* **<log>_relax.log**: Optional per-rank auxiliary log capturing geometry-relaxation output from the MLIP calculator.

======================================
Parallel Replica Execution (MPI)
======================================

PyMC-NMR features built-in, embarrassingly parallel support for launching multiple independent simulation replica trajectories simultaneously using MPI.

To execute 12 parallel copies of your simulation across a cluster or locally:

.. code-block:: bash

   mpirun -np 12 python3 -m pymc_nmr.monte --config config.yaml --structure basin.xyz --log mc.log

If ``mpi4py`` is not installed, the program runs as a single rank and emits a warning when an MPI runtime is detected.

------------------------
Differentiated Outputs
------------------------

When running in parallel mode (where total processors > 1):

* **No Write Conflicts:** Ranks automatically suffix their output filenames with ``_rank<N>`` before the extension (e.g. ``mc_rank0.log``, ``archive_rank0.xyz``, ``accepted_rank0.xyz``) to avoid write corruption.
* **Unique Trajectories:** Seeds for random number generators are dynamically varied by rank (``seed = 42 + rank``) to guarantee that each replica explores a different trajectory of the configuration landscape.
* **Silenced Console Output:** Console outputs are automatically silenced for all ranks except rank 0, keeping your stdout logs clean.

------------------------
Recombining Results
------------------------

Once execution finishes, combine the independent runs and perform unified Boltzmann ensemble weighting using the post-processing script:

.. code-block:: bash

   # 1. Combine all accepted coordinates into a single dataset
   cat accepted_rank*.xyz > combined_accepted.xyz

   # 2. Run post-processing weights calculator
   python3 -m pymc_nmr.weights --input combined_accepted.xyz --output final_weights

The output is a table with one row per retained structure:

.. code-block:: text

   # frame  energy_eV  relative_energy_eV  boltzmann_factor  weight

Rows whose weight falls below 1e-6 are omitted, so the printed weights do not
sum exactly to one; the total over all frames is reported on stdout.

.. _weighting-caveat:

--------------------------------------------
Known limitation: how to weight the ensemble
--------------------------------------------

.. warning::

   **The correct weighting of** ``accepted.xyz`` **is a scientific decision that
   PyMC-NMR does not make for you. Read this before publishing ensemble
   averages.**

``accepted.xyz`` is appended to on every Metropolis-accepted move. It is
therefore a Markov chain already distributed according to
:math:`\exp(-E / k_{B}T_{\mathrm{run}})` at the temperature of the run.

The ``weights`` script applies a Boltzmann factor to those same structures. If
the weighting temperature equals the run temperature, the energy dependence is
counted twice, giving an effective :math:`\exp(-2E / k_{B}T)`.

Two interpretations are self-consistent, and you must decide which one your
analysis assumes:

1. **Chain semantics.** Treat ``accepted.xyz`` as an MC trajectory. Ensemble
   averages then use uniform weights :math:`1/N`, and the ``weights`` script
   should not be used.
2. **Pool semantics.** Treat ``accepted.xyz`` as a pool of candidate orderings
   discovered by the search. Boltzmann reweighting is then appropriate, but the
   pool must first be reduced to symmetry-distinct structures, and each must be
   weighted by its configurational multiplicity. The MC temperature is then only
   a search parameter, not the physical temperature of the ensemble.

The ``weights`` script currently performs **neither** the deduplication nor the
multiplicity weighting required by interpretation 2. For site-ordering problems
such as Si/Al distribution, the omitted multiplicity is a leading-order effect,
not a small correction.

--------------------------------------------------
Known limitation: unconverged geometry relaxations
--------------------------------------------------

.. warning::

   In ``basinhop`` and ``airss`` modes, if a geometry relaxation does not reach
   ``relax.tol`` within ``relax.steps``, the energy is replaced by a sentinel
   value of ``1.0e6`` eV rather than raising an error.

That sentinel is a finite number, so it propagates through the Metropolis test
and into the reported statistics like any other energy. The consequences are:

* If the **new** state is unconverged but the old one is fine, the enormous
  energy difference causes the move to be rejected. This is the intended
  behaviour.
* If **both** states are unconverged, their difference is exactly zero, so
  :math:`\exp(-\beta \Delta E) = 1` and the move is accepted unconditionally.
  Because the current energy is carried forward, a single early failure can
  poison the rest of the run.

The symptom is a swap acceptance ratio close to 100% together with a reported
average energy of ``1.0e+06``. PyMC-NMR now writes an explicit warning to
``mc.log`` every time this happens — **grep your log for**
``WARNING: relaxation failed`` before trusting a run:

.. code-block:: bash

   grep -c "WARNING: relaxation failed" mc.log

If it appears, increase ``relax.steps``, loosen ``relax.tol``, or check that
the starting structure is sensible.
