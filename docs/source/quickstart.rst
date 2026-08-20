==========
Quickstart
==========

This page takes you from a fresh machine to a finished simulation. It assumes
no prior experience with Python packaging or YAML configuration files. If a
step fails, see :ref:`troubleshooting` at the bottom.

What you will end up with
=========================

PyMC-NMR searches for low-energy arrangements of atoms in a crystal — for
example, which of the silicon sites in a zeolite are replaced by aluminium.
It produces a set of candidate structures and their energies, which you can
then use to weight computed NMR parameters.

You need three things:

1. **A structure file** describing your starting crystal (an ``.xyz`` file).
2. **A potential file** — the machine-learned model that computes energies.
3. **A configuration file** (``config.yaml``) saying what to do.


Step 1: Install
===============

PyMC-NMR needs Python 3.10 or newer. We recommend `uv
<https://docs.astral.sh/uv/>`_, which manages both Python and the packages for
you. Install uv first (once per machine):

.. code-block:: bash

   curl -LsSf https://astral.sh/uv/install.sh | sh

Now download PyMC-NMR and install it into its own isolated environment:

.. code-block:: bash

   git clone https://github.com/CCP-NC/pymc-nmr.git
   cd pymc-nmr
   uv venv                  # creates a .venv folder holding a private Python
   source .venv/bin/activate    # on Windows: .venv\Scripts\activate
   uv pip install -e .

.. note::

   The ``source .venv/bin/activate`` line must be repeated **every time you
   open a new terminal**. You can tell it worked because your prompt gains a
   ``(pymc-nmr)`` prefix. Without it, the ``pymc-nmr`` command will not be
   found.

Check the installation:

.. code-block:: bash

   pymc-nmr --help

You should see a list of options. This downloads and installs PyTorch, so the
first install may take several minutes and a few GB of disk.


Step 2: Get a potential
=======================

Energies are computed by a machine-learned interatomic potential (MLIP), not
by PyMC-NMR itself. The examples use a MACE model. Download one from the
`mace-mp project <https://github.com/ACEsuit/mace-mp>`_ and put the
``.model`` file somewhere you can find it — the examples assume it sits in the
same folder you are running from.

These files are large (often 50–100 MB), so they are deliberately not stored
in this repository.


.. _hardware:

Step 2b: Check your hardware — you probably need a GPU
======================================================

Almost all of the runtime is spent evaluating the potential. Every basin-hopping
step performs a geometry relaxation, and every relaxation step is another
evaluation, so a single move costs hundreds of potential calls.

**A GPU is roughly one to two orders of magnitude faster than a CPU for this.**
On a CPU, one basin-hopping step of a ~100-atom cell takes minutes, so a few
hundred steps becomes an overnight job. On a GPU the same run is typically
minutes to hours. For anything beyond a tutorial, use a GPU.

Select the device in the ``potential`` block of your config:

.. code-block:: yaml

   potential:
     model: MACE-matpes-r2scan-omat-ft.model
     device: cuda        # NVIDIA GPU

Valid values are:

``cuda``
   An NVIDIA GPU. This is what you want for production runs.
``mps``
   Apple Silicon (M1/M2/M3...). Faster than ``cpu``, slower than a real GPU,
   and occasionally less numerically robust.
``cpu``
   Works everywhere. Fine for tutorials and for checking that a config is
   valid; too slow for real work.

Before setting ``device: cuda``, confirm that PyTorch can actually see the GPU:

.. code-block:: bash

   python -c "import torch; print(torch.cuda.is_available())"

If this prints ``True`` you are ready. If it prints ``False`` you have a
CPU-only build of PyTorch and ``device: cuda`` will fail; install a CUDA-enabled
build following the `PyTorch instructions <https://pytorch.org/get-started/locally/>`_.

.. tip::

   On a shared HPC machine you usually also need to request a GPU from the
   scheduler (for example ``--gres=gpu:1`` under Slurm) and load the
   appropriate CUDA module. Ask your local support team which is right for
   your cluster.

If you are stuck on CPU, keep runs tractable by lowering ``monte.steps``,
lowering ``relax.steps``, and using the smallest cell that answers your
question. Setting ``potential.threads_per_rank`` to the number of physical
cores you actually have will also help — but do not set it higher, as
oversubscription makes things slower.


Step 3: Write a configuration file
==================================

Create a file called ``config.yaml`` next to your structure file. Here is a
complete, minimal example — every line is explained below:

.. code-block:: yaml

   job:
     structure_method: basinhop
     temperature: 1200.0

   potential:
     model: MACE-matpes-r2scan-omat-ft.model

   monte:
     steps: 50

   moves:
     swaps:
       - pairs:
           - [Al, Si]

What each part means:

``job.structure_method: basinhop``
   Which search algorithm to use. ``basinhop`` swaps two atoms and then relaxes
   the geometry — the usual choice for site-ordering problems. The alternatives
   are ``monte`` (no relaxation) and ``airss`` (randomise everything, then
   relax).

``job.temperature: 1200.0``
   The sampling temperature in Kelvin. This controls how readily the search
   accepts a move that *raises* the energy: higher temperature explores more
   widely, lower temperature stays near the minimum. **There is no default —
   you must state it.**

``potential.model``
   The path to the ``.model`` file from Step 2.

``monte.steps: 50``
   How many moves to attempt. Start small (say 50) to check everything runs,
   then increase.

``moves.swaps.pairs: [[Al, Si]]``
   Swap the positions of an Al atom and an Si atom. This is the move that
   explores different site orderings.

Everything else has a sensible default. Notably, **you do not list the chemical
elements anywhere** — PyMC-NMR reads them from your structure file.

.. note::

   YAML is indentation-sensitive, like Python. Use spaces, never tabs, and keep
   the indentation consistent. If you mistype a key name, PyMC-NMR will stop
   and tell you, rather than silently ignoring it.


Step 4: Run it
==============

.. code-block:: bash

   pymc-nmr --config config.yaml --structure basin.xyz

Both options default to ``config.yaml`` and ``basin.xyz`` in the current
directory, so if your files use those names you can simply run:

.. code-block:: bash

   pymc-nmr

While it runs, progress is written to ``mc.log``. To watch it live, open a
second terminal and run ``tail -f mc.log``.


Step 5: Look at the results
===========================

A completed run leaves several files behind:

``mc.log``
   The human-readable log: energies, which moves were accepted, and running
   averages. **Read this first.**
``accepted.xyz``
   Every structure that was accepted during the search. This is the main
   scientific output.
``downhill.xyz``
   Only those structures that lowered the energy (if ``save_downhill`` is on).
``restart.xyz``
   The final structure. Copy it over ``basin.xyz`` to continue the run.

A healthy run shows the energy in ``mc.log`` falling and then fluctuating
around a plateau, with a swap acceptance rate that is neither ~0% (temperature
too low, or the moves are all unfavourable) nor ~100% (temperature too high to
discriminate).

**Always check for failed relaxations before trusting a basin-hopping run:**

.. code-block:: bash

   grep -c "WARNING: relaxation failed" mc.log

If that count is not zero, some energies in your run are a placeholder value of
1e6 eV rather than a real energy, and the acceptance statistics are
meaningless. Increase ``relax.steps`` or loosen ``relax.tol``. A reported
average energy of ``1.0e+06`` together with a 100% acceptance rate is the
tell-tale sign.

.. warning::

   Before converting ``accepted.xyz`` into ensemble-averaged NMR parameters,
   read :ref:`weighting-caveat` in the User Guide. How these structures should
   be weighted is a genuine scientific decision, and the supplied
   ``weights`` script does not make it for you.


Worked examples
===============

The ``example/`` folder contains complete, runnable cases with Jupyter
notebooks that build the structure, write the config, and run the simulation
step by step:

* ``example/CHA`` — chabazite with Brønsted protons. ``tutorial1.ipynb`` is the
  gentlest starting point; ``tutorial2.ipynb`` compares the three modes.
* ``example/STA30`` — a silicoaluminophosphate zeotype.

To run the notebooks you also need Jupyter:

.. code-block:: bash

   uv pip install -e ".[notebooks]"


.. _troubleshooting:

Troubleshooting
===============

``command not found: pymc-nmr``
   The virtual environment is not active. Run ``source .venv/bin/activate``
   from the project folder.

``Failed to validate config schema``
   There is a mistake in ``config.yaml``. The message names the offending key.
   Common causes: a tab instead of spaces, wrong indentation, a misspelled key,
   or a missing ``temperature``.

``Field required [type=missing] ... temperature``
   ``job.temperature`` has no default and must be set. See Step 3.

``Extra inputs are not permitted``
   A key in your config is not recognised — usually a typo. Check it against
   the keyword table in the User Guide.

The run is extremely slow
   Energy evaluation dominates the cost, and PyMC-NMR is far faster on a GPU.
   See :ref:`hardware`. Failing that, reduce ``monte.steps`` and
   ``relax.steps``, or use a smaller cell.

``CUDA out of memory`` or no GPU present
   Set ``potential.device: cpu`` in your config.
