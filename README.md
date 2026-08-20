# pymc-nmr

Generate ensembles of candidate crystal structures — for example the possible
Si/Al orderings in a zeolite — and their energies, so that computed NMR
parameters can be averaged over a realistic set of configurations.

Energies come from a machine-learned interatomic potential (MACE, via
[janus-core](https://github.com/stfc/janus-core)), so a search that would be
prohibitive with DFT becomes routine.

## Install

Requires Python 3.10+. Using [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/CCP-NC/pymc-nmr.git
cd pymc-nmr
uv venv
source .venv/bin/activate
uv pip install -e .
```

## Run

You need a structure file, a MACE model file, and a `config.yaml`:

```yaml
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
```

Then:

```bash
pymc-nmr --config config.yaml --structure basin.xyz
```

Progress goes to `mc.log`; accepted structures accumulate in `accepted.xyz`.

### Use a GPU

Runtime is dominated by potential evaluations — each basin-hopping step runs a
full geometry relaxation — so **a GPU is typically 1–2 orders of magnitude
faster than a CPU**. Set it in the config:

```yaml
potential:
  device: cuda    # NVIDIA GPU; "mps" on Apple Silicon; "cpu" otherwise
```

Check PyTorch can see it first, or `device: cuda` will fail:

```bash
python -c "import torch; print(torch.cuda.is_available())"
```

CPU runs are fine for tutorials and for validating a config, but too slow for
production.

New users should follow the
[Quickstart](docs/source/quickstart.rst), which explains each step, or the
worked notebooks in `example/CHA/`.

## Modes of operation

All three are driven by **swap moves**, in which two atoms chosen at random
exchange positions. The move is accepted with the Metropolis probability
*acc(o → n) = min(1, exp(−β[E<sub>n</sub> − E<sub>o</sub>]))*.

| Mode | What it does | When to use it |
| --- | --- | --- |
| `monte` | Standard Monte Carlo in the canonical (*NVT*) or isobaric-isothermal (*NPT*) ensemble, with particle translation and volume moves as well as swaps. | Sampling at a well-defined thermodynamic state point. |
| `basinhop` | Each swap is followed by a geometry relaxation, and the Metropolis test is applied to the relaxed energies. Optionally a short Langevin MD "shake" helps escape a basin. | Site-ordering problems, where swapping two species costs a lot of energy before relaxation. The usual choice. |
| `airss` | All atoms of the chosen types are randomised, then relaxed — *ab initio* random structure search. | Exploring broadly when no good starting ordering is known. |

## Units

ASE conventions throughout: energies in eV, lengths in Å, temperature in K, and
therefore pressure in eV/Å³ (**not** GPa).

## Weighting the results — read this first

Converting `accepted.xyz` into ensemble-averaged NMR parameters requires a
decision that this code does not make for you: whether that file is a Markov
chain (weight uniformly) or a pool of candidate orderings (deduplicate to
symmetry-distinct structures and weight by Boltzmann factor **and**
multiplicity). Applying `pymc_nmr.weights` naively to a Metropolis-sampled
trajectory double-counts the energy. See "Known limitation" in the
[User Guide](docs/source/UserGuide.rst) before publishing averages.

## Tests

```bash
uv pip install -e ".[test]"
pytest
```

Most tests are physics regression tests that need a MACE model in `tests/data/`;
they skip automatically when it is absent. See [tests/README.md](tests/README.md).

## Documentation

```bash
uv pip install -e ".[docs]"
sphinx-build -b html docs/source docs/build/html
```

## Licence

MIT. See [LICENSE](LICENSE).
