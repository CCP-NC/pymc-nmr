# pymc-nmr test suite

Run the whole suite from anywhere in the repository:

```bash
pytest
```

## What needs the MACE model

Most tests are **physics regression tests**: they assert an energy, or an energy
difference, against a value produced by a specific potential. They therefore
need `MACE-matpes-r2scan-omat-ft.model` (obtainable from the mace-mp project)
placed in `tests/data/`:

```
tests/data/MACE-matpes-r2scan-omat-ft.model
```

That file is ~76 MB and is not stored in the repository. When it is absent
those tests **skip** with a message telling you this, rather than failing, so
`pytest` is green on a fresh checkout. Use `pytest -rs` to list what skipped.

Tests run on the CPU. Shared fixtures live in `conftest.py`, which also pins the
NumPy seed so that "the atom the RNG happens to pick" is reproducible.

## The tests

| Test | What it checks |
| --- | --- |
| `test_1` | ASE, janus-core and the MACE model load, and give the reference single-point energy for `basin.xyz`. |
| `test_2` | An ASE LBFGS cell+geometry relaxation converges to the reference relaxed energy. |
| `test_3` | A single-point energy computed *through pymc-nmr* (`Field`, `Config`, `Species`) equals the one from plain ASE. |
| `test_4` | As `test_3`, for a relaxed energy. |
| `test_5` | Energy change of a Si/Al swap, and that swapping back restores the original energy exactly. |
| `test_6` | As `test_5`, with an LBFGS relaxation after each swap (the basin-hopping path). |
| `test_7` | Energy change of a single-atom translation, and that undoing it restores the original energy. |
| `test_8` | Isotropic (cubic) volume move, and that restoring the cell restores the energy. |
| `test_9` | As `test_8`, for an anisotropic cell distortion. |
| `test_11` | The DFT-D3 dispersion correction agrees between janus-core and mace-torch. Skips unless `torch-dftd` is installed. |
| `test_write_statistics` | The CSV statistics output has the right header and rows. Needs no model. |

In tests 5–9 the "undo" assertion is the scientifically important half: a move
that is not exactly reversible breaks detailed balance in the Metropolis loop.
The forward energy change is a regression value tied to the model above.
