# CHA example — chabazite, Si:Al = 5, with Brønsted protons

A worked example of searching Al/Si orderings in the zeolite chabazite, with a
charge-compensating proton relocated alongside each Al.

## Order to work through

1. **`gen_initial_CHA.ipynb`** — builds the starting structures from `CHA.cif`:
   substitutes Al to a Si:Al ratio of 5, enforces Löwenstein's rule (no
   Al–O–Al), then adds one Brønsted H per Al. Writes
   `CHA_SAR5_with_H.extxyz` plus 20 numbered variants.
   *Optional* — `CHA_SAR5_with_H.extxyz` is already provided.
2. **`tutorial1.ipynb`** — the place to start. Creates `basin.xyz`, writes a
   short `config_tutorial.yaml`, and runs a basin-hopping search.
3. **`tutorial2.ipynb`** — compares the `monte`, `basinhop` and `airss` modes.

## Files

| File | Tracked | Notes |
| --- | --- | --- |
| `CHA.cif` | yes | Pristine all-silica chabazite, 36 T-sites. |
| `CHA_SAR5_with_H.cif` / `.extxyz` | yes | Al-substituted, protonated starting structure. Input to `tutorial1.ipynb`. |
| `config_production.yaml` | yes | Realistic long run (200 steps). **Needs a GPU** — see below. |
| `CHA_SAR5_with_H_0NN.extxyz` | no | 20 further candidates from the generator notebook. Regenerate as needed. |
| `basin.xyz`, `config_tutorial.yaml`, `config_{mc,bh,airss}.yaml` | no | Written by the notebooks. |
| `MACE-matpes-r2scan-omat-ft.model` | no | ~76 MB. Download from the mace-mp project and place it here. |

## Two things that surprise people

**The protons look unbonded in the starting structure.** In
`CHA_SAR5_with_H.extxyz` each H sits ~2.5 Å from the nearest oxygen, not the
~0.98 Å of a real O–H bond, because the generator places them on a sphere
around each Al without attaching them to anything. This is an intentionally
crude guess. The first geometry relaxation seats all six protons on Al–O–Si
bridges at 0.97–1.00 Å — genuine Brønsted acid sites — and lowers the energy by
about 23 eV. Nothing is wrong.

**You really want a GPU.** Each basin-hopping step runs a full geometry
relaxation, so one move costs hundreds of potential evaluations. Set
`device: cuda` in the config (`mps` on Apple Silicon) and check PyTorch can see
it first:

```bash
python -c "import torch; print(torch.cuda.is_available())"
```

On CPU, `config_production.yaml` is an overnight job; the tutorial config is
deliberately tiny so it finishes quickly.
