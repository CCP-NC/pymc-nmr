"""Dispersion correction agrees between janus-core and mace-torch directly.

Both routes should produce the same DFT-D3 corrected energy; a difference means
one of the two is applying the correction differently.
"""

import pytest
from ase.io import read

from conftest import ARCH, DEVICE, PRECISION, requires_model


@requires_model
def test_dispersion_consistent_between_janus_and_mace(model_path):
    # janus-core routes the D3 correction through torch-dftd (the "d3" extra),
    # which is what actually has to be installed here.
    pytest.importorskip(
        "torch_dftd",
        reason="torch-dftd (janus-core's 'd3' extra) is required for the dispersion test",
    )

    from janus_core.helpers.mlip_calculators import choose_calculator
    from mace.calculators import mace_mp

    atoms = read("basin.xyz")

    atoms.calc = choose_calculator(
        arch=ARCH, dispersion=True, model=model_path, precision=PRECISION, device=DEVICE
    )
    janus_energy = atoms.get_potential_energy()

    atoms.calc = mace_mp(
        model=model_path, dispersion=True, default_dtype=PRECISION, device=DEVICE
    )
    mace_energy = atoms.get_potential_energy()

    assert mace_energy - janus_energy == pytest.approx(0.0)
