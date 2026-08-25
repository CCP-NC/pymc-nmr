"""Shared fixtures for the pymc-nmr test suite.

Most tests are *physics regression tests*: they assert an energy or an energy
difference against a value produced by a specific MACE model. That model is a
~76 MB binary and is not stored in the repository, so those tests skip
automatically when it is absent (see ``README``). Tests that need no model
always run.

All tests run with the working directory set to this folder, so the fixture
files ``basin.xyz`` and ``data/`` resolve regardless of where pytest is invoked
from.
"""

import os
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).parent
MODEL_PATH = TEST_DIR / "data" / "MACE-matpes-r2scan-omat-ft.model"
BASIN_XYZ = TEST_DIR / "basin.xyz"

# The reference energies below were generated with this model at float64
# precision. Changing the model, the precision, or the structure invalidates
# them, so they are defined once here rather than repeated in each test.
ARCH = "mace_mp"
PRECISION = "float64"
DEVICE = "cpu"

requires_model = pytest.mark.skipif(
    not MODEL_PATH.exists(),
    reason=(
        f"MACE model not found at {MODEL_PATH}. Download "
        "MACE-matpes-r2scan-omat-ft.model from the mace-mp project and place "
        "it in tests/data/ to run the physics regression tests."
    ),
)


@pytest.fixture(autouse=True)
def _run_in_test_dir(monkeypatch):
    """Run every test from the tests directory so relative paths resolve."""
    monkeypatch.chdir(TEST_DIR)


@pytest.fixture(autouse=True)
def _fixed_seed():
    """Pin the RNG so move-selection is reproducible across runs.

    Several tests assert an exact energy difference for "the swap that numpy
    happens to choose", which is only meaningful with a fixed seed.
    """
    np.random.seed(0)


@pytest.fixture
def model_path():
    return str(MODEL_PATH)


@pytest.fixture
def calculator(model_path):
    """A bare janus-core/ASE calculator, used as the oracle for pymc-nmr."""
    from janus_core.helpers.mlip_calculators import choose_calculator

    return choose_calculator(
        arch=ARCH, model=model_path, precision=PRECISION, device=DEVICE
    )


@pytest.fixture
def field(model_path):
    """A configured pymc-nmr Field wrapping the same potential."""
    from pymc_nmr.field import Field

    fld = Field()
    fld.arch = ARCH
    fld.model = model_path
    fld.device = DEVICE
    fld.precision = PRECISION
    fld.dispersion = False
    fld.species = True
    fld.setup()
    return fld


@pytest.fixture
def basin(tmp_path):
    """The test structure loaded into a pymc-nmr Config object."""
    from pymc_nmr.config import Config
    from pymc_nmr.species import Element, Species

    cfg = Config()
    with open(BASIN_XYZ, "r") as instream:
        cfg.read_config(instream)

    # Species are derived from the structure, mirroring what monte.py does.
    symbols = []
    for s in cfg.symbol:
        if s not in symbols and s != "ghost":
            symbols.append(s)
    spec = Species()
    spec.number_of_elements = len(symbols)
    spec.ele_data = [Element(name=s) for s in symbols]

    with open(os.path.join(tmp_path, "setup.log"), "w") as log:
        cfg.setup_configuration(spec, log)
    return cfg
