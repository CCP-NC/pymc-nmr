"""As test_5, but the swap is followed by an LBFGS relaxation (basin hopping).

Reversibility is only approximate here: relaxing, swapping back and relaxing
again need not land on bit-identical coordinates, so an absolute tolerance is
used.
"""

import pytest

from conftest import requires_model


@requires_model
def test_swap_relax_energy_and_reversibility(field, basin):
    def relaxed():
        return field.calculate_energy_relax(
            basin.create_atoms_object(), "lbfgs", 1000, 1.0e-3, "conp", False
        )

    old_energy = relaxed()

    atm1 = basin.select_atom_of_type("Si")
    atm2 = basin.select_atom_of_type("Al")

    basin.swap_atom_positions(atm1, atm2)
    delta_swap = old_energy - relaxed()

    basin.swap_atom_positions(atm1, atm2)
    delta_back = old_energy - relaxed()

    assert delta_swap == pytest.approx(-0.12121373476838926)
    assert delta_back == pytest.approx(0.0, rel=1.0e-6, abs=1.0e-6)
