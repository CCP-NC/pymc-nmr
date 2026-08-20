"""Energy change of a Si/Al swap move, and its reversibility.

The reversibility check (``delta_back == 0``) is the important one: undoing a
swap must restore the original energy exactly, or detailed balance in the
Metropolis loop is broken.
"""

import pytest

from conftest import requires_model


@requires_model
def test_swap_energy_and_reversibility(field, basin):
    old_energy = field.calculate_energy(basin.create_atoms_object(), False)

    atm1 = basin.select_atom_of_type("Si")
    atm2 = basin.select_atom_of_type("Al")

    basin.swap_atom_positions(atm1, atm2)
    swapped_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_swap = old_energy - swapped_energy

    basin.swap_atom_positions(atm1, atm2)
    restored_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_back = old_energy - restored_energy

    assert delta_swap == pytest.approx(0.17007792199865435)
    assert delta_back == pytest.approx(0.0)
