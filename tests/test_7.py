"""Energy change of a single-atom translation move, and its reversibility."""

import pytest

from conftest import requires_model


@requires_model
def test_atom_move_energy_and_reversibility(field, basin):
    old_energy = field.calculate_energy(basin.create_atoms_object(), False)

    atm = basin.select_atom_of_type("Si")
    old_pos = basin.make_atom_move(atm, 0.1)

    moved_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_move = old_energy - moved_energy

    basin.reject_atom_move(atm, old_pos)
    restored_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_back = old_energy - restored_energy

    assert delta_move == pytest.approx(0.0268727707330072)
    assert delta_back == pytest.approx(0.0)
