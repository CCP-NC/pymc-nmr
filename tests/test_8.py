"""Isotropic (cubic) volume move: energy change and restoration of the cell."""

import numpy as np
import pytest

from conftest import requires_model


@requires_model
def test_cubic_volume_move_and_reversibility(field, basin):
    old_energy = field.calculate_energy(basin.create_atoms_object(), False)

    old_vec = basin.get_vectors()
    old_pos = basin.get_positions()

    bulks = np.ones(3, dtype=np.float64)
    basin.expand_cell_cubic(bulks, 0.1)

    expanded_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_expand = old_energy - expanded_energy

    basin.set_positions(old_pos)
    basin.set_vectors(old_vec)
    restored_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_back = old_energy - restored_energy

    assert delta_expand == pytest.approx(0.558811339897602)
    assert delta_back == pytest.approx(0.0)
