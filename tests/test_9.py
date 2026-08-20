"""Anisotropic cell distortion: energy change and restoration of the cell."""

import numpy as np
import pytest

from conftest import requires_model


@requires_model
def test_cell_distortion_and_reversibility(field, basin):
    old_energy = field.calculate_energy(basin.create_atoms_object(), False)

    old_vec = basin.get_vectors()
    old_pos = basin.get_positions()

    bulks = np.zeros(6, dtype=np.float64)
    max_vol = np.full(6, 0.1, dtype=np.float64)
    indx = int(6.0 * np.random.random())
    basin.distort_cell(indx, bulks, max_vol)

    distorted_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_distort = old_energy - distorted_energy

    basin.set_positions(old_pos)
    basin.set_vectors(old_vec)
    restored_energy = field.calculate_energy(basin.create_atoms_object(), False)
    delta_back = old_energy - restored_energy

    assert delta_distort == pytest.approx(-0.5764823374752268)
    assert delta_back == pytest.approx(0.0)
