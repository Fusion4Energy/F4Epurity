import pytest
import numpy as np

from f4epurity.reaction_rate import calculate_reaction_rate


def test_calculate_reaction_rate():
    # Line source case
    delta_impurity = 10
    sigma_eff = np.array([1, 2, 3])
    flux_spectrum = np.array([4, 5, 6])
    expected_output = np.array([4e-23, 1e-22, 1.8e-22])
    output = calculate_reaction_rate(delta_impurity, sigma_eff, flux_spectrum)
    assert np.allclose(output, expected_output, rtol=1e-6)

    # Point source case
    delta_impurity = 10
    sigma_eff = 2
    flux_spectrum = np.array([4, 5, 6])
    expected_output = 1.5e-22
    output = calculate_reaction_rate(delta_impurity, sigma_eff, flux_spectrum)
    assert np.isclose(output, expected_output, rtol=1e-6)
