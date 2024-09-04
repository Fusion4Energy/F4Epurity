import numpy as np
import os
import pytest
import pandas as pd

from f4epurity.dose import (
    dose_from_line_source,
    convert_to_dose,
    extract_dose_factors,
    write_vtk_file,
)


@pytest.fixture
def dose_df():
    # Create a dummy DataFrame
    data = {"Nuclide": ["X1", "X2", "X3"], "Dose Factor": [0.5, 0.7, 0.9]}
    df = pd.DataFrame(data)
    return df


# Test test_dose_from_line_source
def test_dose_from_line_source():
    # Define the input parameters for the test case
    dose = [1, 2, 3, 4, 5]
    x1, y1, z1 = 0, 0, 0
    x2, y2, z2 = 0, 0, 10
    x, y, z = 10, 0, 8

    # I have calculated analytically the expected result for the the given coordinates
    # Dose = Dose * theta / (4 * pi * x)
    expected_result = 1.5 * 0.87266 / (4 * np.pi * np.sqrt(100))

    result = dose_from_line_source(dose, x1, y1, z1, x2, y2, z2, x, y, z)

    # Set a 1% tolerance for the test
    assert np.isclose(result, expected_result, rtol=0.01)


def test_extract_dose_factors(dose_df):
    # Test with a nuclide that is in the DataFrame
    nuclide = "X1"
    result = extract_dose_factors(nuclide, dose_df)
    assert result == 0.5

    # Test with a nuclide that is not in the DataFrame
    nuclide = "X4"
    result = extract_dose_factors(nuclide, dose_df)
    assert result == f"Dose conversion factor not found for nuclide {nuclide}."


def test_convert_to_dose(dose_df):
    # Define the input parameters for the test case
    nuclide = "X1"
    activity = [1, 2, 3, 4, 5]

    # Calculate the expected result for the given inputs
    expected_result = [0.5 * act * 1e6 for act in activity]

    result = convert_to_dose(nuclide, activity, dose_df)

    assert np.allclose(result, expected_result)
