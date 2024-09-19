import os
import numpy as np
import pandas as pd
import pytest

from f4epurity.utilities import (
    get_molar_mass,
    get_isotopes,
    add_user_irrad_scenario,
    normalise_nuclide_name,
    convert_names,
    sum_vtr_files,
)


@pytest.fixture
def nist_df():
    # Create a sample DataFrame for testing
    data = {
        "Atomic Symbol": ["co", "co", "ta", "ta"],
        "Mass Number": [59, 60, 180, 181],
        "Standard Atomic Weight": [58.933195, 58.933195, 180.94788, 180.94788],
        "Isotopic Composition": [1.0, np.nan, 0.00012, 0.99988],
    }
    return pd.DataFrame(data)


def test_get_molar_mass(nist_df):
    # Test case 1: Valid isotope
    isotope = "co59"
    expected_mass = 58.933195
    assert get_molar_mass(isotope, nist_df) == expected_mass

    # Test case 2: Invalid isotope
    isotope = "nb94"
    with pytest.raises(SystemExit):
        get_molar_mass(isotope, nist_df)

    # Test case 3: Valid isotope with different mass number
    isotope = "ta181"
    expected_mass = 180.94788
    assert get_molar_mass(isotope, nist_df) == expected_mass


def test_get_isotopes(nist_df):
    # Test case 1: Valid element with one naturally occurring isotope
    element = "Co"
    expected_isotopes = ["co59"]
    assert get_isotopes(element, nist_df) == expected_isotopes

    # Test case 2: Valid element with 2 naturally occuring isotopes
    element = "Ta"
    expected_isotopes = ["ta180", "ta181"]
    assert get_isotopes(element, nist_df) == expected_isotopes

    # Test case 3: Invalid element
    element = "X"
    expected_isotopes = []
    assert get_isotopes(element, nist_df) == expected_isotopes


def test_add_user_irrad_scenario():
    # Test case 1: Valid file and data
    day_to_sec = 24 * 60 * 60
    my_scenario = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data", "my_irrad_scenario.txt"
    )
    irrad_scenarios = {}
    add_user_irrad_scenario(my_scenario, irrad_scenarios)
    expected_scenario = {
        "times": [0 * day_to_sec, 86400 * day_to_sec, 172800 * day_to_sec],
        "fluxes": [1.0, 2.0, 3.0],
    }
    assert irrad_scenarios[my_scenario] == expected_scenario

    # Test case 2: Invalid file
    filename = "invalid_file.txt"
    irrad_scenarios = {}
    with pytest.raises(FileNotFoundError):
        add_user_irrad_scenario(filename, irrad_scenarios)
    assert filename not in irrad_scenarios


def test_normalise_nuclide_name():
    # Test case 1: Nuclide with no suffix
    nuclide = "co59"
    expected_result = "Co59"
    assert normalise_nuclide_name(nuclide) == expected_result

    # Test case 2: Nuclide with 'm' suffix
    nuclide = "co59m"
    expected_result = "Co59m"
    assert normalise_nuclide_name(nuclide) == expected_result

    # Test case 3: Nuclide with 'n' suffix
    nuclide = "co59n"
    expected_result = "Co59n"
    assert normalise_nuclide_name(nuclide) == expected_result

    # Test case 4: Nuclide with leading zeros in the mass number
    nuclide = "co059"
    expected_result = "Co59"
    assert normalise_nuclide_name(nuclide) == expected_result

    # Test case 5: Nuclide with lowercase element name
    nuclide = "co59"
    expected_result = "Co59"
    assert normalise_nuclide_name(nuclide) == expected_result

    # Test case 6: Nuclide with mixed case element name
    nuclide = "Co59"
    expected_result = "Co59"
    assert normalise_nuclide_name(nuclide) == expected_result


def test_convert_names():
    # Test case 1: Empty dictionary
    nuc_dict = {}
    expected_result = {}
    assert convert_names(nuc_dict) == expected_result

    # Test case 2: Dictionary with one parent and one daughter
    nuc_dict = {"co59": {"atoms": 10, "reactions": {"Co60m": 5}}}
    expected_result = {"Co059": {"atoms": 10, "reactions": {"Co060m": 5}}}
    assert convert_names(nuc_dict) == expected_result

    # Test case 3: Dictionary with multiple parents and daughters
    nuc_dict = {
        "ta180": {"atoms": 10, "reactions": {"Ta181": 5, "Ta181m": 3}},
        "TA181": {"atoms": 20, "reactions": {"ta182": 5, "ta182m": 3}},
    }
    expected_result = {
        "Ta180": {"atoms": 10, "reactions": {"Ta181": 5, "Ta181m": 3}},
        "Ta181": {"atoms": 20, "reactions": {"Ta182": 5, "Ta182m": 3}},
    }
    assert convert_names(nuc_dict) == expected_result


@pytest.mark.parametrize(
    ["masses", "expected"],
    [
        [None, np.ones((2, 2, 2)) * 2],  # simple sum
        [[1, 4], np.ones((2, 2, 2)) * 5],  # weighted sum
    ],
)
def test_sum_vtr_files(masses, expected, tmpdir):
    steps = 3
    x = np.linspace(0, 1, num=steps)
    y = np.linspace(0, 1, num=steps)
    z = np.linspace(0, 1, num=steps)
    dose_array = np.ones((len(x) - 1, len(y) - 1, len(z) - 1))
    dose_arrays = [dose_array] * 2

    assert (
        sum_vtr_files(dose_arrays, x, y, z, tmpdir, masses=masses) == expected
    ).all()
