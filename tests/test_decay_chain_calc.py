import numpy as np
import os
import json
import pytest

from f4epurity.decay_chain_calc import calculate_total_activity
from f4epurity.utilities import normalise_nuclide_name


def test_calculate_total_activity():
    decay_data_path = os.path.join(os.path.dirname(__file__), "data", "Decay2020.json")
    with open(decay_data_path, "r") as f:
        decay_data = json.load(f)
    # Define dummy nuclide dictionary
    nuclide_dict = {
        "co59": {
            "atoms": 1.02186e25,
            "reactions": {
                "co60m": [1.08e-14],
                "co60": [3.05e-14],
            },
        },
        "Nb093": {
            "atoms": 6.48194e24,
            "reactions": {
                "Nb094m": [2.49e-15],
                "Nb094": [1.16e-15],
                "Nb092m": [4.02e-18],
                "Nb092": [6.72e-18],
                "Nb093m": [3.93e-17],
            },
        },
    }
    irrad_scenario = "DT1"
    decay_time = 1e6  # 12 days

    # Define the expected output from FISACT-II calculation
    expected_output = {
        "Co060": np.array([2.53773e08]),
        "Co060m": np.array([0.0]),
        "Nb094m": np.array([0.0]),
        "Nb094": np.array([5.05530e03]),
        "Nb092m": np.array([2.14198e04]),
        "Nb092": np.array([5.32794e-03]),
        "Nb093m": np.array([6.07715e04]),
    }

    # Convert the numpy arrays to floats
    output = calculate_total_activity(
        nuclide_dict, irrad_scenario, decay_time, decay_data
    )
    output = {
        key: [float(value[0]) for value in values] for key, values in output.items()
    }

    assert output == pytest.approx(expected_output, rel=0.05, abs=1e-5)


# Test the time correction factors are the same as in D1SUNED
def test_correction_factors():

    decay_data_path = os.path.join(os.path.dirname(__file__), "data", "Decay2020.json")
    with open(decay_data_path, "r") as f:
        decay_data = json.load(f)

    irrad_sch = "SA2"
    decay_time = 1e6

    time_factors_path = os.path.join(
        os.path.dirname(__file__), "data", "d1s_time_correction_SA2"
    )
    with open(time_factors_path, "r") as f:
        time_correction_factors = {}

        for line in f:
            stripped_line = line.strip()
            if (
                stripped_line.startswith("C")
                or stripped_line.startswith("nsc")
                or not stripped_line
            ):
                continue

            words = line.split()

            time_fact_0 = float(words[2])

            nuclide = words[3]

            time_correction_factors[nuclide] = time_fact_0

    for nuclide in time_correction_factors:
        nuclide_dict = {"co59": {"atoms": 1e25, "reactions": {nuclide: [1e-14]}}}
        activity = calculate_total_activity(
            nuclide_dict, irrad_sch, decay_time, decay_data
        )

        updated_activity = {}
        for nuclide, nuclide_activity in activity.items():
            nuclide = normalise_nuclide_name(nuclide)
            updated_activity[nuclide] = nuclide_activity

    # Set a relative tolerance of 5%
    assert updated_activity[nuclide][0] / 1e25 / 1e-14 == pytest.approx(
        time_correction_factors[nuclide], rel=0.05
    )
