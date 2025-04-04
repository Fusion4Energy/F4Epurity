import math
from copy import deepcopy

import numpy as np

from f4epurity.utilities import add_user_irrad_scenario, convert_names

# conversion from days to seconds
DAY_TO_SEC = 24 * 60 * 60
YEAR_TO_SEC = 365.25 * DAY_TO_SEC

# Define the known irradiation scenarios
# (times is a list of length of irradiation periods, fluxes is a list of the relative (to the nominal source strength) source strength for the given irradiation period)
# REF : ITER_D_8WK64Y
IRRAD_SCENARIOS = {
    "DT1": {
        "times": [
            730.5 * DAY_TO_SEC,
            730.5 * DAY_TO_SEC,
            730.5 * DAY_TO_SEC,
            730.5 * DAY_TO_SEC,
            730.5 * DAY_TO_SEC,
            600,
        ],
        "fluxes": [
            1.70427e-06,
            1.17412e-04,
            3.87673e-04,
            9.95946e-04,
            1.58543e-03,
            5.01221e-01,
        ],
    },
    "SA2": {
        "times": [
            2 * YEAR_TO_SEC,
            10 * YEAR_TO_SEC,
            0.667 * YEAR_TO_SEC,
            1.325 * YEAR_TO_SEC,
        ]
        + [3920, 400] * 17
        + [3920, 400] * 3,
        "fluxes": [0.00536, 0.0412, 0, 0.083] + [0, 1] * 17 + [0, 1.4] * 3,
    },
    "Y1": {
        "times": [365 * DAY_TO_SEC],
        "fluxes": [1],
    },
}


def irradiate(
    time: float,
    parent_atoms: int,
    nuclide: str,
    decay_data_dic: dict,
    lambda_temp: list,
    parent: str,
) -> dict:
    """Function to work out the number of atoms for a given nuclide, for a
    given pulse (note this is a recursive function i.e. it calls itself until
    it reaches the end of a decay chain)

    Parameters
    ----------
    time : float
        decay time in seconds
    parent_atoms : int
        number of parent atoms
    nuclide : str
        nuclide name
    decay_data_dic : dict
        dictionary containing decay data
    lambda_temp : list
        list of decay constants
    parent : str
        parent nuclide name

    Returns
    -------
    dict
        dictionary containing the number of atoms for each nuclide in the decay chain
    """

    # Set up a temporary dictionary to record the number of each nuclide in the decay chain
    nuclides = {}
    atoms = 0

    # Calculated atoms for given nuclide
    for i, lambda_i in enumerate(lambda_temp):
        alpha_i = 1
        for j, lambda_j in enumerate(lambda_temp):
            if i != j:
                alpha_i = alpha_i * lambda_j / (lambda_j - lambda_i)

        atoms += lambda_i * alpha_i * np.exp(-lambda_i * time)

    if lambda_temp[-1] == 0.0:
        atoms = parent_atoms
    else:
        atoms = parent_atoms * atoms / lambda_temp[-1]

    # Add in the atoms to the nuclide list
    nuclides[nuclide] = atoms

    # If this is an unstable nuclide loop over each of the decays and work out the number of atoms for those daughters
    for decay_type, lambd, daughter in zip(
        decay_data_dic[nuclide]["Decay_name"][1:],
        decay_data_dic[nuclide]["lambda"][1:],
        decay_data_dic[nuclide]["Decay_daughter_names"][1:],
    ):
        # Tool cannot handle fission or very long lived >1e19 s half-lives, therefore treat these as stable.
        if decay_type == "f" or ("BR" in decay_data_dic[nuclide] and lambd < 1e-20):
            continue

        # Change the lambda to include the branching ratio for this decay
        lambda_temp[-1] = lambd

        # Add the lambda for the daugher in this decay (if it is stable we will continue as we do not need to work out the number of atoms as it does not contribute to the dose rate)
        if "lambda" in decay_data_dic[daughter] and daughter != parent:
            lambda_temp.append(decay_data_dic[daughter]["lambda"][0])
        else:
            continue

        # Irradiate this nuclide and work out the number of each nuclide in the decay chain
        nuclides_out = irradiate(
            time, parent_atoms, daughter, decay_data_dic, lambda_temp, parent
        )

        # Add the nuclides to a total list of atoms for each nuclide
        for nuclide_out in nuclides_out:
            if nuclide_out in nuclides:
                nuclides[nuclide_out] += nuclides_out[nuclide_out]
            else:
                nuclides[nuclide_out] = nuclides_out[nuclide_out]

        # Remove last lambda and move back up the chain.
        lambda_temp.pop()

    return nuclides


def get_nuclides(scenario: dict, parent: str, decay_data_dic: dict, atoms: int) -> dict:
    """Function to get the number of each nuclide from an irradiation schedule

    Parameters
    ----------
    scenario : dict
        irradiation scenario
    parent : str
        parent nuclide name
    decay_data_dic : dict
        dictionary containing decay data
    atoms : int
        number of parent atoms

    Returns
    -------
    dict
        dictionary containing the number of atoms for each nuclide in the decay chain
    """

    # Set up a temp array which is used as you move down the decay chain.
    lambda_temp = [0]

    # Set up an initial nuclides and lambda_temp array
    init_nuclides = {parent: atoms}

    # Get the initial lambda array for the parent (this is the nominal lambda values which will subsequently need to be normalised to the flux for a given pulse)
    initial_lambda = decay_data_dic[parent]["lambda"]

    # Loop over each time step in the irradiation and decay
    for time, flux in zip(scenario["times"], scenario["fluxes"]):
        # initialise a dictionary of nuclides which will occur
        decay_data_dic[parent]["lambda"] = [flux * x for x in initial_lambda]
        sum_nuclides = {}

        # Loop over each of the nuclides that is present at the start of the pulse
        for nuclide in init_nuclides:
            # Set the first lambda in the decay chain to the lambda of this nuclide
            lambda_temp[0] = decay_data_dic[nuclide]["lambda"][0]

            # Irradiate this nuclide and work out the number of each nuclide in the decay chain
            out_nuclides = irradiate(
                time,
                # flux,
                init_nuclides[nuclide],
                nuclide,
                decay_data_dic,
                lambda_temp,
                parent,
            )

            # Add the nuclides to a total list of atoms for each nuclide
            for out_nuclide in out_nuclides:
                if out_nuclide in sum_nuclides:
                    sum_nuclides[out_nuclide] += out_nuclides[out_nuclide]
                else:
                    sum_nuclides[out_nuclide] = out_nuclides[out_nuclide]

        # Set initial nuclides of next pulse to final nuclides of this pulse
        init_nuclides = sum_nuclides

    return sum_nuclides


def create_dictionary(decay_data: dict, parent: str, daughters: list) -> dict:
    """Function to create a dictionary of decay data for a given parent and daughters

    Parameters
    ----------
    decay_data : dict
        dictionary containing decay data
    parent : str
        parent nuclide name
    daughters : list
        list of daughter nuclides

    Returns
    -------
    dict
        dictionary containing decay data for the parent and daughters

    Raises
    ------
    KeyError
        if either parents or daughters are not found in the decay data
    """

    # Create decay dictionary
    decay_data_dic = {}

    # Loop over each decay data entry and add it to the dictionary
    for i in decay_data:
        decay_data_dic[i["name"]] = deepcopy(i)

        # If this is an unstable nuclide work out the decay constants, including a total lambda for decay
        if "half_life_secs" in decay_data_dic[i["name"]]:
            decay_data_dic[i["name"]]["lambda"] = [
                x * math.log(2) / decay_data_dic[i["name"]]["half_life_secs"]
                for x in decay_data_dic[i["name"]]["BR"]
            ]
            decay_data_dic[i["name"]]["lambda"].insert(
                0, sum(decay_data_dic[i["name"]]["lambda"])
            )
            decay_data_dic[i["name"]]["Decay_daughter_names"].insert(0, i["name"])
            decay_data_dic[i["name"]]["Decay_name"].insert(0, "RR")

        # If this is the parent then add the reactions and reaction rates as lambda
        if i["name"] == parent:
            decay_data_dic[parent]["Decay_daughter_names"] = [parent] + [
                daughter for daughter in daughters
            ]
            decay_data_dic[parent]["lambda"] = [
                sum([daughters[daughter] for daughter in daughters])
            ] + [daughters[daughter] for daughter in daughters]
            decay_data_dic[i["name"]]["Decay_name"] = ["RR"] * len(
                decay_data_dic[parent]["lambda"]
            )

    if parent not in decay_data_dic:
        raise KeyError("Parent {} not in decay data".format(parent))
    for daughter in daughters:
        if daughter not in decay_data_dic:
            raise KeyError("Daughter {} not in decay data".format(parent))

    return decay_data_dic


def calculate_total_activity(
    nuclide_dict: dict, irrad_scenario: dict, decay_time: int, decay_data: dict
) -> dict[str, list[np.ndarray]]:
    """Function to calculate the total activity of a given nuclide dictionary

    Parameters
    ----------
    nuclide_dict : dict
        dictionary containing nuclide data
    irrad_scenario : dict
        irradiation scenario
    decay_time : int
        decay time in seconds
    decay_data : dict
        dictionary containing decay data

    Returns
    -------
    dict[str, list[np.ndarray]]
        dictionary containing the total activity for each nuclide in the nuclide dictionary
    """

    if irrad_scenario not in IRRAD_SCENARIOS:
        add_user_irrad_scenario(irrad_scenario, IRRAD_SCENARIOS)

    # Get the irradiation scenario and decay constant for the isotope from the dictionaries
    irrad_scenario = deepcopy(IRRAD_SCENARIOS[irrad_scenario])

    # Decay times should be appended to the end of the irradiation scenario dictionary
    irrad_scenario["times"].append(decay_time)
    irrad_scenario["fluxes"].append(0)

    nuclide_dict = convert_names(nuclide_dict)

    # Store the activities for each nuclide calculated in a dictionary
    activities = {}

    for parent in nuclide_dict:
        # Convert lists to numpy arrays
        for reaction in nuclide_dict[parent]["reactions"]:
            if isinstance(nuclide_dict[parent]["reactions"][reaction], list):
                nuclide_dict[parent]["reactions"][reaction] = np.array(
                    nuclide_dict[parent]["reactions"][reaction]
                )

        # Get the maximum length of the reactions lists
        max_len = max(
            len(values) for values in nuclide_dict[parent]["reactions"].values()
        )

        for i in range(max_len):
            # Create a temporary dictionary for each reaction value
            temp_nuclide_dict = {
                parent: {
                    "atoms": nuclide_dict[parent]["atoms"],
                    "reactions": {
                        reaction: np.array(
                            [nuclide_dict[parent]["reactions"][reaction][i]]
                        )
                        for reaction in nuclide_dict[parent]["reactions"]
                        if i < len(nuclide_dict[parent]["reactions"][reaction])
                    },
                }
            }
            decay_data_dic = create_dictionary(
                decay_data, parent, temp_nuclide_dict[parent]["reactions"]
            )
            final_nuclides = get_nuclides(
                irrad_scenario,
                parent,
                decay_data_dic,
                temp_nuclide_dict[parent]["atoms"],
            )

            for nuclide in final_nuclides:
                # Only output unstable nuclides
                if "half_life_secs" in decay_data_dic[nuclide]:
                    activity = (
                        final_nuclides[nuclide] * decay_data_dic[nuclide]["lambda"][0]
                    )
                    if nuclide in activities:
                        activities[nuclide].append(activity)
                    else:
                        activities[nuclide] = [activity]
    return activities
