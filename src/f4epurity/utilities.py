import os
import re
import sys

import numpy as np
import pandas as pd
import pyevtk


def get_isotopes(element, nist_df):

    # Find element and only retrieve those occuring naturally
    isotopes_df = nist_df[
        (nist_df["Atomic Symbol"].str.strip() == element.lower())
        & (pd.notna(nist_df["Isotopic Composition"]))
    ]

    # Concatenate the 'Atomic Symbol' and 'Mass Number' columns
    isotopes = (
        isotopes_df["Atomic Symbol"].str.strip()
        + isotopes_df["Mass Number"].astype(str)
    ).tolist()

    return isotopes


def calculate_number_of_atoms(nuclide, delta_impurity, nist_df):

    # Get the molar mass for the given isotope
    molar_mass = get_molar_mass(nuclide, nist_df)

    # Calculate the number of atoms in the sample
    number_of_atoms = (delta_impurity / 100) * 6.02214076e23 / molar_mass

    return number_of_atoms


def get_molar_mass(isotope, nist_df):

    # Remove leading and trailing spaces from 'Atomic Symbol' column
    nist_df["Atomic Symbol"] = nist_df["Atomic Symbol"].str.strip()

    # Use a regular expression to split the element string into the atomic symbol and atomic number parts
    match = re.match(r"([a-zA-Z]+)(\d+)", isotope, re.I)

    items = match.groups()

    atomic_symbol = items[0]
    mass_number = int(items[1])

    # Filter the DataFrame to include only rows for the given atomic symbol and atomic number
    df_filtered = nist_df[
        (nist_df["Atomic Symbol"] == atomic_symbol)
        & (nist_df["Mass Number"] == mass_number)
    ]

    # If no matching row is found, return an error
    if df_filtered.empty:
        print(
            f"Error: No matching atomic symbol and mass number found for {atomic_symbol}, {mass_number}"
        )
        sys.exit(1)

    # Get the molar mass for the given atomic symbol and atomic number
    molar_mass = float(df_filtered["Standard Atomic Weight"].values[0])

    return molar_mass


def sum_vtr_files(dose_arrays, x, y, z, run_dir):
    # Sum the dose arrays from multiple point/line sources
    sum_dose = np.sum(dose_arrays, axis=0)

    os.makedirs("output", exist_ok=True)

    # Write the summed dose to a VTK file
    pyevtk.hl.gridToVTK(
        f"{run_dir}/dose_total", x, y, z, cellData={"Dose_Total": sum_dose}
    )


def get_reactions_from_file(filename):
    reactions = set()
    with open(filename, "r") as file:
        for line in file:
            if "(" in line:  # This characterises specific reaction line
                reaction = re.split(r"\s+|\(|\)", line)
                parent = reaction[0]
                product = reaction[4]
                reactions.add((parent, product))
    return reactions


def normalise_nuclide_name(nuclide):
    # Split the nuclide name into the element name, the mass number, and 'm' or 'n' if present
    match = re.match(r"([A-Za-z]+)(\d+)([mn]?)", nuclide)
    element_name, mass_number, suffix = match.groups()

    # Capitalize the first letter of the element name
    element_name = element_name.capitalize()

    # Remove leading zeros from the mass number
    mass_number = mass_number.lstrip("0")

    return element_name + mass_number + suffix


# Return the long nuclide name
def get_name(parent):
    parts = re.split("(\d+)", parent)
    if len(parts[1]) < 3:
        iso = "0" * (3 - len(parts[1])) + parts[1]
    else:
        iso = parts[1]

    name = parts[0].capitalize() + iso + parts[2]

    return name


# Convert names of dictionary to long nuclide names
def convert_names(nuc_dict):

    new_dict = {}
    for parent in nuc_dict:
        new_name = get_name(parent)
        new_dict[new_name] = {}
        new_dict[new_name]["atoms"] = nuc_dict[parent]["atoms"]
        new_dict[new_name]["reactions"] = {}
        for daughter in nuc_dict[parent]["reactions"]:
            new_dict[new_name]["reactions"][get_name(daughter)] = nuc_dict[parent][
                "reactions"
            ][daughter]

    return new_dict


def add_user_irrad_scenario(filename, irrad_scenarios):
    day_to_sec = 24 * 60 * 60

    df = pd.read_csv(filename, sep="\s+", header=None, names=["times", "fluxes"])

    df["times"] = df["times"] * day_to_sec

    # Add the new irradiation scenario to the dictionary
    irrad_scenarios[filename] = {
        "times": df["times"].tolist(),
        "fluxes": df["fluxes"].tolist(),
    }
