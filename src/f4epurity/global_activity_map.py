import argparse
import datetime
import copy
import importlib.resources as pkg_resources
import json
import numpy as np
import pandas as pd
import os
import pyvista as pv

from f4epurity.collapse import extract_flux_values, perform_collapse, extract_xs
from f4epurity.reaction_rate import calculate_reaction_rate
from f4epurity.utilities import get_reactions_from_file, calculate_number_of_atoms
from f4epurity.decay_chain_calc import calculate_total_activity


def write_activity_map(element, filepath_flux, delta_impurity, decay_time, irrad_scenario, run_dir):
        
    # Read the tabulated NIST data
    nist_df = pd.read_excel('resources/NIST_tabulated.xlsx')

    # Read the possible reaction channels for the given element
    xs_file_path = pkg_resources.path('f4epurity.resources.xs', f'{element}_xs')

    with xs_file_path as fp:
        reactions = get_reactions_from_file(fp)

    # Read the neutron spectra VTR file
    neutron_spectra = pv.read(filepath_flux)

    # Load the decay data from the json file
    file = open('resources/Decay2020.json')
    decay_data = json.load(file)

    reaction_rates = {}
    activities = {}

    # Loop over each voxel in the VTR file
    for i in range(neutron_spectra.n_cells):

        reaction_rates = {}

        # Loop over each reaction channel for the given element to output effective cross section for each
        for parent, product in reactions:

            xs_values = extract_xs(parent, product, element)

            if parent not in reaction_rates:
                # Determine the number of atoms of the parent nuclide
                number_of_atoms = calculate_number_of_atoms(parent, delta_impurity, nist_df)

                # Add the parent to the dictionary with the number of atoms and an empty dictionary for the reaction rates
                reaction_rates[parent] = {"atoms": number_of_atoms, "reactions": {}}

            # Extract the flux spectrum for every voxel (cell) in the rad map
            flux_bin_values = extract_flux_values(neutron_spectra, i)

            flux_bin_values = np.array(flux_bin_values)
            
            # Calculate the effective cross-section by collapsing with the flux
            sigma_eff = perform_collapse(xs_values, flux_bin_values)

            # Calculate the reaction rate based on the flux and effective cross section
            reaction_rate = calculate_reaction_rate(delta_impurity, sigma_eff, flux_bin_values)

            reaction_rates[parent]["reactions"][product] = reaction_rate

            print (reaction_rates)

        # Calculate the total activity for this voxel
        activity = calculate_total_activity(reaction_rates, irrad_scenario, decay_time, decay_data)

        print(activity)

        # Append the activity to the list for each nuclide
        for nuclide, value in activity.items():
            if nuclide not in activities:
                activities[nuclide] = []
            activities[nuclide].append(value)

    # Add it to the neutron_spectra as cell data
    for nuclide, values in activities.items():
        neutron_spectra.cell_data[f'{nuclide} Bq per gram'] = np.array(values)

    neutron_spectra.save(f'{run_dir}/activity_NeutronFlux.vtr')

def main():
    # Define the command-line arguments
    parser = argparse.ArgumentParser(description='Produce a global map of the activity each nuclide')
    parser.add_argument('--element', required=True, type=str)
    parser.add_argument('--delta_impurity', required=True, type=float)
    parser.add_argument('--input_flux', required=True)
    parser.add_argument('--irrad_scenario', required=True, type=str)
    parser.add_argument('--decay_time', required=True, type=float)

    args = parser.parse_args()

    # Create a unique directory for this run
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = f'output/F4Epurity-activity_{timestamp}'
    os.makedirs(run_dir, exist_ok=True)

    # Write command line arguments to metadata.json
    with open(f'{run_dir}/metadata.json', 'w') as f:
        json.dump(vars(args), f, indent=4)

    write_activity_map(args.element, args.input_flux, args.delta_impurity, args.decay_time, args.irrad_scenario, run_dir)

if __name__ == '__main__':
    main()