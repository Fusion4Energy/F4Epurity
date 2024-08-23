import argparse
import datetime
import json
import os

import numpy as np
import pyvista as pv

from f4epurity.collapse import extract_flux_values, extract_xs, perform_collapse
from f4epurity.utilities import get_reactions_from_file

def write_effective_xs_map(element, filepath_flux, run_dir):

    # Read the neutron spectra VTR file
    neutron_spectra = pv.read(filepath_flux)

    effective_xs = []

    reactions = get_reactions_from_file(f'resources/xs/{element}_xs')

    # Loop over each reaction channel for the given element to output effective cross section for each
    for parent, product in reactions:

        effective_xs = []

        xs_values = extract_xs(parent, product, element)

        # Loop over each voxel in the VTR file
        for i in range(neutron_spectra.n_cells):

            # Extract the flux spectrum for every voxel (cell) in the rad map
            flux_bin_values = extract_flux_values(neutron_spectra, i)

            flux_bin_values = np.array(flux_bin_values)
            
            # Calculate the effective cross-section by collapsing with the flux
            sigma_eff = perform_collapse(xs_values, flux_bin_values)

            # Append the effective cross-section to the results list
            effective_xs.append(sigma_eff)

        effective_xs = np.array(effective_xs)

        # Add the results to the cell data of the VTR file
        neutron_spectra.cell_data[f'{parent} to {product} effective cross section'] = effective_xs

    # Write the VTR file with the new cell data
    neutron_spectra.save(f'{run_dir}/sigmaeff_NeutronFlux.vtr')

def main():
    # Define the command-line arguments
    parser = argparse.ArgumentParser(description='Collapse the cross sections available for an element with the flux spectrum')
    parser.add_argument('--element', required=True, type=str)
    parser.add_argument('--input_flux', required=True)

    args = parser.parse_args()

    # Create a unique directory for this run
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = f'output/F4Epurity-xs_{timestamp}'
    os.makedirs(run_dir, exist_ok=True)

    # Write command line arguments to metadata.json
    with open(f'{run_dir}/metadata.json', 'w') as f:
        json.dump(vars(args), f, indent=4)

    write_effective_xs_map(args.element, args.input_flux, run_dir)

if __name__ == '__main__':
    main()