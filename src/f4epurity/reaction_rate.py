import numpy as np


def calculate_reaction_rate(delta_impurity, sigma_eff, flux_spectrum):

    # Check if sigma_eff has length greater than 1 (line source case)
    if isinstance(sigma_eff, np.ndarray) and len(sigma_eff) > 1:
        # Initialize a list to store the reaction rates
        reaction_rate_list = []

        # Iterate over each of the cells along the line to return a reaction rate for each
        for i in range(len(sigma_eff)):
            # Get the flux spectrum and effective cross section for the cell
            flux = flux_spectrum[i]
            sigma = sigma_eff[i]

            # Calculate the total flux
            total_flux = np.sum(flux)

            # Calculate the reaction rate
            reaction_rate = (delta_impurity / 100) * sigma * total_flux * 1e-24

            # Add the reaction rate to the list
            reaction_rate_list.append(reaction_rate)

        return reaction_rate_list

    else:  # point source case
        # Calculate the total flux
        total_flux = np.sum(flux_spectrum)

        # Calculate the reaction rate
        reaction_rate = (delta_impurity / 100) * sigma_eff * total_flux * 1e-24

        return reaction_rate
