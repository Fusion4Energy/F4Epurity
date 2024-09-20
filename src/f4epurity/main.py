from jsonargparse import Namespace
import csv
import datetime
import json
import numpy as np
import os
from importlib.resources import files, as_file
import pandas as pd

from f4epurity.decay_chain_calc import calculate_total_activity
from f4epurity.collapse import collapse_flux, extract_xs
from f4epurity.dose import convert_to_dose, write_vtk_file, plot_slice
from f4epurity.maintenance import (
    dose_within_workstation,
    get_dose_at_workstation,
    read_maintenance_locations,
)
from f4epurity.reaction_rate import calculate_reaction_rate
from f4epurity.utilities import (
    calculate_number_of_atoms,
    get_isotopes,
    sum_vtr_files,
    normalise_nuclide_name,
    get_reactions_from_file,
)
from f4epurity.parsing import parse_arguments, parse_isotopes_activities_file


# Main function
def calculate_dose_for_source(
    args: Namespace,
    x1: np.ndarray,
    y1: np.ndarray,
    z1: np.ndarray,
    run_dir: str | os.PathLike,
    nist_df: pd.DataFrame,
    reactions: set[tuple[str, str]] | None,
    decay_data: dict,
    dose_factors_df: pd.DataFrame,
    x2: np.ndarray = None,
    y2: np.ndarray = None,
    z2: np.ndarray = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[float]]:
    """Calculate the dose for a given source

    Parameters
    ----------
    args : Namespace
        parsed arguments
    x1 : np.ndarray
        x location(s) of the source(s)
    y1 : np.ndarray
        y location(s) of the source(s)
    z1 : np.ndarray
        z location(s) of the source(s)
    run_dir : str | os.PathLike
        direcctory where the code is run
    nist_df : pd.DataFrame
        dataframed NIST data
    reactions : set[tuple[str, str]] | None
        set of reactions
    decay_data : dict
        decay data
    dose_factors_df : pd.DataFrame
        dose factors dataframe
    x2 : np.ndarray, optional
        second x location(s) of the line source(s), by default None
    y2 : np.ndarray, optional
        second y location(s) of the line source(s) , by default None
    z2 : np.ndarray, optional
        second z location(s) of the line source(s), by default None

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[float]]
        dose_array, x, y, z, total_dose
    """
    if args.activities_file:
        activities = parse_isotopes_activities_file(args.activities_file)
    else:
        # Expand input element to natural isotopes
        isotopes = get_isotopes(args.element, nist_df)

        # Dictionary to store reaction rates
        reaction_rates = {}

        print("Performing Collapse and Calculating Reaction Rates...")
        # Populate the dictionary with the reaction rates for each possible reaction channel for a given element
        for parent, product in reactions:
            if parent not in isotopes:
                continue

            if parent not in reaction_rates:
                # Determine the number of atoms of the parent nuclide
                number_of_atoms = calculate_number_of_atoms(
                    parent, args.delta_impurity, nist_df
                )

                # Add the parent to the dictionary with the number of atoms and an empty dictionary for the reaction rates
                reaction_rates[parent] = {"atoms": number_of_atoms, "reactions": {}}

            xs_values = extract_xs(parent, product, args.element)

            # Get the flux value from the VTK file and collapse with the 5-group cross section
            sigma_eff, flux_spectrum = collapse_flux(
                xs_values, args.input_flux, x1, y1, z1, x2, y2, z2
            )

            # Calculate the reaction rate based on the flux and effective cross section
            reaction_rate = calculate_reaction_rate(
                args.delta_impurity, sigma_eff, flux_spectrum
            )

            # Store the reaction rate in the dictionary
            reaction_rates[parent]["reactions"][product] = reaction_rate

        # Call the decay_chain_calculator to determine the activity of each nuclide
        print("Calculating Activities...")
        activities = calculate_total_activity(
            reaction_rates, args.irrad_scenario, args.decay_time, decay_data
        )
    # Initialize a list to store the total dose for each element
    total_dose = None

    print("Calculating the Dose...")
    # Determine the Dose for each nuclide
    for nuclide, nuclide_activity in activities.items():

        # Convert to format in dose conversion spreadsheet
        nuclide = normalise_nuclide_name(nuclide)

        # Convert the activity to a dose
        doses = convert_to_dose(nuclide, nuclide_activity, dose_factors_df)

        # Check if doses is a list
        if isinstance(doses, list):
            # If total_dose is None, set it to the doses list
            if total_dose is None:
                total_dose = doses
            else:
                # Otherwise, add the corresponding elements of the doses list to total_dose
                total_dose = [total + dose for total, dose in zip(total_dose, doses)]
                total_dose = [dose.item() for dose in total_dose]
        else:
            # If dose is a single number, add it to the total dose directly
            if total_dose is None:
                total_dose = [doses]
            else:
                total_dose = [total + doses for total in total_dose]

    print("Writing the Dose Map...")
    # Write the dose array and output to a VTR file
    dose_array, x, y, z, plot_bounds = write_vtk_file(
        total_dose,
        x1,
        y1,
        z1,
        run_dir,
        args.input_flux,
        x2,
        y2,
        z2,
        args.output_all_vtr,
    )

    # Run quick plot function to output png image
    if args.plot:
        plt = plot_slice(
            dose_array, x, y, z, args.plot[0], int(args.plot[1]), plot_bounds
        )
        if x2 is not None and y2 is not None and z2 is not None:
            plt.savefig(f"{run_dir}/dose_{x1}_{y1}_{z1}_to_{x2}_{y2}_{z2}.png")
        else:
            plt.savefig(f"{run_dir}/dose_{x1}_{y1}_{z1}.png")

    return dose_array, x, y, z, total_dose


def calculate_dose_at_workstations(
    args: Namespace, dose, x1, y1, z1, run_dir, x2=None, y2=None, z2=None
):
    # If dose is to be calculated at workstations
    if args.workstation:
        # Get the coordinates and workstation names for the given workstation and location
        coordinates_list, workstation_list = read_maintenance_locations(
            args.workstation, args.location
        )

        # Create the output directory if it doesn't exist
        os.makedirs("output", exist_ok=True)

        # Create file with the dose values for each workstation
        filename = (
            f"{run_dir}/dose_{x1}_{y1}_{z1}_{args.location}.csv"
            if x2 is None
            else f"{run_dir}/dose_{x1}_{y1}_{z1}_to_{x2}_{y2}_{z2}_{args.location}.csv"
        )
        with open(filename, "w", newline="", encoding="utf8") as f:
            writer = csv.writer(f)

            writer.writerow(["Workstation", "Delta Dose (micro Sieverts per hour)"])

            # Loop over the list of coordinates and workstation names and calculate the dose for each set of coordinates
            for coordinates, workstation in zip(coordinates_list, workstation_list):
                min_x, max_x, min_y, max_y, min_z, max_z = coordinates

                max_dose, max_dose_coord = get_dose_at_workstation(
                    dose,
                    (x1, y1, z1),
                    (x2, y2, z2),
                    min_x,
                    max_x,
                    min_y,
                    max_y,
                    min_z,
                    max_z,
                    x2 is not None,
                )

                # Convert the list to a string
                max_dose_str = "{:.3e}".format(max_dose[0] if x2 is None else max_dose)

                # Write the output to the CSV file
                writer.writerow([workstation, max_dose_str])


def process_sources(args: Namespace) -> None:
    # Create a unique directory for this run
    root = args.root_output
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = f"{root}/F4Epurity_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)

    # Write command line arguments to metadata.json
    with open(f"{run_dir}/metadata.json", "w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=4)

    # Conditionally load files based on the presence of activities_file
    if not args.activities_file:
        # Read the necessary files once
        nist_file_path = files("f4epurity.resources").joinpath("NIST_tabulated.xlsx")
        with as_file(nist_file_path) as fp:
            nist_df = pd.read_excel(fp)

        xs_file_path = files("f4epurity.resources.xs").joinpath(f"{args.element}_xs")
        with as_file(xs_file_path) as fp:
            reactions = get_reactions_from_file(fp)

        decay_data_path = files("f4epurity.resources").joinpath("Decay2020.json")
        with as_file(decay_data_path) as fp:
            with open(fp, "r", encoding="utf-8") as json_file:
                decay_data = json.load(json_file)
    else:
        # If activities_file is provided, skip loading unnecessary files
        nist_df = None
        reactions = None
        decay_data = None

    # Load the dose matrix file (always needed)
    dose_matrix_file_path = files("f4epurity.resources").joinpath("F4E_dosematrix.xlsx")
    with as_file(dose_matrix_file_path) as fp:
        dose_factors_df = pd.read_excel(fp)

    dose_arrays = []
    # Check if a second point was provided - line source
    if args.x2 is not None and args.y2 is not None and args.z2 is not None:
        # Line source
        print("Line source(s) selected")

        # Handle multiple sets of coordinates being provided (multiple line sources)
        for x1, y1, z1, x2, y2, z2 in zip(
            args.x1, args.y1, args.z1, args.x2, args.y2, args.z2
        ):
            dose_array, x, y, z, dose = calculate_dose_for_source(
                args,
                x1,
                y1,
                z1,
                run_dir,
                nist_df,
                reactions,
                decay_data,
                dose_factors_df,
                x2,
                y2,
                z2,
            )
            dose_arrays.append(dose_array)
            calculate_dose_at_workstations(
                args,
                dose,
                x1,
                y1,
                z1,
                run_dir,
                x2,
                y2,
                z2,
            )
    else:
        # Point source
        print("Point source(s) selected")

        # Handle multiple coordinates being provided
        for x1, y1, z1 in zip(args.x1, args.y1, args.z1):
            dose_array, x, y, z, dose = calculate_dose_for_source(
                args,
                x1,
                y1,
                z1,
                run_dir,
                nist_df,
                reactions,
                decay_data,
                dose_factors_df,
            )
            dose_arrays.append(dose_array)
            calculate_dose_at_workstations(args, dose, x1, y1, z1, run_dir)

    # If more than one dose array is present, sum the dose arrays (multiple sources)
    if len(dose_arrays) > 1:
        sum_vtr_files(dose_arrays, x, y, z, run_dir, masses=args.m)
        if args.workstation:
            # Get the coordinates and workstation names for the given workstation and location
            coordinates_list, workstation_list = read_maintenance_locations(
                args.workstation, args.location
            )

            with open(
                f"{run_dir}/dose_{args.location}_total.csv",
                "w",
                newline="",
                encoding="utf-8",
            ) as f:
                writer = csv.writer(f)
                writer.writerow(["Workstation", "Dose"])

                # Loop over the list of coordinates and workstation names and calculate the dose for each set of coordinates
                for coordinates, workstation in zip(coordinates_list, workstation_list):
                    min_x, max_x, min_y, max_y, min_z, max_z = coordinates
                    box_bounds = (min_x, max_x, min_y, max_y, min_z, max_z)
                    max_dose = dose_within_workstation(
                        f"{run_dir}/dose_total.vtr", box_bounds
                    )
                    max_dose_str = "{:.3e}".format(max_dose)
                    writer.writerow([workstation, max_dose_str])


def main(args_list: list[str] | None = None):
    args = parse_arguments(args_list)
    process_sources(args)


if __name__ == "__main__":
    main()
