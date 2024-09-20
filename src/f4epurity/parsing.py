from typing import List, Optional

import os
import numpy as np
from jsonargparse import ArgumentParser, Namespace, ActionConfigFile
import pandas as pd


def parse_arguments(args_list: Optional[List[str]] = None) -> Namespace:
    """Parse the command-line arguments for the tool.

    Parameters
    ----------
    args_list : Optional[List[str]], optional
        list of arguments for the parser, by default None

    Returns
    -------
    Namespace
        attributes are the command line arguments
    """
    # Define the command-line arguments for the tool
    parser = ArgumentParser(
        description="Approximate the deviation in activity and dose rate based on local DR/NCR"
    )
    # put the config file first so that other command lines can override its contents
    parser.add_argument("--cfg", action=ActionConfigFile, help="path to config file")
    # for testing purposes it is nice to allow an optional argument to change
    # the default run directory
    parser.add_argument(
        "--root_output",
        type=str,
        default="output",
        help="Root directory for output files",
    )

    # Add option to provide a list of isotopes and their activities
    parser.add_argument(
        "--activities_file",
        type=str,
        help="Path to a text file containing isotopes and their pre-computed activities",
    )

    # Add other arguments here
    parser.add_argument(
        "--element",
        type=str,
        help="Element with deviation e.g. Co, Ta, Nb",
    )
    parser.add_argument(
        "--delta_impurity",
        type=float,
        help="Deviation in the elements weight as a percentage",
    )
    parser.add_argument(
        "--input_flux",
        help="Provide a path to the vtr file with the neutron spectrum",
    )
    parser.add_argument(
        "--irrad_scenario",
        type=str,
        help="Irradiation scenario. SA2 and DT1 are available or supply path to file with user defined scenario",
    )
    parser.add_argument(
        "--decay_time",
        type=float,
        help="Decay time for calculating dose in seconds",
    )
    parser.add_argument(
        "--sources_csv",
        type=str,
        help="CSV file containing coordinates for point/line sources. Columns: x1, y1, z1, (x2, y2, z2)",
    )
    parser.add_argument(
        "--x1",
        nargs="+",
        type=float,
        help="x coordinate of point source",
    )
    parser.add_argument(
        "--y1",
        nargs="+",
        type=float,
        help="y coordinate of point source",
    )
    parser.add_argument(
        "--z1",
        nargs="+",
        type=float,
        help="z coordinate of point source",
    )
    parser.add_argument(
        "--x2",
        nargs="+",
        type=float,
        help="x coordinate of the second point for line source",
    )
    parser.add_argument(
        "--y2",
        nargs="+",
        type=float,
        help="y coordinate of the second point for line source",
    )
    parser.add_argument(
        "--z2",
        nargs="+",
        type=float,
        help="z coordinate of the second point for line source",
    )
    parser.add_argument(
        "--m",
        nargs="+",
        type=float,
        help="Mass of the component where the impurity is located in g",
        default=None,
    )
    parser.add_argument(
        "--plot",
        nargs=2,
        metavar=("slice_axis", "slice_location"),
        help="Output image of dose mesh at a given location",
    )
    parser.add_argument(
        "--workstation",
        type=str,
        help="Name of the workstation for which to report the max dose e.g. 1, 3, 4 or 'all'",
    )
    parser.add_argument(
        "--location", type=str, help="Location of the workstation(s) e.g. Nb cell"
    )
    parser.add_argument(
        "--output_all_vtr",
        action="store_true",
        help="Output individual VTR files for each source",
    )

    # Parse the arguments
    args = parser.parse_args(args_list)

    # Manually check required arguments if activities_file is not provided
    if not args.activities_file:
        required_args = [
            "element",
            "delta_impurity",
            "input_flux",
            "irrad_scenario",
            "decay_time",
        ]
        missing_args = [arg for arg in required_args if getattr(args, arg) is None]
        if missing_args:
            parser.error(f"Missing required arguments: {', '.join(missing_args)}")

    # If a location is specified a workstation must also be given and vice versa
    if (args.workstation is None) != (args.location is None):
        parser.error("--workstation and --location must be supplied together")

    # only one between sources_csv and x1, y1, z1 (or m) should be provided
    if args.sources_csv and (args.x1 or args.y1 or args.z1 or args.m):
        parser.error("--sources_csv and --x1, --y1, --z1 are mutually exclusive")
    if not args.sources_csv and not (args.x1 and args.y1 and args.z1):
        parser.error("One between --sources_csv and --x1, --y1, --z1 must be provided")

    # Set output_all_vtr to True if only one set of coordinates is given for point source
    if (
        args.x1
        and len(args.x1) == 1
        and args.y1
        and len(args.y1) == 1
        and args.z1
        and len(args.z1) == 1
    ):
        args.output_all_vtr = True

    # Set output_all_vtr to True if only one line source is specified
    if (
        args.x2
        and len(args.x2) == 1
        and args.y2
        and len(args.y2) == 1
        and args.z2
        and len(args.z2) == 1
    ):
        args.output_all_vtr = True

    del args.cfg  # to avoid issues down the line, job has been done

    # Check if a CSV file was provided
    args = _validate_source_coordinates_input(args)

    # the number of coordinates provided must be the same (and mass)
    # first of all assess that the coordinates have the same length
    if not (len(args.x1) == len(args.y1) == len(args.z1)):
        parser.error("The number of coordinates x1, y1 and z1 must match")
    if args.m is not None:
        if not len(args.m) == len(args.x1):
            parser.error("The number of coordinates must match the number of masses")
    if args.x2 is not None:
        if not (len(args.x2) == len(args.y2) == len(args.z2) == len(args.x1)):
            parser.error("The number of coordinates x2, y2 and z2 must match")

    return args


# For user supplied list of activities, convert the text file to a dictionary
def parse_isotopes_activities_file(
    file_path: str | os.PathLike,
) -> dict[str, np.ndarray]:
    """Parse a text file containing isotopes and their activities into a dictionary.

    Parameters
    ----------
    file_path : str | os.PathLike
        Path to the text file containing isotopes and their activities.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing isotopes and their activities.
    """
    activities = {}
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            isotope, activity = line.strip().split()
            activities[isotope] = [np.array([float(activity)])]
    return activities


def _validate_source_coordinates_input(args: Namespace) -> Namespace:
    if args.sources_csv:
        try:
            # Read the CSV file
            coordinates = pd.read_csv(args.sources_csv)
            args.x1 = coordinates["x1"].tolist()
            args.y1 = coordinates["y1"].tolist()
            args.z1 = coordinates["z1"].tolist()
            # also the second point may be provided
            if len(coordinates.columns) > 4:
                args.x2 = coordinates["x2"].tolist()
                args.y2 = coordinates["y2"].tolist()
                args.z2 = coordinates["z2"].tolist()
        except KeyError as e:
            raise KeyError(
                "CSV file must contain columns 'x1', 'y1', 'z1' and optionally 'x2', 'y2', 'z2'"
            ) from e

        # masses may also have been provided here
        if "m" in coordinates.columns:
            args.m = coordinates["m"].tolist()

    return args
