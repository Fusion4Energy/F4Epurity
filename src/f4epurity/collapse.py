import logging
import numpy as np
import pyvista as pv
import sys
import os
from importlib.resources import files, as_file


def extract_flux_values(grid: pv.Grid, cell_ids: list[int]) -> np.array:
    """Extract the flux values for each energy bin at the intersected cells.
    Important assumption: arrays defining the differen energy bins should
    start with the same prefix "ValueBin-"

    Parameters
    ----------
    grid : pv.Grid
        pyvista grid object
    cell_ids : list[int]
        list of cell ids

    Returns
    -------
    np.array
        flux values for each energy bin at the extracted cells
    """
    # Check if cell_ids is a single integer, if so, convert to list
    if isinstance(cell_ids, int):
        cell_ids = [cell_ids]

    # for those cell_ids, get the flux values
    flux_spectrum = []
    for cell in cell_ids:
        flux_bin_values = []
        # Determine the number of bins dynamically - tool should be flexible to handle any bin structure
        bin_keys = [key for key in grid.cell_data.keys() if key.startswith("ValueBin-")]
        for bin_key in bin_keys:
            flux_bin = grid.cell_data[bin_key][cell]
            flux_bin_values.append(flux_bin)
        # Flux spectrum is now a list of the flux values for each energy bin for each of the intersected cells
        flux_spectrum.append(flux_bin_values)

    return np.array(flux_spectrum)


def get_flux_from_vtk(
    filepath: os.PathLike | str,
    x1: np.array,
    y1: np.array,
    z1: np.array,
    x2: np.array = None,
    y2: np.array = None,
    z2: np.array = None,
) -> np.array:
    """Extract the flux values for each energy bin at the specified location(s)
    from a VTK file. The function can handle both point and line sources.

    Parameters
    ----------
    filepath : os.PathLike | str
        path to the VTK file
    x1 : np.array
        x coordinate of the point of interest
    y1 : np.array
        y coordinate of the point of interest
    z1 : np.array
        z coordinate of the point of interest
    x2 : np.array, optional
        x coordinate of the end of the line of interest, by default None
    y2 : np.array, optional
        y coordinate of the end of the line of interest, by default None
    z2 : np.array, optional
        z coordinate of the end of the line of interest, by default None

    Returns
    -------
    np.array
        _description_
    """
    # Load the VTR file
    grid = pv.read(filepath)

    # Get the bounds of the grid (flux map)
    bounds = grid.bounds

    # Check if the input position is within the grid bounds
    if not (
        bounds[0] <= x1 <= bounds[1]
        and bounds[2] <= y1 <= bounds[3]
        and bounds[4] <= z1 <= bounds[5]
    ):
        raise ValueError("The supplied point is not within the flux map")

    # If a second point is provided, handle the line source case
    if x2 is not None and y2 is not None and z2 is not None:
        # Check if the second input position is within the grid bounds
        if not (
            bounds[0] <= x2 <= bounds[1]
            and bounds[2] <= y2 <= bounds[3]
            and bounds[4] <= z2 <= bounds[5]
        ):
            raise ValueError(
                "The supplied line source extends beyond or is not within the flux map"
            )

        # Extract the cells that are intersected by the line
        cell_ids = pv.DataSet.find_cells_intersecting_line(
            grid, [x1, y1, z1], [x2, y2, z2]
        )

        # Get the flux values for each energy bin at the intersected cells
        flux_spectrum = extract_flux_values(grid, cell_ids)

    # Point source case
    else:
        # Find the closest cell in the grid to the specified point
        cell_id = grid.find_closest_cell([x1, y1, z1])

        # Get the flux values for each energy bin at the cell closest to the location of interest
        flux_spectrum = extract_flux_values(grid, cell_id)

    return np.array(flux_spectrum)


def extract_xs(parent: str, product: str, element: str) -> list[float]:
    """Extract the cross section data for a given parent-daughter pair and element.

    Parameters
    ----------
    parent : str
        parent isotope name
    product : str
        daughter isotope name
    element : str
        element name

    Returns
    -------
    list[float]
        available cross section data
    """
    xs_values = []
    # Get the file path to the xs data
    with as_file(
        files("f4epurity.resources.xs").joinpath(f"{element}_xs")
    ) as xs_data_path:
        # Read the text file containing the xs data
        with open(xs_data_path, "r", encoding="utf-8") as file:
            lines = file.readlines()
    found = False
    for line in lines:
        split_line = line.split()
        if len(split_line) == 0:
            # If found the correct cross section, stop collecting values
            if found:
                break
        # Check found the correct cross section
        elif found:
            if len(split_line) == 3:
                xs_values.append(float(split_line[2]))
        elif len(split_line) == 3:
            # Look for the parent and daughter isotopes requested
            reaction = split_line[0] + " " + split_line[2]
            if reaction == parent + " " + product:
                found = True
    return xs_values


# Mathematical operation for collapsing the cross section with the flux
def perform_collapse(xs_values: np.array, flux: np.array) -> np.array:
    """Perform the collapse operation to calculate the effective cross section.

    Parameters
    ----------
    xs_values : np.array
        array of xs values for each energy bin
    flux : np.array
        array of flux values for each energy bin

    Returns
    -------
    np.array
        array of effective cross section values
    """
    xs_group = np.array(xs_values)

    total_flux = np.sum(flux)

    # Check the flux spectrum is non-zero
    if total_flux == 0:
        logging.warning(
            "Warning: The flux is zero in all energy bins at the selected location."
        )

    # Calculate the effective cross section
    sigma_eff = np.sum(xs_group * flux) / total_flux

    return sigma_eff


# Function to collapse the flux spectrum with the cross section data
def collapse_flux(
    xs_values: np.array,
    filepath_flux: str | os.PathLike,
    x1: np.array,
    y1: np.array,
    z1: np.array,
    x2: np.array = None,
    y2: np.array = None,
    z2: np.array = None,
) -> tuple[np.array, np.array]:
    """Collapse the flux spectrum with the cross section data to calculate
    the effective cross section.

    Parameters
    ----------
    xs_values : np.array
        array of xs values for each energy bin
    filepath_flux : str | os.PathLike
        path to the VTK file containing the flux spectrum
    x1 : np.array
        x coordinate of the point of interest
    y1 : np.array
        y coordinate of the point of interest
    z1 : np.array
        z coordinate of the point of interest
    x2 : np.array, optional
        x coordinate of end of line, by default None
    y2 : np.array, optional
        y coordinate of end of line, by default None
    z2 : np.array, optional
        z coordinate of end of line, by default None

    Returns
    -------
    tuple[np.array, np.array]
        sigma_eff, flux_spectrum
    """

    # Load the flux spectrum from the VTK file
    flux_spectrum = get_flux_from_vtk(filepath_flux, x1, y1, z1, x2, y2, z2)

    # Check if flux_spectrum is a 2D array (line source case)
    # if flux_spectrum.ndim == 2:
    sigma_eff_list = [perform_collapse(xs_values, flux) for flux in flux_spectrum]
    sigma_eff = np.array(sigma_eff_list)

    # # If flux_spectrum is a 1D array (point source case)
    # else:
    #     sigma_eff = perform_collapse(xs_values, flux_spectrum)
    #     print(sigma_eff)
    #     # get data type
    #     print(type(sigma_eff))

    return sigma_eff, flux_spectrum
