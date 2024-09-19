import os
import sys
from math import pi

from importlib.resources import files, as_file
import matplotlib.pyplot as plt
import numpy as np
import pyevtk
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter

from f4epurity.stl_plot import STL


# Extract the dose factors from the database and return the dose factor for the given nuclide
def extract_dose_factors(nuclide, df):
    # Find the row for the given nuclide
    nuclide_row = df[df.iloc[:, 0] == nuclide]

    # Capture if the specified nuclide is not found in the database
    if nuclide_row.empty:
        return f"Dose conversion factor not found for nuclide {nuclide}."

    # Return the dose factor for the nuclide
    return nuclide_row.iloc[:, 1].values[
        0
    ]  # assuming dose factor is in the second column


# Function to convert the input values to a dose
def convert_to_dose(nuclide, activity, df):
    dose_conversion_factor = extract_dose_factors(nuclide, df)

    if hasattr(activity, "__iter__") and len(activity) > 1:
        # Initialize a list to store the doses
        dose_list = []

        # Iterate over each of the activities along the line to return a dose for each
        for i in range(len(activity)):
            # Get the activity for the point
            act = activity[i]

            # Calculate the dose for the point - convert to microsieverts
            # User to fold in total mass of material
            dose = dose_conversion_factor * act * 1e6

            # Add the dose to the list
            dose_list.append(dose)

        return dose_list

    else:  # point source case
        # Calculate the dose for the point - convert to microsieverts
        # User to fold in total mass of material
        dose = dose_conversion_factor * activity[0].item() * 1e6
        return dose


def dose_from_line_source(dose, x1, y1, z1, x2, y2, z2, x, y, z):
    # Mathematical model for the dose calculation from a line source
    # Dose_at_point = Dose * theta / w, where theta is the angle subtended by the line source at the point of interest and w the perpendicular distance from the line source to the point of interest

    # Convert coordinates to numpy arrays
    source_start = np.array([x1, y1, z1])
    source_end = np.array([x2, y2, z2])
    point = np.array([x, y, z])

    # Calculate the vectors from the point to the ends of the line source
    vector_to_start = source_start - point
    vector_to_end = source_end - point

    # Calculate the angle subtended by the line at the point
    cos_angle = np.dot(vector_to_start, vector_to_end) / (
        np.linalg.norm(vector_to_start) * np.linalg.norm(vector_to_end)
    )
    cos_angle = np.clip(cos_angle, -1, 1)  # To avoid errors due to numerical precision
    angle_subtended = np.arccos(cos_angle)

    # Calculate the perpendicular distance from the point to the line
    distance_perpendicular = np.linalg.norm(
        np.cross(source_end - source_start, point - source_start)
    ) / np.linalg.norm(source_end - source_start)

    # Calculate the length of the line source
    length_source = np.linalg.norm(source_start - source_end)
    dose_per_length = sum(dose) / length_source

    # Calculate the dose at the point from the line source
    if distance_perpendicular == 0:
        dose_at_point = dose_per_length
    else:
        dose_at_point = (
            dose_per_length * angle_subtended / (distance_perpendicular * 4 * pi)
        )

    return dose_at_point


def is_within_bounds(x1, y1, z1, x2=None, y2=None, z2=None):

    # Path to stl files for ITER B1 rooms
    folder_path = files("f4epurity.resources").joinpath("building_stl_files")

    # Get a list of all STL files in the folder
    stl_files = [f for f in folder_path.iterdir() if f.suffix == ".stl"]

    # Initialize a list to store the bounds of the STLs
    bounds_list = []

    # Loop over all STL files for ITER rooms
    for stl_file in stl_files:
        # Load the STL file
        stl = pv.read(str(stl_file))

        # Get the bounds of the STL file
        bounds = list(stl.bounds)

        # Increase the bounds by 50 cm for the purposes of plotting (plot not clipped exactly to edge of room)
        for i in range(0, len(bounds), 2):
            bounds[i] -= 50  # min values
            bounds[i + 1] += 50  # max values

        # Check if both sets of coordinates are given - line source case
        if x2 is not None and y2 is not None and z2 is not None:
            # Check if both points lie within the bounds
            if (
                bounds[0] <= x1 <= bounds[1]
                and bounds[2] <= y1 <= bounds[3]
                and bounds[4] <= z1 <= bounds[5]
                and bounds[0] <= x2 <= bounds[1]
                and bounds[2] <= y2 <= bounds[3]
                and bounds[4] <= z2 <= bounds[5]
            ):
                bounds_list.extend(bounds)
                bounds_list.append(stl_file)
        else:
            # Check if the point lies within the bounds - point source case
            if (
                bounds[0] <= x1 <= bounds[1]
                and bounds[2] <= y1 <= bounds[3]
                and bounds[4] <= z1 <= bounds[5]
            ):
                bounds_list.extend(bounds)
                bounds_list.append(stl_file)

    # Check if the bounds_list is empty
    if not bounds_list:
        # If the list is empty, the points are not within any of the files
        print(
            "WARNING: The point(s) of the source are not within the bounds of any of the available ITER room stl files."
        )
        return None

    # Print the name of the stl file in which the points are located
    # print(f"The source is within {bounds_list[6]}")

    return bounds_list


# Function to write the dose values to a vtk file
def write_vtk_file(
    dose, x1, y1, z1, run_dir, x2=None, y2=None, z2=None, output_all_vtr=False
):

    # Get the bounds in which to make the plot
    plot_bounds = is_within_bounds(x1, y1, z1, x2, y2, z2)

    # If the coordinates are outside one of the stls, use an arbitrary volume for writing the mesh
    if plot_bounds is None:
        if x2 is not None and y2 is not None and z2 is not None:  # Line source case
            x_center = (x1 + x2) / 2
            y_center = (y1 + y2) / 2
            z_center = (z1 + z2) / 2

            plot_bounds = [
                x_center - 500,
                x_center + 500,
                y_center - 500,
                y_center + 500,
                z_center - 500,
                z_center + 500,
            ]
        else:  # Point source case
            plot_bounds = [x1 - 500, x1 + 500, y1 - 500, y1 + 500, z1 - 500, z1 + 500]

    # Calculate the number of steps for x, y, and z
    num_steps_x = int((plot_bounds[1] - plot_bounds[0]) / 50)
    num_steps_y = int((plot_bounds[3] - plot_bounds[2]) / 50)
    num_steps_z = int((plot_bounds[5] - plot_bounds[4]) / 50)

    # Generate values from plot_bounds[0] to plot_bounds[1] with num_steps steps
    x = np.linspace(plot_bounds[0], plot_bounds[1], num=num_steps_x)
    y = np.linspace(plot_bounds[2], plot_bounds[3], num=num_steps_y)
    z = np.linspace(plot_bounds[4], plot_bounds[5], num=num_steps_z)

    # Create a 3D numpy array to store the dose values
    dose_array = np.zeros((len(x) - 1, len(y) - 1, len(z) - 1))

    # Check if dose is a list with more than one value (line source case)
    if hasattr(dose, "__iter__") and len(dose) > 1:
        # Iterate over each point in the grid
        for i in range(len(x) - 1):
            for j in range(len(y) - 1):
                for k in range(len(z) - 1):
                    # Calculate the dose at the point from the line source
                    dose_point = dose_from_line_source(
                        dose, x1, y1, z1, x2, y2, z2, x[i], y[j], z[k]
                    )
                    dose_array[i, j, k] = dose_point

    # Create a 3D numpy array filled with the dose values based on 1/r^2
    else:
        for i in range(len(x) - 1):
            for j in range(len(y) - 1):
                for k in range(len(z) - 1):
                    distance = np.sqrt(
                        (x[i] - x1) ** 2 + (y[j] - y1) ** 2 + (z[k] - z1) ** 2
                    )
                    # Avoid division by zero
                    if distance == 0:
                        dose_array[i, j, k] = dose[0]
                    else:
                        dose_array[i, j, k] = dose[0] / (4 * pi * distance**2)

    # Get the indices of the max value
    indices = np.unravel_index(np.argmax(dose_array, axis=None), dose_array.shape)

    # Get the coordinates of the max value
    x_max = x[indices[0]]
    y_max = y[indices[1]]
    z_max = z[indices[2]]

    # Create the output directory if it doesn't exist
    os.makedirs("output", exist_ok=True)

    if output_all_vtr:
        if x2 is not None and y2 is not None and z2 is not None:
            filename = f"{run_dir}/dose_{x1}_{y1}_{z1}_to_{x2}_{y2}_{z2}"
        else:
            filename = f"{run_dir}/dose_{x1}_{y1}_{z1}"
        pyevtk.hl.gridToVTK(
            filename, x, y, z, cellData={"$\Delta$ Dose ($\mu$Sv/hr)": dose_array}
        )

    return dose_array, x, y, z, plot_bounds


def plot_slice(dose_array, x, y, z, slice_axis, slice_location, plot_bounds):

    # Map axis names to indices
    axis_dict = {"x": 0, "y": 1, "z": 2}

    # Map axis names to arrays
    array_dict = {"x": x, "y": y, "z": z}

    # Find the index of slice_location in the appropriate array
    if slice_axis in axis_dict:
        slice_index = (
            np.abs([x, y, z][axis_dict[slice_axis]] - slice_location)
        ).argmin()
    else:
        raise ValueError(
            f"Invalid slice axis: {slice_axis}. Expected one of: 'x', 'y', 'z'."
        )

    # Create a 2D mask for zero values
    mask = np.ma.masked_where(
        dose_array.take(slice_index, axis=axis_dict[slice_axis]) == 0,
        dose_array.take(slice_index, axis=axis_dict[slice_axis]),
    )

    transposed_mask = np.transpose(mask)

    fig, ax = plt.subplots(figsize=(8, 8))

    # Determine the other two axes
    other_axes = [axis for axis in ["x", "y", "z"] if axis != slice_axis]

    # Calculate the minimum and maximum values of the data
    min_value = np.min(transposed_mask)
    max_value = np.max(transposed_mask)
    # Ensure min_value is not zero to avoid division by zero error
    min_value = max(min_value, 1e-30)

    # Create 2D grids for the the sets of coordinates
    coord1, coord2 = np.meshgrid(array_dict[other_axes[0]], array_dict[other_axes[1]])

    # Ensure the shapes of coords1 and coords2 and transposed_mask match
    coord1 = coord1[: transposed_mask.shape[0], : transposed_mask.shape[1]]
    coord2 = coord2[: transposed_mask.shape[0], : transposed_mask.shape[1]]

    levels_plot = np.logspace(np.log10(min_value), np.log10(max_value), 12)

    # Plot the data with the mask applied - log plot
    im = ax.contourf(
        coord1,
        coord2,
        transposed_mask,
        levels=levels_plot,
        cmap="rainbow",
        norm=LogNorm(),
        extent=[x.min(), x.max(), y.min(), y.max()],
    )
    cbar = fig.colorbar(im, ax=ax, label="$\Delta$ Dose ($\mu$Sv/hr)")

    # Set the colorbar ticks to be at the contour levels
    cbar.set_ticks(levels_plot)
    cbar.set_ticklabels(["{:.1e}".format(level) for level in levels_plot])

    # Create 6 logarithmically spaced levels between the minimum and maximum values
    levels_contour = np.logspace(np.log10(min_value), np.log10(max_value), 6)

    contours = ax.contour(
        coord1, coord2, transposed_mask, levels=levels_contour, colors="black"
    )

    def format_func(value, tick_number):
        return f"{value:.1e}"

    ax.clabel(contours, inline=True, fontsize=8, fmt=FuncFormatter(format_func))

    # Handle case when point(s) is not within one of the available ITER STLs i.e no STL plotted
    if len(plot_bounds) == 7:
        # The last entry in the bounds list is the stl file
        stl_to_plot = plot_bounds[6]

        stl_path = os.path.join("src", "f4epurity", "resources", "building_stl_files")

        stl = STL(os.path.join(stl_path, stl_to_plot), format="ascii")

        axis_pairs = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}

        # Loop over each solid in the stl file and plot the slice
        for solid in stl.solids:
            slice = solid.slice(plane=slice_axis, intercept=slice_location)
            pairs = getattr(slice, axis_pairs[slice_axis][0] + "_pairs"), getattr(
                slice, axis_pairs[slice_axis][1] + "_pairs"
            )
            for pair in zip(*pairs):
                ax.plot(*pair, color="black", linewidth=1)

    ax.set_title("Estimated Dose Deviation")
    ax.set_xlabel("Horizontal Coordinate")
    ax.set_ylabel("Vertical Coordinate")
    ax.set_aspect("equal")

    return fig
