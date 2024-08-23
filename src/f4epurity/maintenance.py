import sys

import importlib.resources as pkg_resources
import numpy as np
import pandas as pd
import pyvista as pv

from f4epurity.dose import dose_from_line_source

def read_maintenance_locations(workstation, location):

    # Read the excel document containing the workstations
    with pkg_resources.path('f4epurity.resources', 'workstations.xlsx') as file_path:
        df = pd.read_excel(file_path)

    # Convert the 'workstation' column to string
    df['workstation'] = df['workstation'].astype(str)

    # Check if the location exists in the DataFrame
    if location not in df['location'].values:
        print(f"Error: No matching location found for {location}")
        sys.exit(1)

    if workstation.lower() == 'all':
        # If 'all' is passed as the workstation, return a list of coordinates and workstation names for all workstations
        df_filtered = df[df['location'] == str(location)]
        coordinates = df_filtered[['x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']].values.tolist()
        workstations = df_filtered['workstation'].values.tolist()
    else:
        # Filter the DataFrame to include only rows for the given workstation and location
        df_filtered = df[(df['workstation'] == str(workstation)) & (df['location'] == str(location))]

        # If workstation doesnt exist, return an error
        if df_filtered.empty:
            print(f"Error: No matching workstation and location found for {workstation}, {location}")
            sys.exit(1)

        # Get the coordinates and workstation name for the given workstation and location
        coordinates = [df_filtered[['x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']].values[0].tolist()]
        workstations = [workstation]

    return coordinates, workstations

def get_dose_at_workstation(dose, source_start, source_end, min_x, max_x, min_y, max_y, min_z, max_z, is_line_source):
    # Initialize max_dose and max_dose_coord
    max_dose = 0
    max_dose_coord = None

    # Check if the source is within the workstation volume - in which case just return the dose
    #TODO should print which source is in what volume (multiple source case)
    if min_x <= source_start[0] <= max_x and min_y <= source_start[1] <= max_y and min_z <= source_start[2] <= max_z:
        print('The source is within the workstation volume')
        return dose, source_start

    if is_line_source:
        # Line source mode
        # Check if any part of the line source is within the workstation volume
        if (min_x <= source_end[0] <= max_x and min_y <= source_end[1] <= max_y and min_z <= source_end[2] <= max_z) or \
           (source_start[0] <= min_x <= source_end[0] or source_start[0] >= max_x >= source_end[0]) and \
           (source_start[1] <= min_y <= source_end[1] or source_start[1] >= max_y >= source_end[1]) and \
           (source_start[2] <= min_z <= source_end[2] or source_start[2] >= max_z >= source_end[2]):
            print('The line source intersects the workstation volume')
            return dose, source_start
        num_samples_volume = 10

        # Create arrays of x, y, z coordinates within the volume
        volume_x = np.linspace(min_x, max_x, num_samples_volume)
        volume_y = np.linspace(min_y, max_y, num_samples_volume)
        volume_z = np.linspace(min_z, max_z, num_samples_volume)

        # Calculate the dose at each point in the volume from the line source
        doses = np.empty((num_samples_volume, num_samples_volume, num_samples_volume))
        for i in range(num_samples_volume):
            for j in range(num_samples_volume):
                for k in range(num_samples_volume):
                    dose_at_point = dose_from_line_source(
                                    dose, 
                                    source_start[0], 
                                    source_start[1], 
                                    source_start[2], 
                                    source_end[0], 
                                    source_end[1], 
                                    source_end[2], 
                                    volume_x[i], 
                                    volume_y[j], 
                                    volume_z[k]
                                )
                    doses[i, j, k] = dose_at_point
                    # Update max_dose and max_dose_coord if the dose at this point is greater than max_dose
                    if dose_at_point > max_dose:
                        max_dose = dose_at_point
                        max_dose_coord = (volume_x[i], volume_y[j], volume_z[k])

        # Return the maximum dose in the volume and its coordinates
        return max_dose, max_dose_coord
    else:
        # Point source
        # Calculate the coordinates of the point within the workstation volume that is closest to the source
        closest_x = min(max(min_x, source_start[0]), max_x)
        closest_y = min(max(min_y, source_start[1]), max_y)
        closest_z = min(max(min_z, source_start[2]), max_z)

        # Calculate the distance from the source to the closest point
        distance = np.sqrt((closest_x - source_start[0])**2 + (closest_y - source_start[1])**2 + (closest_z - source_start[2])**2)

        # Calculate the dose at the closest point using 1/r^2
        max_dose = dose / (4 * np.pi * distance**2)

        return max_dose, (closest_x, closest_y, closest_z)
    
def dose_within_workstation(grid_file, box_bounds):
    # Uses pyvista to determine the max value in a given volume on a grid
    # Read the total dose vtr file 
    grid = pv.read(grid_file)
    box = pv.Box(box_bounds)
    # Using the workstation bounds, find the max dose value
    clipped = grid.clip_box(box.bounds, invert=False)
    max_dose = clipped.cell_data['Dose_Total'].max()

    return max_dose
