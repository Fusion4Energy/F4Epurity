import csv
import os

import actigamma as ag
import numpy as np

# Global counter for the number of times mcnp_main is called
mcnp_main_call_count = 1
current_path = os.getcwd()


def convert_to_actigamma_format(data):
    # Setup the DB
    db = ag.Decay2012Database()

    inv_data = []

    for point_data in data:
        activities = point_data["activities"]
        for nuclide, activity_str in activities.items():
            # Strip brackets and other characters, then split by spaces
            activity_values = (
                activity_str.replace("array(", "")
                .replace(")", "")
                .replace("[", "")
                .replace("]", "")
                .split()
            )
            for activity in activity_values:
                activity = float(activity)
                if activity > 0:
                    inv_data.append((db.getzai(nuclide.strip()), activity))

    # Format the data for actigamma
    formatted_data = "inv = ag.UnstablesInventory(data=[\n"
    formatted_data += ",\n".join(
        f'    (db.getzai("{nuclide}"), {activity})' for nuclide, activity in inv_data
    )
    formatted_data += "\n])"

    return formatted_data


def dump_mcnp_source(activities: str, x1, y1, z1, run_dir):
    global mcnp_main_call_count

    # Parse all_activities_str to extract coordinates and activities for multiple points
    points = all_activities_str.strip().split("\n\n")
    data = []

    for point_index, point in enumerate(points):
        lines = point.strip().split("\n")
        coordinates_line = lines[0]
        activity_lines = lines[1:]

        # Extract coordinates
        coordinates = {}
        for part in coordinates_line.replace("Coordinates: ", "").split(", "):
            key, value = part.split("=")
            coordinates[key.strip()] = float(value.strip())
        # Extract activities
        activities = {}
        for line in activity_lines:
            nuclide, activity = line.split(": ")
            activities[nuclide.strip()] = activity.strip()

        # Store the parsed data with a unique identifier for each point
        data.append(
            {
                "point": mcnp_main_call_count,
                "coordinates": coordinates,
                "activities": activities,
            }
        )

        coordinates_source_path = os.path.join(os.getcwd(), "coordinates.txt")
        all_coordinates = f"{coordinates['x1']} {coordinates['y1']} {coordinates['z1']}"

        # Second line (coordinates)
        if mcnp_main_call_count == 1:
            second_line = "SI1 L " + all_coordinates + " "
            with open(coordinates_source_path, "w", encoding="utf-8") as f:
                f.write(second_line)
                f.truncate()
        else:
            second_line_update = all_coordinates + " "
            with open(coordinates_source_path, "a", encoding="utf-8") as f:
                f.write(second_line_update)
                f.truncate()

    # Printed parsed data to double check the order is correct
    for i, point_data in enumerate(data):
        print(f"\nPoint {point_data['point']}:")
        print("Coordinates:")
        for key, value in point_data["coordinates"].items():
            print(f"  {key}: {value}")
        print("Activities:")
        for nuclide, activity in point_data["activities"].items():
            print(f"  {nuclide}: {activity}")

    # Convert to actigamma format
    formatted_data = convert_to_actigamma_format(data)
    print("\nData formatted for actigamma")
    # print(formatted_data)

    # Normally gamma but could be something else - "alpha", "beta" if data exists!
    SPECTYPES = ["gamma"]

    # Setup the DB and define an energy grid
    db = ag.Decay2012Database()
    grid = ag.EnergyGrid(bounds=np.linspace(0.0, 4e6, 10000))

    # Ensure the mcnp_source.txt file exists
    mcnp_source_path = os.path.join(current_path, "mcnp_source.txt")
    if not os.path.exists(mcnp_source_path):
        with open(mcnp_source_path, "w", encoding="utf-8") as f:
            f.write("")

    # Read the existing content of the file
    with open(mcnp_source_path, "r", encoding="utf-8") as f:
        existing_content = f.read()

    # Check if the initial line is already present
    if "SDEF PAR=2 POS=d1 ERG FPOS d2" not in existing_content:
        # Otherwise generate the first line
        new_content = "SDEF PAR=2 POS=d1 ERG FPOS d2\n \n"

        # Line 3 Find the CSV file to extract the rel strengths of the sources
        csv_files = [file for file in os.listdir() if file.endswith(".csv")]
        if csv_files:
            source_csv_path = csv_files[0]
            with open(source_csv_path, "r", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                m_values = [float(row["m"]) for row in reader if "m" in row]

                if m_values:
                    total_m = sum(m_values)
                    relative_strengths = [m / total_m for m in m_values]
                    new_content += (
                        "SP1 "
                        + " ".join(f"{strength:.2f}" for strength in relative_strengths)
                        + " $ rel strengths of sources\n"
                    )

        # Generate 4th line (DS2 S line)
        ds2_s_line = "DS2 S " + " ".join(
            str(i) for i in range(3, mcnp_main_call_count + 3)
        )
        # Check if the last value is an integer before adding the comment
        last_value = ds2_s_line.split()[-1]
        if last_value.isdigit():
            next_value = int(last_value) + 1
            while next_value <= mcnp_main_call_count + 3:
                ds2_s_line += f" {next_value}"
                next_value += 1
            ds2_s_line += f" {next_value} $ energy distributions\n"
        else:
            ds2_s_line += "\n"

        new_content += ds2_s_line

        # Write the new content followed by the existing content back to the file
        with open(mcnp_source_path, "w", encoding="utf-8") as f:
            f.write(new_content + existing_content)

    # Generate SIx L and SPx lines for each source point
    with open(mcnp_source_path, "a", encoding="utf-8") as f:
        for i, point_data in enumerate(data):
            activities = point_data["activities"]
            inv_data = []
            for nuclide, activity_str in activities.items():
                activity_values = (
                    activity_str.replace("array(", "")
                    .replace(")", "")
                    .replace("[", "")
                    .replace("]", "")
                    .split()
                )
                for activity in activity_values:
                    activity = float(activity)
                    if activity > 0:
                        inv_data.append((db.getzai(nuclide.strip()), activity))

            inv = ag.UnstablesInventory(data=inv_data)
            lc = ag.MultiTypeLineAggregator(db, grid)
            hist, bin_edges = lc(inv)
            X, Y = ag.getplotvalues(bin_edges, hist)

            # Filter out zero count values
            X_filtered = [x for x, y in zip(X, Y) if y != 0]
            Y_filtered = [y for y in Y if y != 0]

            # Write SIx L line
            si_line = (
                f"SI{i + mcnp_main_call_count + 2} L "
                + " ".join(map(str, X_filtered))
                + f" \n"
            )
            f.write(si_line)

            # Write SPx line
            sp_line = (
                f"SP{i + mcnp_main_call_count + 2} "
                + " ".join(map(str, Y_filtered))
                + f" \n"
            )
            f.write(sp_line)
            mcnp_main_call_count += 1

    with open(coordinates_source_path, "r+", encoding="utf-8") as f:
        second_line = f.readline()
        second_line_update = f.readline()
        combined_line = second_line.strip() + " " + second_line_update.strip() + "\n"
    with open(mcnp_source_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lines[1] = combined_line

    with open(mcnp_source_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
