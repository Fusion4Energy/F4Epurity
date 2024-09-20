import numpy as np
import pandas as pd
from jsonargparse import ArgumentParser, ActionConfigFile
from scipy.interpolate import interp1d, interp2d
import os
import json
import datetime
from importlib.resources import files, as_file
import urllib.request

livechart = "https://nds.iaea.org/relnsd/v0/data?"


def isnumber(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def pos_to_char(pos):
    if pos > 0:
        return chr(pos + 108)
    else:
        return ""


# Class to extract the decay data from the IAEA livechart
class Nuclide:
    def __init__(self, name, activity, filter):
        self.name = name
        self.activity = activity
        self.filter = filter
        self.energies = []
        self.intensities = []

    def load(self):
        parts = self.name.split("-")
        if self.name[-1].isdigit():
            nuclide_string = parts[1] + parts[0]
            isomer = ""
        else:
            nuclide_string = parts[1][:-1] + parts[0]
            isomer = self.name[-1]
        url = (
            livechart + "fields=decay_rads&nuclides=" + nuclide_string + "&rad_types=g"
        )
        req = urllib.request.Request(url)
        req.add_header(
            "User-Agent",
            "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0",
        )
        df = pd.read_csv(urllib.request.urlopen(req))
        if ("energy" in df) and ("intensity" in df) and ("p_energy" in df):
            levels = list(sorted(set(df["p_energy"])))
            for energy, intensity, p_energy in zip(
                df.energy, df.intensity, df.p_energy
            ):
                if isnumber(energy) and isnumber(intensity) and isnumber(p_energy):
                    if isomer == pos_to_char(levels.index(p_energy)):
                        energy = float(energy) / 1000  # convert to MeV
                        intensity = self.activity * float(intensity) / 100.0
                        if intensity > self.filter:
                            self.energies.append(energy)
                            self.intensities.append(intensity)

    @property
    def print(self):
        outstring = ""
        for energy, intensity in zip(self.energies, self.intensities):
            outstring += self.name + "\t" + str(energy) + "\t" + str(intensity) + "\n"
        return outstring


def load_table(file_path):
    """Load tabulated data from a CSV file."""
    return pd.read_csv(file_path)


def interpolate_value(table, energy, column):
    """Interpolate the value for the given energy."""
    energies = table["energy"].values
    values = table[column].values
    interpolator = interp1d(energies, values, kind="linear", fill_value="extrapolate")
    return interpolator(energy)


def calculate_intensity(I0, B, mu, x):
    """Calculate the intensity through the shield."""
    return I0 * B * np.exp(-mu * x)


def main():
    # Parse the command line arguments
    parser = ArgumentParser(
        description="Simple estimator of shielding for a photon source"
    )
    parser.add_argument("--cfg", action=ActionConfigFile, help="Path to config file")
    parser.add_argument(
        "--nuclide",
        type=str,
        required=True,
        help="Name of the nuclide (e.g., Co-60, Cs-137)",
    )
    parser.add_argument(
        "--activity",
        type=float,
        required=True,
        help="Activity of the photon source (Bq)",
    )
    parser.add_argument(
        "--shield_thickness",
        type=float,
        required=True,
        help="Thickness of the shield (cm)",
    )
    parser.add_argument(
        "--material",
        type=str,
        required=True,
        help="Material of the shield (concrete, lead, water, iron)",
    )
    parser.add_argument(
        "--root_output",
        type=str,
        default="output",
        help="Root directory for output files",
    )

    args = parser.parse_args()

    # Currently supported materials
    supported_materials = ["lead", "iron", "water", "concrete"]

    # Check that the input material is supported
    if args.material not in supported_materials:
        raise ValueError(f"Material must be one of {supported_materials}.")

    with as_file(
        files("resources.mass_att_coeff").joinpath("density_table.csv")
    ) as density_table_path:
        density_table = load_table(density_table_path)

    with as_file(
        files("resources.mass_att_coeff").joinpath(f"{args.material}_mu.csv")
    ) as mu_table_path:
        mu_table = load_table(mu_table_path)

    # Get the density of the material
    density = density_table.loc[
        density_table["material"] == args.material, "density"
    ].values[0]

    total_intensity = 0

    # Instantiate the Nuclide class with the provided nuclide name
    nuclide = Nuclide(name=args.nuclide, activity=args.activity, filter=0)
    nuclide.load()

    # Iterate over each energy and intensity
    for energy, intensity in zip(nuclide.energies, nuclide.intensities):
        # Interpolate the mass attenuation coefficient for the given energy
        mu_mass = interpolate_value(mu_table, energy, "mu_mass")

        # Calculate the linear attenuation coefficient
        mu = mu_mass * density

        # Calculate the number of mean free paths
        mfp = mu * args.shield_thickness

        # Buildup =1 #TODO: Implement buildup factor interpolation
        B = 1

        # Calculate the intensity for this energy and add to the total intensity
        intensity = calculate_intensity(intensity, B, mu, args.shield_thickness)
        total_intensity += intensity

    # Create a unique directory for this run
    root = args.root_output
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = f"{root}/ShieldingEstimator_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)

    # Save the metadata
    metadata = vars(args)
    metadata_path = os.path.join(run_dir, "metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)

    # Save the result
    result_path = os.path.join(run_dir, "result.txt")
    with open(result_path, "w") as f:
        f.write(f"Nuclide: {args.nuclide}\n")
        f.write(f"Initial Activity: {args.activity}\n")
        f.write(f"Shield Thickness (x): {args.shield_thickness}\n")
        f.write(f"Material: {args.material}\n")
        f.write(f"Density (rho): {density}\n")
        f.write(f"Total Transmitted Intensity: {total_intensity}\n")

    print(f"Results saved to {run_dir}")


if __name__ == "__main__":
    main()
