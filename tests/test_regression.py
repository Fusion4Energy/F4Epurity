import os

from f4epurity.main import process_sources
from f4epurity.parsing import parse_arguments


# path to the test files
# get the path to this file
cp = os.path.dirname(os.path.abspath(__file__))
TESTFILES = os.path.join(cp, "data", "regression")


def test_main(tmpdir):

    flux_path = f"{TESTFILES}/11-L3-03N.vtu"
    sources_path = f"{TESTFILES}/11-L3-03N_valves_clipped.csv"

    command = [
        "--element",
        "Ta",
        "--delta_impurity",
        "0.1",  # weight percentage
        "--input_flux",
        str(flux_path),
        "--irrad_scenario",
        "SA2",
        "--decay_time",
        "1e6",  # seconds
        "--root_output",
        str(tmpdir),
        "--input_flux",
        str(flux_path),
        "--sources_csv",
        str(sources_path),
    ]

    # Run the command
    args = parse_arguments(command)
    process_sources(args)

    # Check the output files have been produced
    for file in os.listdir(tmpdir):
        if "F4Epurity" in file:
            assert os.path.exists(os.path.join(tmpdir, file, "dose_total.vtr"))

    # TODO add more substantial checks for regression
