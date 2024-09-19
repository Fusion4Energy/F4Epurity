import os
from f4epurity.main import parse_arguments, _validate_source_coordinates_input
import pytest


# path to the test files
# get the path to this file
cp = os.path.dirname(os.path.abspath(__file__))
TESTFILES = os.path.join(cp, "data", "cli")


def test_cli_command(tmpdir):
    # Define the CLI command to test
    flux_path = f"{TESTFILES}/flux.vtu"

    command = [
        "--element",
        "Ta",
        "--delta_impurity",
        "0.1",  # weight percentage
        "--input_flux",
        "dummy",
        "--irrad_scenario",
        "SA2",
        "--x1",
        "-835",
        "--y1",
        "1994",
        "--z1",
        "1230",
        "--decay_time",
        "1e6",  # seconds
        "--root_output",
        str(tmpdir),
        "--input_flux",
        str(flux_path),
    ]

    # Run the command
    args = parse_arguments(command)
    assert args.x1 == [-835]


@pytest.mark.parametrize("format", ["json", "yaml"])
def test_cli_cfg(tmpdir, format):
    # Define the CLI command to test
    config_path = f"{TESTFILES}/cfg_test.{format}"
    flux_path = f"{TESTFILES}/flux.vtu"
    command = [
        "--cfg",
        str(config_path),
        # override the dummy input flux
        "--input_flux",
        str(flux_path),
        "--root_output",
        str(tmpdir),
    ]

    args = parse_arguments(command)
    assert args.input_flux == flux_path  # check override works
    assert args.irrad_scenario == "SA2"  # check a random value


def test_csv_source_parsing():

    base_commands = [
        "--element",
        "Ta",
        "--delta_impurity",
        "0.1",  # weight percentage
        "--input_flux",
        "dummy",
        "--irrad_scenario",
        "SA2",
        "--decay_time",
        "1e6",  # seconds
    ]

    # provide both csv and coordinates should cause an error
    additional_commands = [
        "--x1",
        "-835",
        "--sources_csv",
        "dummy",
    ]

    commands = base_commands + additional_commands
    with pytest.raises(SystemExit):
        parse_arguments(commands)

    # provide only one or two of the three coordinates should cause an error
    additional_commands = [
        "--x1",
        "-835",
        "--y1",
        "511",
    ]

    commands = base_commands + additional_commands
    with pytest.raises(SystemExit):
        parse_arguments(commands)

    # provide only the .csv is ok
    additional_commands = [
        "--sources_csv",
        "dummy",
    ]

    commands = base_commands + additional_commands
    parse_arguments(commands)


def test_validate_source_coordinates_input():
    flux_path = f"{TESTFILES}/flux.vtu"
    base_commands = [
        "--element",
        "Ta",
        "--delta_impurity",
        "0.1",  # weight percentage
        "--input_flux",
        flux_path,
        "--irrad_scenario",
        "SA2",
        "--decay_time",
        "1e6",  # seconds
    ]

    # provide a point
    csv_path = f"{TESTFILES}/sources1.csv"
    additional_commands = ["--sources_csv", csv_path]
    commands = base_commands + additional_commands
    args = parse_arguments(commands)
    args = _validate_source_coordinates_input(args)
    assert len(args.x1) == 2

    # provide lines
    csv_path = f"{TESTFILES}/sources2.csv"
    additional_commands = ["--sources_csv", csv_path]
    commands = base_commands + additional_commands
    args = parse_arguments(commands)
    _validate_source_coordinates_input(args)
    assert len(args.x2) == 2

    # incorrect format
    csv_path = f"{TESTFILES}/sources3.csv"
    additional_commands = ["--sources_csv", csv_path]
    commands = base_commands + additional_commands
    args = parse_arguments(commands)
    with pytest.raises(KeyError):
        _validate_source_coordinates_input(args)
