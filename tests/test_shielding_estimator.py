import pytest
from unittest.mock import patch, mock_open, MagicMock
import pandas as pd
import numpy as np
from f4epurity.shielding_estimator import (
    isnumber,
    pos_to_char,
    Nuclide,
    load_table,
    interpolate_value,
    calculate_intensity,
    main,
)


def test_isnumber():
    assert isnumber("123") is True
    assert isnumber("abc") is False


def test_pos_to_char():
    assert pos_to_char(1) == "m"
    assert pos_to_char(0) == ""
    assert pos_to_char(-1) == ""


@patch("f4epurity.shielding_estimator.urllib.request.urlopen")
@patch("f4epurity.shielding_estimator.pd.read_csv")
def test_nuclide_load(mock_read_csv, mock_urlopen):
    mock_urlopen.return_value = MagicMock()
    mock_read_csv.return_value = pd.DataFrame(
        {"energy": ["1000", "2000"], "intensity": ["50", "30"], "p_energy": ["0", "0"]}
    )
    nuclide = Nuclide(name="Co-60", activity=1000, filter=0)
    nuclide.load()
    assert len(nuclide.energies) == 2
    assert len(nuclide.intensities) == 2


@patch(
    "builtins.open", new_callable=mock_open, read_data="material,density\nlead,11.34\n"
)
@patch("f4epurity.shielding_estimator.pd.read_csv")
def test_load_table(mock_read_csv, mock_open):
    mock_read_csv.return_value = pd.DataFrame(
        {"material": ["lead"], "density": [11.34]}
    )
    table = load_table("dummy_path")
    assert table.loc[0, "material"] == "lead"
    assert table.loc[0, "density"] == 11.34


def test_interpolate_value():
    table = pd.DataFrame({"energy": [0.1, 0.2, 0.3], "mu_mass": [1.0, 2.0, 3.0]})
    value = interpolate_value(table, 0.15, "mu_mass")
    assert np.isclose(value, 1.5)


def test_calculate_intensity():
    intensity = calculate_intensity(1000, 1, 0.1, 10)
    assert np.isclose(intensity, 367.87944117144233)
