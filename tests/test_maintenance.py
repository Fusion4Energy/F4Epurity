import pytest
import pandas as pd
from unittest.mock import patch
from f4epurity.maintenance import read_maintenance_locations, get_dose_at_workstation


@pytest.fixture
def mock_df():
    return pd.DataFrame(
        {
            "workstation": ["1", "2", "3"],
            "location": ["loc1", "loc1", "loc1"],
            "x_min": [0, 0, 0],
            "x_max": [1, 1, 1],
            "y_min": [0, 0, 0],
            "y_max": [1, 1, 1],
            "z_min": [0, 0, 0],
            "z_max": [1, 1, 1],
        }
    )


@patch("pandas.read_excel")
def test_read_maintenance_locations(mock_read_excel, mock_df):
    mock_read_excel.return_value = mock_df

    # Test with valid workstation and location specified
    coordinates, workstations = read_maintenance_locations("1", "loc1")
    assert coordinates == [[0, 1, 0, 1, 0, 1]]
    assert workstations == ["1"]

    # Test with 'all' set for the workstation
    coordinates, workstations = read_maintenance_locations("all", "loc1")
    assert coordinates == [[0, 1, 0, 1, 0, 1], [0, 1, 0, 1, 0, 1], [0, 1, 0, 1, 0, 1]]
    assert workstations == ["1", "2", "3"]

    # Test with invalid workstation and location
    with pytest.raises(SystemExit):
        read_maintenance_locations("invalid", "invalid")


def test_get_dose_at_workstation():
    # Test with source within workstation volume
    dose, coord = get_dose_at_workstation(
        10, [0.5, 0.5, 0.5], None, 0, 1, 0, 1, 0, 1, False
    )
    assert dose == 10
    assert coord == [0.5, 0.5, 0.5]

    # Test with source outside workstation volume
    dose, coord = get_dose_at_workstation(
        10, [1.5, 1.5, 1.5], None, 0, 1, 0, 1, 0, 1, False
    )
    assert dose == pytest.approx(1.06, rel=0.01)
    assert coord == (1, 1, 1)
