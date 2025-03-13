import difflib

import actigamma as ag
import matplotlib.pyplot as plt

# global_counter = 0
import numpy as np
import pytest

# from f4epurity.main import calculate_total_activity  # Import the function from main.py
# from f4epurity.mcnp_source_calc import (
#     convert_to_actigamma_format,  # Import the function from mcnp_source_calc.py
# )
from f4epurity.psource import (
    GlobalPointSource,
    PointSource,
    _convert_elem_name,
    insert_wrapped_values,
    insert_wrapped_values_2,
)

DB = ag.Decay2012Database()
GRID = ag.EnergyGrid(bounds=ag.linspace(0, 14e6, 10000))  # TODO
LC = ag.MultiTypeLineAggregator(DB, GRID)


class TestPointSource:
    def test_actigamma_import(self):
        try:
            db = ag.Decay2012Database()
            assert db is not None
        except ImportError:
            assert False, "Failed to import actigamma"

    def test_compute_lines(self):
        activities = {
            "Co60": [np.array([100])],
        }
        x1 = 0  # , -1107.091, -953.078
        y1 = 0  # , 1720.756, 1537.211
        z1 = 0  # , 1241.0, 1378.307
        coord = [x1, y1, z1]
        mass = 1.0
        psource = PointSource(activities, coord, mass)
        X_filtered, Y_filtered, t = psource._compute_lines()

        y = np.array(Y_filtered)
        x = np.array(X_filtered)
        newy = y[y > 1e-30]
        newx = x[y > 1e-30]
        xmin = newx.min()
        xmax = newx.max()

        # assert pytest.approx(xmin, 1e-2) == 1.1732e6
        # assert pytest.approx(xmax, 1e-2) == 1.3325e6
        # assert pytest.approx(newy.max(), 1e-2) == newy.min()
        X2_filtered, Y2_filtered, t = psource._compute_lines()
        assert X2_filtered == X_filtered
        assert Y2_filtered == Y_filtered
        assert y.sum() == 199.86007030599998

    # def test_compute_lines_multi_unstable(self):
    #     activities = {
    #         "Co60": [np.array([1])],
    #         "Co60m": [np.array([0.5])],
    #     }
    #     x1 = -277.391  # , -1107.091, -953.078
    #     y1 = 2086.437  # , 1720.756, 1537.211
    #     z1 = 1241.0  # , 1241.0, 1378.307
    #     coord = [x1, y1, z1]
    #     mass = 1.0
    #     psource = PointSource(activities, coord, mass)
    #     X_filtered, Y_filtered = psource._compute_lines()


class TestGlobalSource:
    def test_to_sdef(self, tmpdir):
        delta_impurity = 0.05
        points = [
            {
                "activities": {
                    "Co60m": [np.array([1])],
                    "Co60": [np.array([1])],
                },
                "coord": [0, 0, 0],
                "mass": 2.0,
            },
            {
                "activities": {
                    "Co60m": [np.array([0.5])],
                    "Co60": [np.array([0.5])],
                },
                "coord": [0, 0, 1],
                "mass": 10.0,
            },
        ]

        # Create PointSource instances for each point
        point_sources = [
            PointSource(point["activities"], point["coord"], point["mass"])
            for point in points
        ]

        global_source = GlobalPointSource([point_sources[0], point_sources[1]])
        global_source.to_sdef(tmpdir.join("mcnp_source.txt"))

        if not tmpdir.join("mcnp_source.txt").exists():
            assert False, "Output file not created"
        else:
            with open(tmpdir.join("mcnp_source.txt")) as f:
                lines = f.readlines()


def test_insert_wrapped_values():
    # Test with a simple example
    init_str = "SI L"
    values = np.array([100, 200000, 300000000, 40, 0.5, 6000, 7000, 800, 90, 10])
    wrapped_values = insert_wrapped_values(init_str, values, 60)
    assert (
        wrapped_values
        == "SI L\n       1.00e+02 2.00e+05 3.00e+08 4.00e+01 5.00e-01 6.00e+03\n       7.00e+03 8.00e+02 9.00e+01 1.00e+01"
    )
    wrapped_values2 = insert_wrapped_values_2(init_str, values, 60)
    assert (
        wrapped_values2
        == "SI L 100.000 200000.000 300000000.000 40.000 0.500 6000.000\n       7000.000 800.000 90.000 10.000"
    )


@pytest.mark.parametrize(
    "inp, expected",
    [
        ("Co060m", "Co60m"),
        ("Co060", "Co60"),
        ("H001", "H1"),
        ("Xe162", "Xe162"),
    ],
)
def test_convert_elem_name(inp: str, expected: str):
    assert _convert_elem_name(inp) == expected
