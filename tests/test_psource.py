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
    insert_wrapped_values,
)

DB = ag.Decay2012Database()
GRID = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))  # TODO
LC = ag.MultiTypeLineAggregator(DB, GRID)


class TestPointSource:
    def test_compute_inventory(self):
        activities = {
            "Ta182": [np.array([0.00010949])],
            "Ta182m": [np.array([0.1])],
            "Ta182n": [np.array([0.0])],
        }
        coord = [277.391, 2086.437, 1241.0]
        psource = PointSource(activities, coord)
        psource._compute_inventory()

    def test_actigamma_import(self):
        try:
            db = ag.Decay2012Database()
            assert db is not None
        except ImportError:
            assert False, "Failed to import actigamma"

    def test_compute_lines(self):
        activities = {
            "Co60": [np.array([1])],
        }
        x1 = -277.391  # , -1107.091, -953.078
        y1 = 2086.437  # , 1720.756, 1537.211
        z1 = 1241.0  # , 1241.0, 1378.307
        coord = [x1, y1, z1]
        mass = 1.0
        psource = PointSource(activities, coord, mass)
        X_filtered, Y_filtered = psource._compute_lines()
        plt.plot(X_filtered, Y_filtered, marker="x")
        plt.semilogx()
        plt.savefig(r"D:\DATA\digiama\Desktop\test\line.png")

        y = np.array(Y_filtered)
        x = np.array(X_filtered)
        newy = y[y > 1e-3]
        newx = x[y > 1e-3]
        xmin = newx.min()
        xmax = newx.max()

        assert pytest.approx(xmin, 1e-2) == 1.1732e6
        assert pytest.approx(xmax, 1e-2) == 1.3325e6
        assert pytest.approx(newy.max(), 1e-2) == newy.min()


class TestGlobalSource:
    def test_to_sdef(self, tmpdir):
        points = [
            {
                "activities": {"Co60m": [np.array([0.00456255])]},
                "coord": [0, 0, 0],
                "mass": 1.0,
            },
            {
                "activities": {"Co60": [np.array([0.00053245])]},
                "coord": [0, 0, 0],
                "mass": 1.0,
            },
            # {
            #     "activities": {"Ta182": [np.array([0.00010949])]},
            #     "coord": [-1107.091, 1720.756, 1241.0],
            #     "mass": 2.0,
            # },
            # {
            #     "activities": {"Ta182": [np.array([0.00011723])]},
            #     "coord": [-953.078, 1537.211, 1378.307],
            #     "mass": 3.0,
            # },
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
