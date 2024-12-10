import numpy as np
import pytest

from f4epurity.psource import GlobalSource, PointSource


class TestPointSource:
    @pytest.fixture
    def point_source(self):
        source = PointSource(
            {"Ta182": [np.array([1e-3])], "Ta182m": [np.array([1e-5])]},
            [0.0, 0.0, 0.0],
            1.0,
        )
        return source

    def test_compute_lines(self, point_source: PointSource):
        x, y = point_source._compute_lines()
        assert len(x) == len(y)

    def test_compute_inventory(self, point_source: PointSource):
        inventory = point_source._compute_inventory()
        assert inventory.zais == [731820, 731821]


class TestGlobalSource:
    def test_to_sdef(self):
        pass
