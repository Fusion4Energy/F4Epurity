from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import actigamma as ag
import numpy as np

DB = ag.Decay2012Database()
GRID = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))  # TODO
LC = ag.MultiTypeLineAggregator(DB, GRID)


class GlobalPointSource:
    def __init__(self, sources: list[PointSource]) -> None:
        self.sources = sources

    def to_sdef(self, outfile: str | Path) -> None:
        sdef_line = "sdef\n"
        si1 = "S11 L "
        sp1 = "SP1 "
        ds2 = "DS2 "
        lines_distributions = []
        for source in self.sources:
            si1 = si1 + f" {source.coord[0]} {source.coord[0]} {source.coord[0]}"
            sp1 = sp1 + f"{source.mass} 0.0 0.0"


class SingleSource(ABC):
    def __init__(
        self,
        activities: dict[str, list[np.ndarray]],
        coord: list | np.ndarray,
        mass: float | None = None,
    ) -> None:
        self.activities = activities
        self.coord = coord
        if mass is None:
            mass = 1.0
        self.mass = mass

    def _compute_lines(self) -> tuple[list, list]:
        inventory = self._compute_inventory()
        hist, bin_edges = LC(inventory)
        X, Y = ag.getplotvalues(bin_edges, hist)
        return X, Y

    @abstractmethod
    def _compute_inventory(self) -> ag.UnstablesInventory:
        pass


class PointSource(SingleSource):
    def _compute_inventory(self) -> ag.UnstablesInventory:
        data = []
        for key, value in self.activities.items():
            act = value[0][0]
            data.append((DB.getzai(key), act))
        inventory = ag.UnstablesInventory(data)
        return inventory


class LineSource(SingleSource):
    pass
