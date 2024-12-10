from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np


class GlobalSource:
    def __init__(self, sources: list[SingleSource]) -> None:
        self.sources = sources

    def to_sdef(outfile: str | Path) -> None:
        pass


class SingleSource(ABC):
    def __init__(
        self,
        activities: dict[str, float],
        coord: list | np.ndarray,
        mass: float | None = None,
    ) -> None:
        self.activities = activities
        self.coord = coord
        if mass is None:
            mass = 1.0
        self.mass = mass

    @abstractmethod
    def _compute_lines(self) -> tuple[np.ndarray, np.ndarray]:
        pass


class PointSource(SingleSource):
    def _compute_lines(self) -> tuple[np.ndarray, np.ndarray]:
        pass


class LineSource(SingleSource):
    pass
