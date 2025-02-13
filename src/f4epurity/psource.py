from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import actigamma as ag
import numpy as np

DB = ag.Decay2012Database()
GRID = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))  # TODO
LC = ag.MultiTypeLineAggregator(DB, GRID)
CHAR_LIMIT = 120


class GlobalPointSource:
    def __init__(self, sources: list[PointSource]) -> None:
        self.sources = sources

    def to_sdef(self, outfile: str | Path) -> None:
        """
        Generate an MCNP SDEF (source definition) card for the global point sources and
        write it to a file named mcnp_source.txt.

        Parameters
        ----------
        outfile : str | Path
            The path to the output file where the SDEF card will be written.

        Returns
        -------
        None
        """
        # compute total mass
        t_mass = 0
        for psource in self.sources:
            t_mass += psource.mass

        sdef_line = "sdef PAR=2 POS=d1 ERG FPOS d2\n"
        si1 = "SI1 L "
        sp1 = "SP1 "
        ds2 = "DS2 S "
        intial_ds = 3
        lines_distributions = []
        start_SI = "SI{} L "
        start_SP = "SP{} "
        for idx, psource in enumerate(self.sources):
            # Adjounr global distributions
            si1 = si1 + f" {psource.coord[0]} {psource.coord[1]} {psource.coord[2]} "
            sp1 = sp1 + f"{psource.mass/t_mass:.3f} "
            p_source_idx = intial_ds + idx
            ds2 = ds2 + f" {p_source_idx} "

            # add specific point source distribution
            X, Y = psource._compute_lines()
            SI_line = start_SI.format(p_source_idx)
            SP_line = start_SP.format(p_source_idx)
            SI_line = insert_wrapped_values(SI_line, X, CHAR_LIMIT)
            SP_line = insert_wrapped_values(SP_line, Y, CHAR_LIMIT)
            lines_distributions.append(SI_line)
            lines_distributions.append(SP_line)

        with open(outfile, "w") as f:
            f.write(sdef_line)
            f.write(si1 + "\n")
            f.write(sp1 + "\n")
            f.write(ds2 + "\n")
            for line in lines_distributions:
                f.write(line + "\n")


def insert_wrapped_values(
    initial_str: str, values: list | np.ndarray, limit: int
) -> str:
    """
    Insert values into a string with a limit of characters per line

    Parameters
    ----------
    initial_str : str
        initial string
    values : list | np.ndarray
        values to be inserted
    limit : int
        character limit per line
    Returns
    -------
    str
        The resulting string with values indented according to the character limit.
    """
    new_str = initial_str + "\n      "
    line = 0
    for value in values:
        value = f"{value:.2e}"
        if len(new_str.split("\n")[-1]) + len(value) > limit:
            new_str += "\n      "
            line += 1
        new_str += f" {value}"
    return new_str


class SingleSource(ABC):
    """
    Generate a single source definition for MCNP simulation.

    Parameters
    ----------
    activities : dict
        A dictionary where keys are nuclide names and values are their activities.
    coord : list[float]
        List of coordinates [x1, y1, z1] for point and also [x2, y2, z2] for line source
    mass : float
        The mass of the source.

    Returns
    -------
    str
        A string representing the MCNP source definition.
    """

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
        # Filter out zero count values
        X_filtered = [x for x, y in zip(X, Y) if y != 0]
        Y_filtered = [y for y in Y if y != 0]
        return X_filtered, Y_filtered

    @abstractmethod
    def _compute_inventory(self) -> ag.UnstablesInventory:
        pass


class PointSource(SingleSource):
    """
    A class representing a point source

    Parameters
    ----------
    SingleSource : _type_
        Compute the inventory of unstable nuclides for the point source.
    """

    def _compute_inventory(self) -> ag.UnstablesInventory:
        data = []
        for key, value in self.activities.items():
            act = float(value[0][0])
            if act == 0:
                continue
            data.append((DB.getzai(key), act))
        inventory = ag.UnstablesInventory(data)
        return inventory


class LineSource(SingleSource):
    """
    A class representing a line source

    Parameters
    ----------
    SingleSource : _type_
       Compute the inventory of unstable nuclides for the line source.
    """

    pass
