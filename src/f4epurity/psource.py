from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import actigamma as ag
import numpy as np

GRID = ag.EnergyGrid(bounds=ag.linspace(0, 20e6, 100000))  # TODO

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

        sdef_line = "sdef PAR=2 POS=d1 ERG FPOS d2\n"
        si1 = "SI1 L "
        sp1 = "SP1 "
        ds2 = "DS2 S "
        intial_ds = 3
        lines_distributions = []
        start_SI = "SI{} L "
        start_SP = "SP{} "
        t_gamma = 0
        for idx, psource in enumerate(self.sources):
            X, Y, p_activity = psource._compute_lines()
            t_gamma += np.array(Y).sum()
        fmesh_line = "FMESH4:P  GEOM=XYZ $ Modify the fmesh as you wish\n       ORIGIN 0 0 0\n       IMESH 1 IINTS 1\n       JMESH 1 JINTS 1\n       KMESH 1 KINTS 1\n"
        gamma_per_second_line = (
            f"FM4 {t_gamma:.2E} $ To be multiplied by 3.6E-09 to obtain Sv/hr\n"
        )
        conv_coeff_line = "DE4\n       0.01 0.015 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.10 0.15 0.20 0.3 0.4 0.5 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0\nDF4\n       0.0485 0.1254 0.2050 0.2999 0.3381 0.3572 0.3780 0.4066 0.4399 0.5172 0.7523 1.0041 1.5083 1.9958\n       2.4657 2.9082 3.7269 4.4834 7.4896 12.0153 15.9873 19.9191 23.7600"

        for idx, psource in enumerate(self.sources):
            X, Y, p_activity = psource._compute_lines()
            # Adjounr global distributions
            si1 = insert_wrapped_values_2(si1, psource.coord, CHAR_LIMIT)
            p_gamma = np.array(Y).sum()
            r_gamma = p_gamma / t_gamma
            sp1 = insert_wrap_values_3(sp1, r_gamma, CHAR_LIMIT)
            p_source_idx = intial_ds + idx
            ds2 = insert_wrap_values_3(ds2, p_source_idx, CHAR_LIMIT)

            # add specific point source distribution

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
            f.write(fmesh_line)
            f.write(gamma_per_second_line)
            f.write(conv_coeff_line)


def insert_wrap_values_3(
    initial_str: str, values: float | np.ndarray, limit: int
) -> str:
    """
    Insert values into a string with a limit of characters per line (when not iterable)

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
    new_str = initial_str + ""
    line = 0
    value = f"{values:.3f}"
    if len(new_str.split("\n")[-1]) + len(value) > limit:
        new_str += "\n      "
        line += 1
    new_str += f" {value}"
    return new_str


def insert_wrapped_values_2(
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
    new_str = initial_str + ""
    line = 0
    for value in values:
        value = f"{value:.3f}"
        if len(new_str.split("\n")[-1]) + len(value) > limit:
            new_str += "\n      "
            line += 1
        new_str += f" {value}"
    return new_str


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
    mass : float, optional
        The mass of the source in grams of component times delta impurity.
        If not provided, the mass is assumed to be 1.0.

    Returns
    -------
    str
        A string representing the MCNP source definition.
    """

    def __init__(
        self,
        activities: dict[str, list[np.ndarray]],
        coord: list | np.ndarray,
        mass: float = 1,
    ) -> None:
        self.activities = activities
        self.coord = coord
        self.mass = mass

    def _compute_lines(self) -> tuple[list, list, float]:
        inventory, t_activity, db = self._compute_inventory()
        LC = ag.MultiTypeLineAggregator(db, GRID)
        Y, X = LC(inventory)
        # Filter out zero count values & convert to MeV
        X_filtered = [x for x, y in zip(X, Y) if y != 0]
        Y_filtered = [y for y in Y if y != 0]
        X_filtered = np.array(X_filtered) * 1e-6
        return X_filtered.tolist(), Y_filtered, t_activity

    @abstractmethod
    def _compute_inventory(
        self,
    ) -> tuple[ag.UnstablesInventory, float, ag.DefaultDatabase]:
        pass


class PointSource(SingleSource):
    """
    A class representing a point source

    Parameters
    ----------
    SingleSource : _type_
        Compute the inventory of unstable nuclides for the point source.
    """

    def _compute_inventory(
        self,
    ) -> tuple[ag.UnstablesInventory, float, ag.DefaultDatabase]:
        db = ag.Decay2012Database()
        data = []
        t_activity = 0
        for key, value in self.activities.items():
            act = float(value[0][0]) * self.mass
            if act == 0:
                continue
            key = _convert_elem_name(key)
            data.append((db.getzai(key), act))
            t_activity += act
        inventory = ag.UnstablesInventory(data)
        return inventory, t_activity, db


def _convert_elem_name(f4epurity_name: str) -> str:
    """
    Convert element name from a specific format to the expected format.

    Parameters
    ----------
    inp : str
        The input element name to be converted.

    Returns
    -------
    str
        The converted element name.
    """
    # Remove leading zeros from the numeric part of the element name
    import re

    match = re.match(r"([A-Za-z]+)(0*)(\d+)([A-Za-z]*)", f4epurity_name)
    if match:
        element, zeros, number, suffix = match.groups()
        return f"{element}{int(number)}{suffix}"
    return f4epurity_name


class LineSource(SingleSource):
    """
    A class representing a line source

    Parameters
    ----------
    SingleSource : _type_
       Compute the inventory of unstable nuclides for the line source.
    """

    pass
