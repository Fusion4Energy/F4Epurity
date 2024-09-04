"""
Author: Ross Worrall, UKAEA
Creation Date: 14.08.2024
Description: Module for reading stl files and plotting slices along the x,y and z axes. 
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def vector_plane_intersection(p0, p1, plane_coordinate, plane_normal, tol=1e-4):
    """
    p0, p1 are the two vertices of each edge of the facet
    plane_coordinate is a single x,y,z point on the slice plane
    plane_normal is the normal of the slice plane
    tol is the tolerance
    """

    # set vector to origin
    vec = p1 - p0

    # calculate dot product of plane formal to line vector
    dot_plane_norm_vec = np.dot(plane_normal, vec)

    if abs(dot_plane_norm_vec) > tol:
        # case where intercept found
        magnitude = -np.dot(plane_normal, p0 - plane_coordinate) / dot_plane_norm_vec
        vec = vec * magnitude
        return p0 + vec
    else:
        # case that line does not intercept plane
        return None


class STLFacet:
    normal: float
    vertex_1: float
    vertex_2: float
    vertex_3: float

    def __init__(self, normal, vertex_1, vertex_2, vertex_3):
        self.normal = normal
        self._vertex_1 = np.array([float(v) for v in vertex_1])
        self._vertex_2 = np.array([float(v) for v in vertex_2])
        self._vertex_3 = np.array([float(v) for v in vertex_3])

    @property
    def vertices(self):
        return (self._vertex_1, self._vertex_2, self._vertex_3)

    @property
    def x_min(self):
        return min(self._vertex_1[0], self._vertex_2[0], self._vertex_3[0])

    @property
    def x_max(self):
        return max(self._vertex_1[0], self._vertex_2[0], self._vertex_3[0])

    @property
    def y_min(self):
        return min(self._vertex_1[1], self._vertex_2[1], self._vertex_3[1])

    @property
    def y_max(self):
        return max(self._vertex_1[1], self._vertex_2[1], self._vertex_3[1])

    @property
    def z_min(self):
        return min(self._vertex_1[2], self._vertex_2[2], self._vertex_3[2])

    @property
    def z_max(self):
        return max(self._vertex_1[2], self._vertex_2[2], self._vertex_3[2])

    def coordinate_in_bounding_box(self, coordinate, tol=1e-6):
        if coordinate is None:
            return None
        else:
            return (
                (self.x_min - tol <= coordinate[0] <= self.x_max + tol)
                & (self.y_min - tol <= coordinate[1] <= self.y_max + tol)
                & (self.z_min - tol <= coordinate[2] <= self.z_max + tol)
            )


class STLSolid:
    def __init__(self, facets, name=None):
        self._name = name
        self._facets = facets
        self._vertices_loop = None

    @property
    def name(self):
        return self._name

    @property
    def facets(self):
        return self._facets

    @property
    def vertices_in_loop(self):
        if self._vertices_loop:
            # carry out method to find loop
            all_points = []

    def slice(self, plane="x", intercept=0.0, create_loops=False):

        xs = []
        ys = []
        zs = []

        if plane.lower() == "x":
            vec_intercept = (intercept, 0, 0)
            vec_projects = (1, 0, 0)
        elif plane.lower() == "y":
            vec_intercept = (0, intercept, 0)
            vec_projects = (0, 1, 0)
        elif plane.lower() == "z":
            vec_intercept = (0, 0, intercept)
            vec_projects = (0, 0, 1)
        else:
            msg = f"plane of slice should be either x, y, or z not {plane}"
            raise Exception(msg)

        for i, facet in enumerate(self.facets, start=1):
            # line 1 - v1 v2
            inter_1 = vector_plane_intersection(
                facet.vertices[0], facet.vertices[1], vec_intercept, vec_projects
            )
            inter_1 = inter_1 if facet.coordinate_in_bounding_box(inter_1) else None

            # line 2 - v1 v3
            inter_2 = vector_plane_intersection(
                facet.vertices[0], facet.vertices[2], vec_intercept, vec_projects
            )
            inter_2 = inter_2 if facet.coordinate_in_bounding_box(inter_2) else None

            # line 3 - v2 v3
            inter_3 = vector_plane_intersection(
                facet.vertices[1], facet.vertices[2], vec_intercept, vec_projects
            )
            inter_3 = inter_3 if facet.coordinate_in_bounding_box(inter_3) else None

            # calculate intersections between plotting plane and three endges of facet
            inters = [
                inter for inter in [inter_1, inter_2, inter_3] if inter is not None
            ]

            # check for case that there are 2 intercepts
            if len(inters) == 2:
                if plane.lower() == "x":
                    ys.append([inters[0][1], inters[1][1]])
                    zs.append([inters[0][2], inters[1][2]])
                elif plane.lower() == "y":
                    xs.append([inters[0][0], inters[1][0]])
                    zs.append([inters[0][2], inters[1][2]])
                elif plane.lower() == "z":
                    xs.append([inters[0][0], inters[1][0]])
                    ys.append([inters[0][1], inters[1][1]])
                else:
                    msg = f"plane of slice should be either x, y, or z not {plane}"
                    raise Exception(msg)

        if plane.lower() == "x":
            return STLSlicePX(ys, zs, self._name, create_loops=create_loops)
        elif plane.lower() == "y":
            return STLSlicePY(xs, zs, self._name, create_loops=create_loops)
        elif plane.lower() == "z":
            return STLSlicePZ(xs, ys, self._name, create_loops=create_loops)


class STL:
    def __init__(self, file_path, format="ascii"):
        self._format = format
        self._file_path = file_path
        self._solids = []

        PARSE_METHODS = {
            "ascii": self._parse_stl_ascii,
            "binary": self._parse_stl_binary,
        }

        if format in PARSE_METHODS:
            PARSE_METHODS[format](file_path)
        else:
            msg = f"format provided should be ascii or binary, not {format}"
            raise Exception(msg)

    @property
    def solids(self):
        return self._solids

    def _parse_stl_binary(self, file_path):
        with open(file_path, "rb") as input_file:
            header = input_file.read(80).decode("ascii")
            solid_name = header.split("<")[1]
            n_facets = np.fromfile(input_file, count=1, dtype=np.uint32)[0]
            facets = []

            for i in range(n_facets):
                facet_normal = np.fromfile(input_file, count=3, dtype=np.float32)
                vertex_1 = np.fromfile(input_file, count=3, dtype=np.float32)
                vertex_2 = np.fromfile(input_file, count=3, dtype=np.float32)
                vertex_3 = np.fromfile(input_file, count=3, dtype=np.float32)
                input_file.read(2)

                facet = STLFacet(facet_normal, vertex_1, vertex_2, vertex_3)
                facets.append(facet)

            solid = STLSolid(facets=facets, name=solid_name)
            self._solids.append(solid)

    def _parse_stl_ascii(self, file_path):
        # read stl file
        with open(file_path, "r") as input_file:
            rows = [r for r in input_file]

            solid_indices_starts = [
                i for i, row in enumerate(rows) if row[:5] == "solid"
            ]
            solid_indices_ends = [
                i for i, row in enumerate(rows) if row[:8] == "endsolid"
            ]

            # iterate over solids in file
            for i_start, i_end in zip(solid_indices_starts, solid_indices_ends):
                # split rows into solid sections
                section_solid = rows[i_start : i_end + 1]

                if len(section_solid[1:-1]) % 7 == 0:
                    sections_facets = [
                        section_solid[1:-1][i : i + 7]
                        for i in range(0, len(section_solid[1:-1]), 7)
                    ]
                else:
                    msg = f"number of lines in facet subsection not equal to 7"
                    raise Exception(msg)

                facets = []
                for section_facet in sections_facets:
                    # extract normal and vertices from file
                    normal = np.array(
                        [float(val) for val in section_facet[0].split()[-3:]]
                    )
                    vertex_1 = np.array(
                        [float(val) for val in section_facet[2][:-1].split()[-3:]]
                    )
                    vertex_2 = np.array(
                        [float(val) for val in section_facet[3][:-1].split()[-3:]]
                    )
                    vertex_3 = np.array(
                        [float(val) for val in section_facet[4][:-1].split()[-3:]]
                    )

                    facet = STLFacet(normal, vertex_1, vertex_2, vertex_3)
                    facets.append(facet)

                solid_name = rows[i_start][7:-2].split("<")[0]
                solid = STLSolid(facets, name=solid_name)
                self._solids.append(solid)

    @property
    def facets(self):
        return self._facets

    def plot(self):
        fig = plt.figure()
        ax = plt.axes(projection="3d")

        p1 = (0, 0, 0)
        p2 = (0, 1, 0)
        p3 = (1, 0, 0)

        x_min, x_max = 0.0, 0.0
        y_min, y_max = 0.0, 0.0
        z_min, z_max = 0.0, 0.0

        # ax.add_collection3d(Poly3DCollection((p1, p2, p3)))

        for i, facet in enumerate(self._facets):
            for vertex in facet.vertices:
                if vertex[0] < x_min:
                    x_min = vertex[0]
                if vertex[0] > x_max:
                    x_max = vertex[0]
                if vertex[1] < y_min:
                    y_min = vertex[0]
                if vertex[1] > y_max:
                    y_max = vertex[0]
                if vertex[2] < z_min:
                    z_min = vertex[0]
                if vertex[2] > z_max:
                    z_max = vertex[0]

            ax.add_collection3d(Poly3DCollection(facet.vertices))

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        plt.show()


def check_2d_points_are_equal(p0, p1, tol=1e-4):
    return (p0[0] - tol <= p1[0] <= p0[0] + tol) & (p0[1] - tol <= p1[1] <= p0[1] + tol)


class STLSlice:
    def __init__(self):
        pass


class STLSliceLoopPX:
    def __init__(self, ys, zs, solid_name):
        self._solid_name = solid_name
        self._ys = ys
        self._zs = zs

    @property
    def solid_name(self):
        return self._solid_name

    @property
    def ys(self):
        return self._ys

    @property
    def zs(self):
        return self._zs


class STLSliceLoopPY:
    def __init__(self, xs, zs, solid_name):
        self._solid_name = solid_name
        self._xs = xs
        self._zs = zs

    @property
    def solid_name(self):
        return self._solid_name

    @property
    def xs(self):
        return self._xs

    @property
    def zs(self):
        return self._zs


class STLSliceLoopPZ:
    def __init__(self, xs, ys, solid_name):
        self._solid_name = solid_name
        self._xs = xs
        self._ys = ys

    @property
    def solid_name(self):
        return self._solid_name

    @property
    def xs(self):
        return self._xs

    @property
    def ys(self):
        return self._ys


class STLSlicePX(STLSlice):
    def __init__(self, y_pairs, z_pairs, solid_name, create_loops=False):
        # object parameters
        self._solid_name = solid_name

        # store vertices of intesection
        self._ys = []
        self._zs = []

        self._y_pairs = [y_pair for y_pair in y_pairs]
        self._z_pairs = [z_pair for z_pair in z_pairs]

        self._loops = []

        if len(y_pairs) > 0 and create_loops:
            self._sort_intercepts_into_loop(y_pairs, z_pairs)

    @property
    def y_pairs(self):
        return self._y_pairs

    @property
    def z_pairs(self):
        return self._z_pairs

    @property
    def loops(self):
        return self._loops


class STLSlicePY(STLSlice):
    def __init__(self, x_pairs, z_pairs, solid_name, create_loops=False):
        # object parameters
        self._solid_name = solid_name

        # store vertices of intesection
        self._xs = []
        self._zs = []

        self._x_pairs = [x_pair for x_pair in x_pairs]
        self._z_pairs = [z_pair for z_pair in z_pairs]

        self._loops = []

        if len(x_pairs) > 0 and create_loops:
            self._sort_intercepts_into_loop(x_pairs, z_pairs)

    @property
    def x_pairs(self):
        return self._x_pairs

    @property
    def z_pairs(self):
        return self._z_pairs

    @property
    def loops(self):
        return self._loops


class STLSlicePZ(STLSlice):
    def __init__(self, x_pairs, y_pairs, solid_name, create_loops=False):
        # object parameters
        self._solid_name = solid_name

        # store vertices of intesection
        self._xs = []
        self._ys = []

        self._x_pairs = [x_pair for x_pair in x_pairs]
        self._y_pairs = [y_pair for y_pair in y_pairs]

        self._loops = []

        if len(x_pairs) > 0 and create_loops:
            self._sort_intercepts_into_loop(x_pairs, y_pairs)

    @property
    def x_pairs(self):
        return self._x_pairs

    @property
    def y_pairs(self):
        return self._y_pairs

    @property
    def loops(self):
        return self._loops

    def _sort_intercepts_into_loop(self, x_pairs, y_pairs):
        # set starting points
        p_start = (x_pairs[0][0], y_pairs[0][0])
        p_next = (x_pairs[0][1], y_pairs[0][1])
        x_pairs = x_pairs[1:]
        y_pairs = y_pairs[1:]

        # create stores for points in currnet loop
        xs = [p_start[0]]
        ys = [p_start[1]]

        count = 0

        while len(x_pairs) > 0:
            # loop through all other pairs of points to find matching within tolerance
            for i, (x_pair, y_pair) in enumerate(zip(x_pairs, y_pairs)):
                p0 = (x_pair[0], y_pair[0])
                p1 = (x_pair[1], y_pair[1])

                # check if either of the two points in the next vertex match
                if check_2d_points_are_equal(p_next, p0):
                    xs.append(p0[0])
                    ys.append(p0[1])
                    p_next = p1
                    x_pairs.remove(x_pair)
                    y_pairs.remove(y_pair)

                    # case for loop being closed
                    if check_2d_points_are_equal(p_next, p_start):
                        xs.append(p_next[0])
                        ys.append(p_next[1])

                        loop = STLSliceLoopPZ(xs, ys, self._solid_name)
                        self._loops.append(loop)

                        # check for remaining points and set up new loop
                        if len(x_pairs) > 0:
                            # set starting points
                            p_start = (x_pairs[0][0], y_pairs[0][0])
                            p_next = (x_pairs[0][1], y_pairs[0][1])

                            x_pairs = x_pairs[1:]
                            y_pairs = y_pairs[1:]

                            # create stores for points in currnet loop
                            xs = [p_start[0]]
                            ys = [p_start[1]]
                        break

                elif check_2d_points_are_equal(p_next, p1):
                    xs.append(p1[0])
                    ys.append(p1[1])
                    p_next = p0
                    x_pairs.remove(x_pair)
                    y_pairs.remove(y_pair)

                    # case for loop being closed
                    if check_2d_points_are_equal(p_next, p_start):
                        xs.append(p_next[0])
                        ys.append(p_next[1])

                        loop = STLSliceLoopPZ(xs, ys, self._solid_name)
                        self._loops.append(loop)

                        # check for remaining points and set up new loop
                        if len(x_pairs) > 0:
                            # set starting points
                            p_start = (x_pairs[0][0], y_pairs[0][0])
                            p_next = (x_pairs[0][1], y_pairs[0][1])

                            x_pairs = x_pairs[1:]
                            y_pairs = y_pairs[1:]

                            # create stores for points in currnet loop
                            xs = [p_start[0]]
                            ys = [p_start[1]]
                        break

            # # error handling case
            # # end of loop reached and no finished loop found
            # print('x_pairs')
            # print(x_pairs)
            # print('y_pairs')
            # print(y_pairs)

            # msg = f'error end of points and no end to loop found.'
            # raise Exception(msg)

        # once all pairs of points joined, loop back to starting point
        self._xs.append(p_start[0])
        self._ys.append(p_start[1])
