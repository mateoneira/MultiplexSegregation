"""
Module: processGeom.py
Description: Series of functions to clean geometry from blocks, streets, and transport lines.
License: MIT, see full license in LICENSE.txt
Web: https://github.com/mateoneira/MultiplexSegregation
"""

import geopandas as gpd
import shapely.geometry as geometry
import numpy as np
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize, linemerge
import math


def get_vertex_of_polygons(geom):
    """
    Get list of vertices of all polygons in geoseries and return as list of points.
    If no polygons are supplied in geoseries empty list is returned.

    Parameters
    ----------
    :param geom: geopandas.GeoSeries
        geometries of city blocks.

    Returns
    -------
    :return: list
        list of vertices of polygons
    """

    if type(geom) != gpd.geoseries.GeoSeries:
        raise TypeError("geom should be a *geopandas.GeoSeries* type.")

    # get vertex of polygons.
    points = []
    for poly in geom:
        # check if geometry is polygon of multipolygon
        # if polygon add vertices to points list
        if poly.type == 'Polygon':
            for pnt in poly.exterior.coords:
                points.append(geometry.Point(pnt))
        elif poly.type == 'MultiPolygon':
            for parts in poly:
                for pnt in parts.exterior.coords:
                    points.append(geometry.Point(pnt))
    return points


def alpha_shape(points, alpha):
    """
    Calculate alpha shape from set of points and alpha value.

    Parameters
    ----------
    :param points: list
        list containing shapely.Geometry.Point objects
    :param alpha: float
        alpha value greater than 0

    Returns
    -------
    :return: shapely.geometry
    """
    if not all(isinstance(x, geometry.point.Point) for x in points):
        raise TypeError("points list must contain only *geometry.Point* type.")
    if alpha <= 0:
        raise ValueError("alpha must be greater than zero.")
    if len(points) < 3:
        raise TypeError("points list must have at least 3 items.")

    # create Delaunay triangulation
    coords = np.array([point.coords[0] for point in points])
    tri = Delaunay(coords)

    # create empty edge set and point list
    edges = set()
    edge_points = []

    # helper function to calculate which edges to keep
    def add_edge(i, j):
        if (i, j) in edges or (j, i) in edges:
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    for ia, ib, ic in tri.simplices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # calculate length of side of triangles
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        # calculate semiperimeter of triangle
        s = (a + b + c) / 2.0

        # calculate area of triangle
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))

        if area == 0:
            circum_r = 0
        elif area > 0:
            circum_r = a * b * c / (4.0 * area)
        else:
            pass

        # radius filter
        if circum_r < 1.0 / alpha:
            add_edge(ia, ib)
            add_edge(ib, ic)
            add_edge(ic, ia)

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    res = cascaded_union(triangles)

    return res


def boundary_from_areas(blocks, alpha=1, buffer_dist=0):
    """
    Create spatial boundary given unconnected block area geometries of
    city through an alpha shape.

    Parameters
    ----------
    :param blocks: geopandas.GeoDataFrame
        city block geometry
    :param alpha: float
        alpha value for alpha shape calculation
    :param buffer_dist: float
        distance to buffer alpha shape in meters.

    :return: geopandas.GeoSeries
    """

    if type(blocks) != gpd.geodataframe.GeoDataFrame:
        raise TypeError("blocks must be a *geopandas.GeoDataFrame*.")
    if alpha <= 0:
        raise ValueError("alpha must be an float greater than 0.")
    if buffer_dist < 0:
        raise ValueError("buffer_dist must be a float greater than 0.")

    # subset geometry from geodataframe.
    geom = blocks.geometry
    points = get_vertex_of_polygons(geom)

    # calculate alpha shape
    boundary = alpha_shape(points, alpha)

    # buffer alpha shape
    if buffer_dist > 0:
        boundary = boundary.buffer(buffer_dist)

    return gpd.GeoSeries(boundary)


def join_lines(lines, line_list):
    pass


def clean_stops(stops, tolerance=50):
    pass


def clean_lines(line):
    pass


def snap_stops_to_lines(line, stops, area, tolerance=50):
    pass


def snap_lines_to_points(G):
    pass


def cut(line, distance):
    pass


def find_nearest_node(data, nodes, spatial_index, buff=50):
    pass


def intersection_voronoi():
    pass


def area_overlay():
    pass
