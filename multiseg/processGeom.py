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

# Geographical projection of OpenStreetMap data.
crs_osm = {'init': 'epsg:4326'}


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


def join_lines(line, line_list, tolerance=20):
    """
    Join MultiLineStrings and returns SingleLineString through recursion.

    Parameters
    ----------
    :param line: list
        list of coordinates of LineString
    :param line_list: list
        list of list of coordinates of LineStrings
    :param tolerance: float
        tolerance of check if two points are the same point (in meters).

    Return
    ------
    :return: list
    """
    line_list = line_list.copy()
    # get last coordinate of line and make a point
    point_1 = geometry.Point(line[-1])

    # list to store coords list and their reverse
    coord_list = []

    if line_list is not None:
        for coords in line_list:
            # store all lines and reversed lines in one list
            coord_list.append(coords)
            coord_list.append(list(reversed(coords)))

        for coords in coord_list:
            point_2 = geometry.Point(coords[0])
            if point_1.distance(point_2) < tolerance+1:
                line_list.remove(coords)
                for coord in coords:
                    line.append(coord)
                join_lines(line, line_list)
            else:
                return line


def clean_stops(stops, tolerance=50):
    pass


def clean_lines(lines, group_by="None", tolerance=20):
    """
    Creates geodataframe containing geometries of LineString objects.
    MultiLineStrings and LineStringZ is converted to LineStrings.

    Parameters
    ----------
    :param lines: geopandas.GeoDataFrame containing transport line geometries
    :param group_by: str
        column name of group, if "None", the whole dataset is processed as one. Default "None".
    :param tolerance: float
        tolerance of check if two points are the same point (in meters).
    Returns
    -------
    :return: geopandas.GeoDataFrame
    """
    lines = lines.copy()

    # Empty list to store cleaned lines
    lines_list = []

    # Define how data will be subset to process
    if group_by == "None":
        lines['grouped'] = 0
    else:
        lines['grouped'] = lines[group_by]

    # loop through subset of data and join MultiLineString to SingleLineString
    for group in lines.grouped.unique():
        lines_subset = lines[lines.grouped == group]

        # loop through individual geometries
        for i, row in lines_subset.iterrows():
            geom = row.geometry

            # check if line is MultiLineString
            if isinstance(geom, geometry.MultiLineString):
                geom_list = geom.geoms

                # create empty list to store coordinates of line
                lines_coords = []
                for line in geom_list:
                    # if line is not smaller than tolerance meters and not a self-loop
                    if line.length > tolerance and line.coords[0] != line.coords[-1]:
                        if line.has_z:
                            coord_list = []
                            for coord in line.coords:
                                coord_list.append(coord[0:2])
                            lines_coords.append(coord_list)
                        else:
                            coord_list = list(line.coords)
                            lines_coords.append(coord_list)

                # choose first line and look for continuation
                line_coord = lines_coords[0]
                line_list = lines_coords[1:]

                line_joined = join_lines(line_coord, line_list)
                line_joined = join_lines(list(reversed(line_joined)), line_list)

                line_geom = geometry.LineString(coor for coor in line_joined)

            else:
                if geom.has_z:
                    coord_list = []
                    for coord in geom.coords:
                        coord_list.append(coord[0:2])
                    line_geom = geometry.LineString(coor for coor in coord_list)
                else:
                    line_geom = geom

            lines_list.append(line_geom)

    linesGPD = gpd.GeoDataFrame({'geometry': lines_list})

    return(linesGPD)


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
