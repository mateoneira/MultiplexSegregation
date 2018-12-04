"""
Module: processGeom.py
Description: Series of functions to clean geometry from blocks, streets, and transport lines.
License: MIT, see full license in LICENSE.txt
Web: https://github.com/mateoneira/MultiplexSegregation
"""

import geopandas as gpd
from geopandas.tools import overlay
import shapely.geometry as geometry
import numpy as np
from scipy.spatial import Delaunay, Voronoi
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


def clean_stops(stops, boundary, group_by=None, tolerance=50, data={}):
    """
    Create geodataframe containing point geometries representing stops in transport network.
    Points are clustered based on tolerance distance, and and centroid is returned as new point.

    Parameters
    ----------
    :param stops: geopandas.GeoDataFrame
        transport stops geometry.
    :param boundary: geopandas.GeoDataFrame
        geodataframe of boundary polygon.
    :param group_by:  str
        column name of group, if None, the whole dataset is processed as one. Default None.
    :param tolerance: float
        tolerance to check if two points are the sme point (in meters).
    :param data: dict
        data that has to be retained and mapping to ne column name.

    Returns
    -------
    :return: geopandas.GeoDataFrame
    """
    temp = []
    stops = stops.copy()
    boundary_geom = boundary.unary_union

    # check if data values need to be conserved
    mapped_data = {new_column: [] for (old_column, new_column) in data.items()}
    if 'geometry' not in mapped_data.keys():
        mapped_data['geometry'] = []

    # Define how data will be subset to process
    if group_by is None:
        stops['grouped'] = 0
    else:
        stops['grouped'] = stops[group_by]

    mapped_data['grouped'] = []

    # loop through groups, buffer, join, and append new point
    for group in stops.grouped.unique():
        stops_subset = stops[stops.grouped == group]

        buffered_stops = stops_subset.buffer(tolerance).unary_union

        # check if new geom is polygon, and convert to list
        if isinstance(buffered_stops, geometry.Polygon):
            buffered_stops = [buffered_stops]

        for geom in buffered_stops:
            mapped_data['grouped'].append(group)
            mapped_data['geometry'].append(geom.centroid)

            # map data from points to centroids
            if data:
                temp = stops_subset[stops_subset.intersects(geom)]

                for column_name, new_column in data.items():
                    val = ', '.join(str(v) for v in temp[column_name].unique())
                    mapped_data[new_column].append(val)

    stopsGPD = gpd.GeoDataFrame(mapped_data)
    stopsGPD = stopsGPD[stopsGPD.intersects(boundary_geom)]
    return stopsGPD


def clean_lines(lines, group_by=None, tolerance=20, data={}):
    """
    Creates geodataframe containing geometries of LineString objects.
    MultiLineStrings and LineStringZ is converted to LineStrings.

    Parameters
    ----------
    :param lines: geopandas.GeoDataFrame
        transport line geometries
    :param group_by: str
        column name of group, if None, the whole dataset is processed as one. Default None.
    :param tolerance: float
        tolerance of check if two points are the same point (in meters).
    :param data: dict
        data that has to be retained and mapping to new column name.
    Returns
    -------
    :return: geopandas.GeoDataFrame
    """
    lines = lines.copy()

    # check if data values need to be conserved
    mapped_data = {new_column: [] for (old_column, new_column) in data.items()}
    if 'geometry' not in mapped_data.keys():
        mapped_data['geometry'] = []

    # Define how data will be subset to process
    if group_by is None:
        lines['grouped'] = 0
    else:
        lines['grouped'] = lines[group_by]

    mapped_data['grouped'] = []

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

            mapped_data['geometry'].append(line_geom)
            mapped_data['grouped'].append(row['grouped'])

            # map values
            for column_name, new_column in data.items():
                mapped_data[new_column].append(row[column_name])

    linesGPD = gpd.GeoDataFrame(mapped_data)

    return linesGPD


def snap_stops_to_lines(lines, stops, tolerance=50):
    """
    Snaps points to lines based on tolerance distance and route.

    Parameters
    ----------
    :param lines: geopandas.GeoDataFrame
        geodataframe containing line geometries.
    :param stops: geopandas.GeoDataFrame
        geodataframe containing stop geometries.
    :param tolerance: float
        distance tolerance for snapping points (in meters).

    Returns
    -------
    :return: geopandas.GeoDataFrame
        geodataframe with point geometries snapped to closest transport route.
    """

    snapped_stops = gpd.GeoDataFrame()

    for group in lines.grouped.unique():
        lines_subset = lines[lines.grouped == group]
        stops_subset = stops[stops.grouped == group]

        # snap points to lines
        for i, line in lines_subset.iterrows():
            geom = line.geometry

            # get only points within buffer and inside area
            buffer = geom.buffer(tolerance)

            stops_inside = stops_subset[stops_subset.intersects(buffer)].copy()

            points_proj = [geom.project(stop) for stop in stops_inside.geometry]
            stops_inside.geometry = [geom.interpolate(point) for point in points_proj]
            stops_inside['at_length'] = points_proj
            stops_inside['line_id'] = [i for point in points_proj]
            snapped_stops = snapped_stops.append(stops_inside, ignore_index=True)

    snapped_stops = snapped_stops.drop_duplicates(subset=[col for col in snapped_stops.columns if col != 'geometry'])
    snapped_stops = snapped_stops.dropna(how="all")
    snapped_stops['stop_id'] = [i for i in range(len(snapped_stops))]
    snapped_stops['x'] = [point.xy[0][0] for point in snapped_stops.geometry]
    snapped_stops['y'] = [point.xy[1][0] for point in snapped_stops.geometry]
    return snapped_stops


def snap_lines_to_points(G):
    pass


def cut_line(line, distance):
    """
    Cuts line at a set distance.

    Parameters
    ----------
    :param line: shapely.LineString
        line geometry to cut.
    :param distance: float
        distance at which to cut line.

    Returns
    -------
    :return: list
        list containing line segments resultant from the cut.
    """

    if distance <= 0.0 or distance >= line.length:
        return [line]

    coords = list(line.coords)

    for i, p in enumerate(coords):
        current_distance =line.project(geometry.Point(p))

        if current_distance == distance:
            return [geometry.LineString(coords[:i+1]), geometry.LineString(coords[i:])]
        elif current_distance>distance:
            cut_point = line.interpolate(distance)

            return [geometry.LineString(coords[:i+1] + [(cut_point.x, cut_point.y)]),
                    geometry.LineString([(cut_point.x, cut_point.y)] + coords[i:])]


def find_nearest_node(data, nodes, spatial_index, buff=50):
    pass


def create_node_voronoi(G, boundary):
    """


    Parameters
    ----------
    :param G: networkx.MultiDiGraph
        Network for which to create a node tessellation.
    :param boundary: geopandas.GeoDataFrame
        Boundary polygon.

    Returns
    -------
    :return:
    """
    node_x = [float(node['x']) for node in G.node.values()]
    node_y = [float(node['y']) for node in G.node.values()]

    points = np.column_stack((node_x,node_y))

    tessellation = Voronoi(points)

    # create polygon from voronoi tessellation
    lines = [geometry.LineString(tessellation.vertices[line])
             for line in tessellation.ridge_vertices if -1 not in line]

    polygons = list(polygonize(lines))
    polygonsGPD = gpd.GeoDataFrame(geometry=polygons)
    polygonsGPD.crs = G.graph['crs']

    # create intersection with urban limit
    polygonsGPD = gpd.overlay(polygonsGPD, boundary)
    polygonsGPD = polygonsGPD[polygonsGPD.is_valid]

    return polygonsGPD


def area_overlay(sources, targets, population, indices = [], groups = []):
    """
    Calculate area overlay given to geometries and initial values.

    Parameters
    ----------
    :param sources: geopandas.GeoDataFrame
    :param targets: geopandas.GeoDataFrame

    Returns
    -------
    :return:
    """
    new_targets = targets.copy()
    population_data = []
    indices_data = {index: [] for index in indices}
    groups_data = {group: [] for group in groups}

    for i, target in new_targets.iterrows():
        temp_population = 0
        temp_indices = {index: 0 for index in indices}
        temp_groups = {group: 0 for group in groups}
        count = 0
        weight = 0
        target_geom = target.geometry

        # create spatial index and find geometry within polygon
        sindex = sources.sindex
        matches_index = list(sindex.intersection(target_geom.bounds))

        for matched_index in matches_index:
            intersection = overlay(sources.iloc[[matched_index]], new_targets.iloc[[i]], how='intersection')

            if len(intersection) is not 0:
                count += 1
                source_area = sum(sources.iloc[[matched_index]].area)
                inters_area = sum(intersection.geometry.area)
                inters_ratio = inters_area/source_area

                # for population use weighted sum
                population_value = inters_ratio * sources.iloc[[matched_index]][population].values[0]

                temp_population += population_value
                # for indices use weighted mean
                weight += inters_ratio
                for index in temp_indices.keys():
                    temp_indices[index] += inters_ratio * sources.iloc[[matched_index]][index].values[0]

                # groups must be recalculated based weighted population and percentage
                for group in temp_groups.keys():
                    temp_groups[group] += population_value * sources.iloc[[matched_index]][group].values[0]

        if temp_population != 0:
            for group in groups_data.keys():
                groups_data[group].append(temp_groups[group]/temp_population)
        else:
            for group in groups_data.keys():
                groups_data[group].append(0)

        if count != 0:
            for index in indices_data.keys():
                indices_data[index].append(temp_indices[index]/weight)
        else:
            for index in indices_data.keys():
                indices_data[index].append(np.nan)
        population_data.append(temp_population)

    # append values to target geometry
    new_targets[population] = population_data

    for index in indices_data.keys():
        new_targets[index] = indices_data[index]

    for group in groups_data.keys():
        new_targets[group] = groups_data[group]

    return new_targets




