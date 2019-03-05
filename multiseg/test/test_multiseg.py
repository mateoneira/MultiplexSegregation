from ..multiplexSeg import *
from ..processGeom import *

from pytest import raises, approx
import geopandas as gpd
import random


def test_import():
    """
    Test all module imports
    :return:
    """
    # test all ms module imports
    import pandas as pd
    # import geopandas as gpd
    import numpy as np
    import networkx as nx
    import osmnx as ox
    import shapely.geometry as geometry
    from shapely.ops import cascaded_union, polygonize, linemerge
    from shapely.wkt import loads
    import scipy as sp
    from scipy.spatial import Delaunay, Voronoi, voronoi_plot_2d
    import math
    import matplotlib.pyplot as plt
    from descartes import PolygonPatch
    import matplotlib.cm as cm
    from geopandas.tools import overlay
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    import matplotlib as mpl
    import matplotlib.colors as colors
    import seaborn as sns
    from scipy.sparse import identity, spdiags, linalg
    import time


def test_get_vertex_of_polygons():
    """
    Test get_vertex_of_polygons function takes in right values and creates correct output

    :return:
    """

    # check error is raised with incorrect input
    with raises(TypeError) as exception:
        get_vertex_of_polygons([1, 2, 1])

    # create geometries for testing
    block_1 = geometry.Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
    block_2 = geometry.Polygon([(2, 0), (2, 1), (3, 1), (3, 0)])
    block_3 = geometry.Polygon([(0, 2), (0, 3), (1, 3), (1, 2)])
    block_4 = geometry.Polygon([(2, 2), (2, 3), (3, 3), (3, 2)])
    multi_block = geometry.MultiPolygon([block_3, block_4])
    line_geom = geometry.LineString([(1.5, 0), (1.5, 3)])

    simple_box = gpd.GeoSeries(block_1)
    multi_box = gpd.GeoSeries([block_1, block_2, block_3, block_4])
    multi_box_2 = gpd.GeoSeries(multi_block)
    mixed_box = gpd.GeoSeries([block_1, block_2, multi_block])
    mixed_geom = gpd.GeoSeries([block_1, block_2, multi_block, line_geom])
    line = gpd.GeoSeries(line_geom)

    # test returns right type
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(simple_box))
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(multi_box))
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(multi_box_2))
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(mixed_box))
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(line))
    assert all(isinstance(x, geometry.point.Point) for x in get_vertex_of_polygons(mixed_geom))

    # test returns right number of elements
    assert len(get_vertex_of_polygons(simple_box)) == 5
    assert len(get_vertex_of_polygons(multi_box)) == 20
    assert len(get_vertex_of_polygons(multi_box_2)) == 10
    assert len(get_vertex_of_polygons(mixed_box)) == 20
    assert len(get_vertex_of_polygons(line)) == 0


def test_alpha_shape():
    """
    Test alpha_shape function takes in right values and creates correct output.

    :return:
    """

    # generate random points to test
    points = []
    for number in range(1,10):
        pnt = geometry.Point(random.uniform(0,10), random.uniform(0,10))
        points.append(pnt)

    # check error is raised with incorrect_input
    with raises(TypeError) as exception:
        alpha_shape([1, 2, 3], 1)

    with raises(ValueError) as exception:
        alpha_shape(points, 0)

    with raises(TypeError) as exception:
        alpha_shape(points[:1], 1)

    # check output is valid


def test_boundary_from_areas():
    """
    Test boundary_from_area function takes in right values and creates correct output.

    :return:
    """

    # create block gpd to test
    blocks = []
    ids = []
    x_min = 0
    n = 5
    for i in range(n):
        x_max = x_min + 100
        y_min = 0
        for j in range(0, n):
            y_max = y_min + 100
            blck = geometry.Polygon([(x_min, y_min),
                                     (x_max, y_min),
                                     (x_max, y_max),
                                     (x_min, y_max)
                                     ])
            blocks.append(blck)
            ids.append(i * n + j)
            y_min = y_max + 20
        x_min = x_max + 20

    blocks_gpd = gpd.GeoDataFrame({"id": ids, "geometry": blocks})

    # test inputs
    with raises(TypeError) as exception:
        boundary_from_areas([1,2,3])
    with raises(ValueError) as exception:
        boundary_from_areas(blocks_gpd, 0)
    with raises(ValueError) as exception:
        boundary_from_areas(blocks_gpd, 1, -100)

    assert len(boundary_from_areas(blocks_gpd, 0.01)) == 1






