import multiseg as ms
from pytest import raises, approx
import geopandas as gpd

def test_import():
    #test all ms module imports
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
    Check get vertex function takes in right values and
    creates correct output
    :return:
    """

    ##check error is raised with incorrect input
    with raises(TypeError) as exception:
        ms.get_vertex_of_polygons([1,2,1])

    ##test only one element in geoseries

    ##test multiple geom types

    ##test only multipolygon

    ##test multipolygon and polygon





