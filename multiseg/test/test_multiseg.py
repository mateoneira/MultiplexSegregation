import multiseg as ms

def test_import():
    #test all ms module imports
    import pandas as pd
    import geopandas as gpd
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

