import geopandas as gpd
import shapely.geometry as geometry
import numpy as np
from scipy.spatial import Delaunay

def get_vertex_of_polygons(geom):
    """
    Get list of vertices of all polygons in polygon list
    :param geom: geopandas.GeoSeries
    :return: list
        list of vertices of polygons
    """

    if type(geom)!= gpd.geoseries.GeoSeries:
        raise TypeError("geom should be a *geopandas.GeoSeries* type.")

    #get vertex of polygons.
    points = []
    for poly in geom:
        #check if geometry is polygon of multipolygon
        # if polygon add vertices to points list
        if poly.type=='Polygon':
            for pnt in poly.exterior.coords:
                points.append(geometry.Point(pnt))
        elif poly.type=='MultiPolygon':
            for parts in poly:
                for pnt in parts.exterior.coords:
                    points.append(geometry.Point(pnt))
    return points


# def add_edge(edge, edge_points, coords, i, j):
#     if (i, j) in edges or (j, i) in edges:
#         return
#     edge.add((i, j))
#     edge_points.append(coords[[i, j]])

def boundary_from_areas(blocks, alpha=1, buffer_dist=0):
    """
    Create spatial boundary given unconnected block area geometries of
    city through an alpha shape.

    Parameters
    ----------
    :param blocks: geopandas.GeoDataFrame
        city block geometry
    :param alpha: int
        alpha value for alpha shape calculation
    :param buffer_dist: int
        distance to buffer alpha shape in meters.

    :return: geopandas.GeoSeries
    """

    #subset geometry from geodataframe.
    geom =blocks.geometry
    points=get_vertex_of_polygons(geom)

    #create delaunay triangulation of points
    coords=np.array([point.coords[0] for point in points])
    tri=Delaunay(coords)
    edges(set)
    edge_points=[]


