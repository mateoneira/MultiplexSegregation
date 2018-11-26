import geopandas as gpd
import shapely.geometry as geometry
import numpy as np
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize, linemerge

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

def alpha_shape(points, alpha):
    """
    Calculate alpha shape from set of points and alpha value.
    :param points: list
        list containing shapely.Geometry.Point objects
    :param alpha: int
        alpha value
    :return: shapely.geometry.Polygon
    """
    #create Delaunay triangulation
    coords=np.array([point.coords[0] for point in points])
    tri=Delaunay(coords)

    #create empty edge set and point list
    edges=set()
    edge_points=[]

    ##helper function to calculate which edges to keep
    def add_edge(edge, edge_points, coords, i, j):
        if (i, j) in edges or (j, i) in edges:
            return
        edge.add((i, j))
        edge_points.append(coords[[i, j]])

    for ia, ib, ic in tri.vertices:
        pa=coords[ia]
        pb=coords[ib]
        pc=coords[ic]

        #calculate length of side of triangles
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        #calculate semiperimeter of triangle
        s=(a+b+c)/2.0

        #calculate area of triangle
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))

        if area ==0:
            circum_r=0
        elif area>0:
            circum_r = a * b * c / (4.0 * area)

        #radius filter
        if circum_r<1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)

    m=geometry.MultiLineString(edge_points)
    triangles=list(polygonize(m))
    res=cascade_union(triangles)

    return res





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


