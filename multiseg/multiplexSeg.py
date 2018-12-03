from .processGeom import *
import pandas as pd
import networkx as nx
import osmnx as ox
from descartes import PolygonPatch
import matplotlib.pyplot as plt


def transport_graph_from_lines_stops(linesGPD, stopsGPD, boundary, speed=30, group_lines_by=None, group_stops_by=None,
                                     lines_data={}, stops_data={}):
    """=
    Create a transport graph from line and stop shapefile data within
    a spatial boundary.

    The transport network is created as a time weighted multidigraph, all bus lines
    are modelled as bi-directional edges to represent the network.

    Transfers between stations are modelled as time weighted edges.

    Parameters
    ----------
    :param linesGPD:
    :param stopsGPD:
    :param boundary:
    :param speed:
    :param group_lines_by:
    :param group_stops_by:
    :param lines_data:
    :param stops_data:

    Returns
    -------
    :return: networkx.MultiDiGraph
    """
    # clean input data
    lines = clean_lines(linesGPD, group_by=group_lines_by, data=lines_data)
    stops = clean_stops(stopsGPD, group_by=group_stops_by, data=stops_data)

    stopsGPD = snap_stops_to_lines(linesGPD, stopsGPD, boundary)

    UJT = 1/(speed * 16.6666666667) #turn Km/h to min/meter

    line_list = []


def street_graph_from_boundary(boundary, speed=5, network_type="walk"):
    """
    Create a street network from OSM data within the spatial boundary.

    The street network is created as a time weighted multidigraph,
    all streets are modelled as undirected edges to represent a pedestrian network.

    Parameters
    ----------
        boundary: geopandas.GeoDataFrame
            boundary for which to download data.
        crs: dict
            projection system of blocks geometry.
        speed: float
            speed for calculating weights of edges (km/h).

    Returns
    -------
    :return: networkx.MultiDiGraph

    """

    # re-project boundary to projection system used on OSM and get only geometry
    geoms = boundary.to_crs(crs_osm).unary_union

    print('Generating street graph using OSMnx')
    G = ox.graph_from_polygon(geoms, network_type, name='street network')
    G = ox.project_graph(G, to_crs=boundary.crs)
    G = to_time_weighted(G, speed)
    G = add_geometry_to_network(G)

    return G


def add_geometry_to_network(G):
    """
    Convert x, y node attributes to shapely.Point object

    Parameters
    ----------
    :param G: networkx.MultiDiGraph
        Network with x, y in node attributes

    Returns
    -------
    :return: networkx.MultiDiGraph
    """
    # get x y coordinates of nodes
    node_id = [u for u in G.nodes(data=False)]
    geoms = [geometry.Point(data['x'], data['y']) for u, data in G.nodes(data=True)]
    node_dict = dict(zip(node_id, geoms))
    nx.set_node_attributes(G, node_dict, "geometry")
    print("created street network graph")
    return G


def to_time_weighted(G, speed):
    """
    Convert from distance weighted network to time weighted network.

    Parameters
    ----------
    :param G: networkx.MultiDiGraph
        transport network with length attribute.
    :param speed: float
        speed in Km/h.

    Parameters
    ----------
    :return: networkx.MultiDiGraph
        time-weighted network
    """
    for u, v, key, data in G.edges(data='length', keys=True):
        UJT = 1 / (speed * 16.666666667)  # turn km/h to min/meter
        G[u][v][key]['weight'] = data * UJT

    return G


def plot_network(G, boundary):
    """
    Create plot for network and urban boundary.

    Parameters
    ----------
    :param G: networkx.MultiDiGraph
    :param boundary: geopandas.GeoDataFrame

    Returns
    -------
    :return:
    """

    fig, ax = ox.plot_graph(G,
                            fig_height=15,
                            node_size=10,
                            node_zorder=1,
                            edge_color="#333333",
                            edge_linewidth=0.5,
                            save=False,
                            show=False,
                            close=False,
                            bgcolor="w",
                            dpi=1200,
                            equal_aspect=True
                            )

    urban_outline = PolygonPatch(boundary.geometry[0],
                                 fc="w",
                                 ec="r",
                                 linewidth=1,
                                 alpha=0.5,
                                 zorder=-1
                                 )

    ax.add_patch(urban_outline)
    margin = 0.02
    west, south, east, north = boundary.geometry[0].bounds
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))
    plt.show()
