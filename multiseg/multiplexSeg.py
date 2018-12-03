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
    stops = clean_stops(stopsGPD, boundary, group_by=group_stops_by, data=stops_data)
    stops = snap_stops_to_lines(lines, stops)

    UJT = 1/(speed * 16.6666666667) #turn Km/h to min/meter

    lines_dict = {column: [] for column in lines.columns}
    lines_dict['length'] = []
    lines_dict['u'] = []
    lines_dict['v'] = []
    lines_dict['from'] =[]
    lines_dict['to'] = []
    lines_dict['weight'] = []

    for i, line in lines.iterrows():
        geom = line.geometry
        stops_on_line = stops[stops.line_id == i].copy()

        stops_on_line = stops_on_line.sort_values(by = 'at_length').reset_index()

        for i, stop in stops_on_line[:-1].iterrows():
            start_pnt = stop.at_length
            end_pnt = stops_on_line.at_length[i+1]

            distance_between = end_pnt - start_pnt

            # cut line
            first_segments = cut_line(geom, start_pnt)

            if len(first_segments) == 1:
                temp_segment = first_segments[0]
            else:
                temp_segment = first_segments[1]

            # make second cut
            segment = cut_line(temp_segment, distance_between)[0]

            # if stops are endpoints of new lines, save
            if (stop.geometry.distance(geometry.Point(segment.coords[0]))<1) and (stops_on_line.geometry[i+1].distance(geometry.Point(segment.coords[-1]))<1):
                for column in lines.columns:
                    if column != 'geometry':
                        lines_dict[column].extend([line[column], line[column]])
                lines_dict['geometry'].extend([segment, segment])
                lines_dict['length'].extend([distance_between, distance_between])
                lines_dict['u'].extend([stop.stop_id, stops_on_line.stop_id[i+1]])
                lines_dict['v'].extend([stops_on_line.stop_id[i+1], stop.stop_id])
                lines_dict['from'].extend([segment.coords[0], segment.coords[-1]])
                lines_dict['to'].extend([segment.coords[-1], segment.coords[0]])
                lines_dict['weight'].extend([distance_between * UJT, distance_between * UJT])

    edge_list = gpd.GeoDataFrame(lines_dict)
    edge_list['key'] = [key for key in range(len(edge_list))]
    edge_list.crs = linesGPD.crs

    stops.gdf_name = 'node_list'
    edge_list.gdf_name = 'edge_list'

    G = ox.gdfs_to_graph(stops, edge_list)
    
    return G




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
