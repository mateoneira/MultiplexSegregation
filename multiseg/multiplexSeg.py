from .processGeom import *
import pandas as pd
import networkx as nx
import osmnx as ox
from descartes import PolygonPatch
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy.sparse import identity, spdiags, linalg
import scipy as sp


from mpl_toolkits.mplot3d.art3d import Line3DCollection


def transport_graph_from_lines_stops(linesGPD, stopsGPD, boundary, speed=30, name='', group_lines_by=None,
                                     group_stops_by=None,
                                     lines_data=None, stops_data=None):
    """
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
    :param name:

    Returns
    -------
    :return: networkx.MultiDiGraph
    """
    # clean input data
    lines = clean_lines(linesGPD, group_by=group_lines_by, data=lines_data)
    stops = clean_stops(stopsGPD, boundary, group_by=group_stops_by, data=stops_data)
    stops = snap_stops_to_lines(lines, stops)

    ujt = 1/(speed * 16.6666666667)  # turn Km/h to min/meter

    lines_dict = {column: [] for column in lines.columns}
    lines_dict['length'] = []
    lines_dict['u'] = []
    lines_dict['v'] = []
    lines_dict['from'] = []
    lines_dict['to'] = []
    lines_dict['weight'] = []

    for i, line in lines.iterrows():
        geom = line.geometry
        stops_on_line = stops[stops.line_id == i].copy()

        stops_on_line = stops_on_line.sort_values(by='at_length').reset_index()

        for j, stop in stops_on_line[:-1].iterrows():
            start_pnt = stop.at_length
            end_pnt = stops_on_line.at_length[j+1]

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
            if (stop.geometry.distance(geometry.Point(segment.coords[0])) < 1)\
                    and (stops_on_line.geometry[j+1].distance(geometry.Point(segment.coords[-1])) < 1):
                for column in lines.columns:
                    if column != 'geometry':
                        lines_dict[column].extend([line[column], line[column]])
                lines_dict['geometry'].extend([segment, segment])
                lines_dict['length'].extend([distance_between, distance_between])
                lines_dict['u'].extend([stop.stop_id, stops_on_line.stop_id[j+1]])
                lines_dict['v'].extend([stops_on_line.stop_id[j+1], stop.stop_id])
                lines_dict['from'].extend([segment.coords[0], segment.coords[-1]])
                lines_dict['to'].extend([segment.coords[-1], segment.coords[0]])
                lines_dict['weight'].extend([distance_between * ujt, distance_between * ujt])

    edge_list = gpd.GeoDataFrame(lines_dict)
    edge_list['key'] = [key for key in range(len(edge_list))]
    edge_list.crs = linesGPD.crs

    stops.gdf_name = 'node_list'
    edge_list.gdf_name = 'edge_list'

    G = ox.gdfs_to_graph(stops, edge_list)
    G = add_transfer_edges(G, 15)
    G.name = name

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


def find_nearest_node(node, node_list, spatial_index, search_radius=50):
    """

    :param node: dict
    :param node_list:
    :param spatial_index:
    :return:
    """

    polygon = node['geometry'].buffer(search_radius)
    possible_matches_index = list(spatial_index.intersection(polygon.bounds))
    possible_matches = node_list.iloc[possible_matches_index]

    if len(possible_matches) == 0:
        search_radius += search_radius
        v = find_nearest_node(node, node_list, spatial_index, search_radius)
    elif len(possible_matches) == 1:
        v = possible_matches.index[0]
        return v
    elif len(possible_matches) > 1:
        point = node['geometry']
        distance = possible_matches.distance(point)
        v = distance[distance == min(distance)].index[0]
        return v
    else:
        print('no match')

    return v


def map_values_to_nodes(G, spatial_units, boundary, population, indices, groups):
    """

    :param G:
    :param spatial_units:
    :param boundary:
    :return:
    """
    G = G.copy()
    tessellation = create_node_voronoi(G, boundary)
    overlay_poly = area_overlay(spatial_units, tessellation, population, indices, groups)

    # assign values to nodes in graph
    nodes = {node: data for node, data in G.nodes(data=True)}
    nodes_gpd = gpd.GeoDataFrame(nodes).T
    nodes_gpd.crs = G.graph['crs']

    nodes_attrib = gpd.sjoin(nodes_gpd, overlay_poly, how='inner', op='within')
    nodes_dict = nodes_attrib.to_dict(orient='index')
    nx.set_node_attributes(G, nodes_dict)

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
        ujt = 1 / (speed * 16.666666667)  # turn km/h to min/meter
        G[u][v][key]['weight'] = data * ujt

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


def create_bipartite_graph(G_1, G_2, weight=(15, 0.1)):
    """

    :param G_1:
    :param G_2:
    :param weight:
    :return:
    """

    # create an empty graph
    B_1 = nx.DiGraph()
    B_2 = nx.DiGraph()

    # get node dicts
    nodes_1 = gpd.GeoDataFrame(dict(G_1.nodes(data=True))).T
    nodes_2 = gpd.GeoDataFrame(dict(G_2.nodes(data=True))).T

    spatial_index = nodes_1.sindex

    # iterate and check
    for u, node in nodes_2.iterrows():
        v = find_nearest_node(node, nodes_1, spatial_index)

        # add nodes
        B_2.add_node(u, bipartite=1, layer=G_2.name)
        B_2.add_node(v, bipartite=0, layer=G_1.name)
        B_1.add_node(u, bipartite=1, layer=G_2.name)
        B_1.add_node(v, bipartite=0, layer=G_1.name)

        # add edge
        B_2.add_edge(u, v, transfer=True, weight=weight[1])
        B_1.add_edge(v, u, transfer=True, weight=weight[0])

    return B_1, B_2


def add_transfer_edges(G, weight=15):
    """
    Create transfer edge based on proximity of nodes.

    Parameters
    ----------
    :param G: networkx.MultiDiGraph
    :param weight: float

    Returns
    -------
    :return:
    """
    G = G.copy()
    nodes = {node: data for node, data in G.nodes(data=True)}
    nodes_gpd = gpd.GeoDataFrame(nodes).T

    nodes_buffer = nodes_gpd.buffer(15).unary_union
    nodes_buffer = gpd.GeoSeries(list(nodes_buffer))

    # create edge for all points within the buffer
    for area in nodes_buffer:
        nodes_in_buffer = nodes_gpd[nodes_gpd.intersects(area)]
        new_node = area.centroid

        for i in nodes_in_buffer.index:
            G.node[i]['x'] = new_node.xy[0][0]
            G.node[i]['y'] = new_node.xy[1][0]
            G.node[i]['geometry'] = new_node

            for j in nodes_in_buffer.index:
                if i != j:
                    G.add_edge(i,
                               j,
                               way='transfer',
                               key=0,
                               weight=weight)
    return G


def plot_multiplex(M):
    """

    :param M:
    :return:
    """

    node_Xs = [float(node['x']) for node in M.node.values()]
    node_Ys = [float(node['y']) for node in M.node.values()]
    node_Zs = np.array([float(node['z'])*1000 for node in M.node.values()])

    node_size = []
    size = 1
    node_color = []

    for i, d in M.nodes(data=True):
        if d['z'] == 0:
            node_size.append(size)
            node_color.append('#66ccff')
        if d['z'] == 1:
            node_size.append(size*4)
            node_color.append('#fb9a98')
        if d['z'] == 2:
            node_size.append(size*8)
            node_color.append('#9bc37e')

    lines = []
    line_width = []
    lwidth = 0.2

    for u, v, key, data in M.edges(keys=True, data=True):
        if 'geometry' in data:
            xs, ys = data['geometry'].xy
            zs = [M.node[u]['z']*1000 for i in range(len(xs))]

            lines.append([list(a) for a in zip(xs, ys, zs)])

            if 'transfer' in data.keys():
                line_width.append(lwidth/1.5)
            else:
                line_width.append(lwidth)
        else:
            x1 = M.node[u]['x']
            y1 = M.node[u]['y']
            z1 = M.node[u]['z']*1000
            x2 = M.node[v]['x']
            y2 = M.node[v]['y']
            z2 = M.node[v]['z']*1000

            line = [[x1, y1, z1], [x2, y2, z2]]

            lines.append(line)

            if 'transfer' in data.keys():
                line_width.append(lwidth/1.5)
            else:
                line_width.append(lwidth)

    fig_height = 15

    lc = Line3DCollection(lines, linewidths=line_width, alpha=1, color="#999999", zorder=1)

    edges = ox.graph_to_gdfs(M, nodes=False, fill_edge_geometry=True)
    west, south, east, north = edges.total_bounds
    bbox_aspect_ratio = (north - south) / (east - west)
    fig_width = fig_height / bbox_aspect_ratio
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.gca(projection='3d')
    ax.add_collection3d(lc)
    ax.scatter(node_Xs, node_Ys, node_Zs, s=node_size, c=node_color, zorder=2)
    ax.set_ylim(south, north)
    ax.set_xlim(west, east)
    ax.set_zlim(0, 2500)
    ax.axis('off')
    ax.margins(0)
    ax.tick_params(which='both', direction='in')
    fig.canvas.draw()

    ax.set_facecolor('white')
    ax.set_aspect('equal')

    plt.show()


def plot_supra_adjacency(M):
    """

    :param M:
    :return:
    """
    fig = plt.figure(figsize=(15, 10))

    partitions = [[node for node, data in M.nodes(data=True) if data['layer'] == layer] for layer in M.layers()]

    node_order = [node for partition in partitions for node in partition]

    supra_adjacency = M.supra_adjacency(node_order=node_order)

    plt.imshow(supra_adjacency, cmap="Greys", interpolation="none")

    ax = plt.gca()
    current_idx = 0
    for partition in partitions:

        ax.add_patch(patches.Rectangle((current_idx, current_idx),
                                       len(partition), len(partition),
                                       facecolor="none",
                                       edgecolor="Blue",
                                       linewidth=1))
        current_idx += len(partition)


def plot_adjacency(G):
    fig = plt.figure(figsize=(15, 10))

    adjacency = nx.to_numpy_matrix(G, dtype=np.bool)

    plt.imshow(adjacency, cmap="Greys", interpolation="none")


def random_walk_segregation(M, group="Q1", alpha=0.85, walk_type="RW", weight=None):
    """

    :param M:
    :param group:
    :param alpha:
    :return:
    """
    G = M
    for i, data in G.nodes(data=True):
        if 'nPeople' not in data.keys():
            G.node[i]['nPeople'] = 0
        if group not in data.keys():
            G.node[i][group] = 0

    # get initial values as column vectors
    n_i = np.array(list(nx.get_node_attributes(G, 'nPeople').values()))
    n = n_i.sum()

    c_gi = np.array(list(nx.get_node_attributes(G, group).values()))
    n_gi = c_gi * n_i
    n_g = n_gi.sum()
    d_gi = n_gi / n_g

    c_gi.shape = (len(c_gi), 1)
    d_gi.shape = (len(d_gi), 1)

    if walk_type == "RW":
        adj = nx.adjacency_matrix(G, weight=None)
        degree = sp.sparse.spdiags(adj.sum(1).T, 0, *adj.shape)
        identity_matrix = sp.identity(len(G))

        P = linalg.inv(degree) * adj
        Q = (1 - alpha) * (identity_matrix - alpha * P).I * P

    if walk_type == "PN":
        W = nx.adjacency_matrix(G, weight=weight)
        degree = sp.sparse.spdiags(W.sum(1).T, 0, *W.shape)
        identity_matrix = sp.identity(len(G))

        P = linalg.inv(degree) * W
        Q = (1 - alpha) * (identity_matrix - alpha * P).I * P

    if walk_type == "LF":
        node_order = list(G.nodes())
        # create distance matrix
        dist = dict(nx.all_pairs_dijkstra_path_length(G))

        # turn dict of dict into matrix
        list_of_dist = []
        for node_i in node_order:
            distances = []
            for node_j in dist[node_i].keys():
                distances.append(dist[node_i][node_j])
            else:
                distances.append(0)
            list_of_dist.append(distances)

        dst_matrix = np.matrix(list_of_dist, dtype='float')

        # probability transition matrix based on distance
        levy = np.power(dst_matrix, -alpha, out=np.zeros_like(dst_matrix), where=dst_matrix != 0)
        degree = sp.sparse.spdiags(levy.sum(1).T, 0, *levy.shape)
        P = linalg.inv(degree) * sp.sparse.csr_matrix(levy)

    isolation_gi = np.multiply(d_gi, (Q * c_gi))

    norm_isolation_gi = (n_g / n) ** -1 * isolation_gi * len(isolation_gi)
    print('norm seg for group {} = {}'.format(group, norm_isolation_gi.mean()))
    sigma_bar = list(np.array(norm_isolation_gi.flatten())[0])
    res = dict(zip(list(G.nodes()), sigma_bar))





