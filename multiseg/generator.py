import geopandas as gpd
from shapely import geometry
import networkx as nx


# create block gpd to test

def generate_manhattan_blocks(blocksize=(100, 100), num_blocks=(5, 5), seperation=20):
    G = nx.Graph()
    blocks = []
    ids = []
    x_min = 0
    offset = seperation / 2
    for i in range(num_blocks[0]):
        x_max = x_min + blocksize[0]
        y_min = 0
        for j in range(0, num_blocks[1]):
            y_max = y_min + blocksize[1]
            blck = geometry.Polygon([(x_min + offset, y_min + offset),
                                     (x_max - offset, y_min + offset),
                                     (x_max - offset, y_max - offset),
                                     (x_min + offset, y_max - offset)
                                     ])
            blocks.append(blck)
            ids.append(i * num_blocks[0] + j)
            if y_min == 0 and x_min == 0:
                G.add_node("{}-{}".format(x_min, y_min), x=x_min, y=y_min, geometry=geometry.Point((x_min, y_min)))
                G.add_node("{}-{}".format(x_max, y_max), x=x_max, y=y_max, geometry=geometry.Point((x_max, y_max)))
                G.add_node("{}-{}".format(x_min, y_max), x=x_min, y=y_max, geometry=geometry.Point((x_min, y_max)))
                G.add_node("{}-{}".format(x_max, y_min), x=x_max, y=y_min, geometry=geometry.Point((x_max, y_min)))
                G.add_edge("{}-{}".format(x_min, y_min), "{}-{}".format(x_min, y_max))
                G.add_edge("{}-{}".format(x_min, y_max), "{}-{}".format(x_max, y_max))
                G.add_edge("{}-{}".format(x_max, y_max), "{}-{}".format(x_max, y_min))
                G.add_edge("{}-{}".format(x_max, y_min), "{}-{}".format(x_min, y_min))
            elif y_min != 0 and x_min == 0:
                G.add_node("{}-{}".format(x_min, y_max), x=x_min, y=y_max, geometry=geometry.Point((x_min, y_max)))
                G.add_node("{}-{}".format(x_max, y_max), x=x_max, y=y_max, geometry=geometry.Point((x_max, y_max)))
                G.add_edge("{}-{}".format(x_min, y_max), "{}-{}".format(x_max, y_max))
                G.add_edge("{}-{}".format(x_min, y_min), "{}-{}".format(x_min, y_max))
                G.add_edge("{}-{}".format(x_max, y_min), "{}-{}".format(x_max, y_max))
            elif y_min == 0 and x_min != 0:
                G.add_node("{}-{}".format(x_max, y_min), x=x_max, y=y_min, geometry=geometry.Point((x_max, y_min)))
                G.add_node("{}-{}".format(x_max, y_max), x=x_max, y=y_max, geometry=geometry.Point((x_max, y_max)))
                G.add_edge("{}-{}".format(x_max, y_min), "{}-{}".format(x_max, y_max))
                G.add_edge("{}-{}".format(x_min, y_min), "{}-{}".format(x_max, y_min))
                G.add_edge("{}-{}".format(x_min, y_max), "{}-{}".format(x_max, y_max))
            else:
                G.add_node("{}-{}".format(x_max, y_max), x=x_max, y=y_max, geometry=geometry.Point((x_max, y_max)))
                G.add_node("{}-{}".format(x_max, y_max), x=x_max, y=y_max, geometry=geometry.Point((x_max, y_max)))
                G.add_edge("{}-{}".format(x_max, y_max), "{}-{}".format(x_max, y_min))
                G.add_edge("{}-{}".format(x_min, y_max), "{}-{}".format(x_max, y_max))

            y_min = y_max
        x_min = x_max
    G.graph['crs'] = None
    blocks_gpd = gpd.GeoDataFrame({"id": ids, "geometry": blocks})
    boundary = gpd.GeoDataFrame(geometry=gpd.GeoSeries(geometry.Polygon([(0, 0),
                                                                         (0, num_blocks[1] * blocksize[1]),
                                                                         (num_blocks[0] * blocksize[0],
                                                                          num_blocks[1] * blocksize[1]),
                                                                         (num_blocks[0] * blocksize[0], 0)
                                                                         ]).buffer(50)))
    return blocks_gpd, G, boundary
