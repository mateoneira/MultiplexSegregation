from networkx import MultiDiGraph, disjoint_union_all, set_node_attributes, to_numpy_matrix
import pandas as pd
import numpy as np

class MultilayerGraph(MultiDiGraph):
    """
    Multilayer definition
    """

    def __init__(self, layers=None, bipartite_interlayers=None):
        # set basic information
        if layers:
            self._layers = [layer.copy() for layer in layers]
            self._layers_name = [layer.name for layer in layers]

        else:
            self._layers = []

        for i, layer in enumerate(self._layers):
            set_node_attributes(layer, layer.name, 'layer')
            set_node_attributes(layer, {node_id: node_id for node_id in list(layer.nodes)}, 'layer_node_id')
            set_node_attributes(layer, i, 'z')

        if self._layers:
            super().__init__(disjoint_union_all(self._layers))
        else:
            super().__init__()

        self._nodes_pd = pd.DataFrame(dict(self.nodes(data=True))).T

        if bipartite_interlayers:
            self._bipartite_interlayers = [bipartite_interlayer.copy() for bipartite_interlayer in bipartite_interlayers]

        else:
            self._bipartite_interlayers = []

        for B in self._bipartite_interlayers:
            for u, v, data in B.edges(data=True):
                u_layer = B.nodes[u]['layer']
                v_layer = B.nodes[v]['layer']

                u_m = self._nodes_pd[(self._nodes_pd.layer == u_layer) & (self._nodes_pd.layer_node_id == u)].index[0]
                v_m = self._nodes_pd[(self._nodes_pd.layer == v_layer) & (self._nodes_pd.layer_node_id == v)].index[0]

                self.add_edge(u_m, v_m, 1, **data)

    def number_of_layers(self):
        """
        Return the number of layers in the Multilayer Graph
        :return: int
        """

        return len(self._layers)

    def number_of_interlayer_edges(self):
        return sum([B.number_of_edges() for B in self._bipartite_interlayers])

    def number_of_intralayer_edges(self):
        return sum([G.number_of_edges() for G in self._layers])

    def supra_adjacency(self, node_order = None):
        return to_numpy_matrix(self, nodelist=node_order, dtype=np.bool)

    def layers(self):
        return self._layers_name

    # def add_inter



