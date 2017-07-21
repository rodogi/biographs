import networkx as nx
from biographs.lib.bpdb import residue_adjacency


def network(model, cutoff=5, weight=True):

    adjacency_dictionary = residue_adjacency(model, cutoff=cutoff,
                                            weight=weight)

    return nx.Graph(adjacency_dictionary)
