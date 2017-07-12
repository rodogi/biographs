from .bpdb import *
import networkx as nx

def network(model, cutoff=5, weight=True):

    adjacency_dictionary = residue_adjacency(model, cutoff=cutoff,
                                            weight=weight)

    return nx.Graph(adjacency_dictionary)
