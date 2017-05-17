import networkx as nx
from bpdb import *

def network(structure_file, cutoff=5, weight=True):

    adjacency_dictionary = residue_adjacency(structure_file, cutoff=cutoff,
                                            weight=weight)

    return nx.Graph(adjacency_dictionary)
