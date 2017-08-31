import networkx as nx
from biographs.lib.bpdb import residue_adjacency


def network(model, cutoff=5, weight=True):
    """Return the interaction network of a protein structure.

    The interaction network is defined by a distance cutoff.

    Parameters
    ----------
    model: Bio.PDB.model
        The protein structure.
    cutoff: float
        The distance cutoff defining an interaction between two nodes.
    weight: boolean
        True if atomic interactions are to be considered.
    """

    adjacency_dictionary = residue_adjacency(model, cutoff=cutoff,
                                            weight=weight)

    return nx.Graph(adjacency_dictionary)
