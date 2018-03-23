# class to deal with protein structures
# python 2

from __future__ import absolute_import
from biographs.lib.bpdb import pdb_model
from biographs.lib.bgraph import network


class Pmolecule(object):
    """Create a Pmolecule object.

    The Pmolecule calls a number of methods for the analysis of protein
    structure. This includes the contruction of the interaction network of the
    protein.

    Parameters
    ----------
    structure_file : str
        The path to the structure file of the targeted protein. Three
        structure-file formats are supported: `pdb', `cif', and `ent'.
    water : boolean, (default: False)
        If False, water molecules are removed.

    Attributes
    ----------
    model: Bio.PDB.model
        The structural model of the structure. See www.biopython.org.
    network: networkx:Graph
        The amino acid network of the protein based on a distance
        cutoff (default=5 angs). See :Pmolecule:network: for more info.
    path_to_file: str
        The path to the structural file used to instantiate the class.
    """

    def __init__(self, structure_file, water=False):
        self.model = pdb_model(structure_file, water=water)
        self.path_to_file = structure_file

    def network(self, cutoff=5, weight=True):
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

        return network(self.model, cutoff=cutoff, weight=weight)
