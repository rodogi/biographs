# class to deal with protein structures
# python 2

from __future__ import absolute_import
from biographs.lib.bpdb import pdb_model
from biographs.lib.bgraph import network
from biographs.lib.bspace import (void_delaunay, volume_delaunay,
    volume_convex_hull, void_ken_dill, void_convex_hulls)


class Pmolecule(object):
    """Create a Pmolecule object.

    The Pmolecule calls a number of methods for the analysis of protein
    structure. This includes the contruction of the interaction network of the
    protein.

    Parameters
    ----------
    structure_file = str
        The path to the structure file of the targeted protein. Three
        structure-file formats are accepted: `pdb', `cif', and `ent'.
    water: boolean, default is False
        If false, water molecules are ignored.

    Attributes
    ----------
    model: Bio.PDB.model
        The structural model of the structure. See www.biopython.org.
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

    def void(self, cutoff=5, mean=0, sigma=0):
        """Return dictionary with void of each residue

        Parameters
        ----------
            model : Bio.PDB.model
                Structure model of the protein.

            cutoff  : int or float, optional
                Upper bound for distance between nighbors in a
                tetrahedron.

            mean : int, optional

            sigma: int, optional
                If `cutoff` is set by the mean and/or standard deviation of the
                distribution of the edge lengths in the triangulation, the
                formula:: mean * mean + sigma * standard deviation is used.
                Typically if mean != 0, then mean == 1.

        """

        return void_delaunay(self.model, cutoff=cutoff, mean=mean,
                                    sigma=sigma)

    def volume_delaunay(self):
        """Returns dictionary containing volume for each residue in `model`.

        Parameters
        ----------
        model: Bio.PDB.Model.Model
            Model of the protein structure.

        """

        return volume_delaunay(self.model)

    def volume_convex_hull(self):
        """Return dictionary with the volume of each residue in model.

        The volume of each residue is equal to the volume of the convex hull of
        its atoms.

        Parameters
        ----------
        model: Bio.PDB.Model.Model
            Model for the structure of the protein.

        """

        return volume_convex_hull(self.model)

    def void_alpha_shape(self):
        """Return dict with the void of each residue in `model`.

        Atom radii are taken from [1] except for hidrogen taken from [2].

        Parameters
        ----------
        model: Bio.PDB.Model.Model
        The structure model of the protein

        Notes
        -----
        [1]: D. Flatow et al. (Volumes and surface areas: Geometries and
            scaling relationships between coarse grained and atomic
            structures).
        [2]: J.C Gaines et al. (Packing in protein cores).
        """

        return void_ken_dill(self.model)

    def void_convex_hulls(self):
        """ Return void around each residue using convex hulls.

        For each residue r, the void is defined as the difference in volume
        between the convex hull of r and a larger convex hull noted r'.

        Parameters
        ----------
        model : Bio.PDB.Model.Model
        The model of the protein structure defined by its atomic coordinates.

        """

        return void_convex_hulls(self.model)
