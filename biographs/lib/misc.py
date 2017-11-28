"""Miscelaneous functions"""

from __future__ import absolute_import

from biographs.classes.pmolecule import Pmolecule
from scipy.spatial import ConvexHull


def buriedness(pmol):
    """Return amino acids average distance to protein surface.
    
    The surface of the protein is modeled by its convex hull.

    Parameters
    ----------
    pmol : Pmolecule
        :pmolecule:Pmolecule: object
    """

    atoms = np.array([atom for atom in pmol.model.get_atoms()])
    surface = ConvexHull(atoms)
