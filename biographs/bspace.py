from Bio.PDB import Selection
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from bpdb import label_residue

def void_delaunay(model, cutoff=5, mu=0, sigma=0):
    """Return dictionary with void of each residue

    Parameters
    ----------
        model : Bio.PDB.model
            Structure model of the protein.

        cutoff  : int or float, optional
            Upper bound for distance between nighbors in a
            tetrahedron.

        mu : int, optional

        sigma: int, optional
            If `cutoff` is set by the mean and/or standard deviation of the
            distribution of the edge lengths in the triangulation, the
            formula:: mu * mean + sigma * standard deviation is used. Typically
            if mu != 0, then mu == 1.

    """

    def void_residue(residue):
        """Return the void of the residue.

        Variable `atoms` containing a list of atoms of the protein must exist
        before calling this function.  Moreover, variables indices and indptr,
        created with the scipy.spatial.Delaunay method vertex_neighbor_vertices
        must exist before calling this function. Finally, variable `cutoff`
        must exist before calling this function.

        Parameters
        ----------
        residue: Bio.PDB.Residue.Residue
            Target residue.

        """

        void = 0
        atom_indices = set([atoms.index(atom) for atom in residue.get_atom()])

        # v represents a point (atom) of residue
        for v in atom_indices:
            selec_edges = []
            # Loop over the neighbors of v and select only neighbors n such
            # that dist(v, n) <= cutoff.
            for n in indptr[indices[v]:indices[v+1]]:
                if atoms[v] - atoms[n] <= cutoff:
                    selec_edges.append(n)

            selec_tri = []
            # For each pair of neighbors of v, (k, j), select it if k and j
            # are neighbors and dist(k, j) <= cutoff.
            for i_n, k in enumerate(selec_edges[:-1]):
                n_n = indptr[indices[k]:indices[k+1]] # Neighbors of k

                for j in selec_edges[i_n+1:]:
                    if atoms[k] - atoms[j] <= cutoff:
                        if j in n_n:
                            selec_tri.append((k, j))

            # select_tri contains all 3-tuple of pair-wise neighbors such that
            # at least one point in the tuple is a point (atom) of residue.
            # For each such tuple, add a point of the triangulation, if that
            # point is itself a neighbor of each point in the tuple and if
            # the resulting 4-tuple does not contain only points of residue.
            for k, j in selec_tri:
                for n in selec_edges:
                    # Check if all four are neighbors.
                    if ((n,k) in selec_tri) and ((n,j) in selec_tri):
                        simplex = set([v,k,j,n])
                        # Select tuples not containing only points in residue.
                        if len(simplex.intersection(atom_indices)) < 4:
                            cv = ConvexHull([atoms[a].coord for a in simplex])
                            void += cv.volume

        return void

    atoms = [atom for atom in model.get_atoms()]
    DT = Delaunay([a.coord for a in atoms]) # Delaunay tessellation

    if (mu or sigma):
        dis_edges = [atoms[simplex[i]] - atoms[simplex[j]] for simplex in DT.simplices
                     for i in range(3) for j in range(i+1, 4)]
        cutoff = mu * np.mean(dis_edges) + sigma * np.std(dis_edges)
        del dis_edges

    indices, indptr = DT.vertex_neighbor_vertices

    Void = {label_residue(residue): void_residue(residue) for residue in
            model.get_residues()}

    return Void

def volume_delaunay(model):
    """Returns dictionary containing volume for each residue in `model`.

    Parameters
    ----------
    model: Bio.PDB.Model.Model
        Model of the protein structure.

    """

    volume_dict = {}
    atoms = np.array([atom for atom in model.get_atoms()])
    DT = Delaunay([atom.coord for atom in atoms])

    for simplex in DT.simplices:
        parent_residues = Selection.get_unique_parents(atoms[simplex])
        # Simplex is taken into account only if is totally contained in one
        # residue.
        if len(parent_residues) is 1:
            unique_parent = label_residue(parent_residues[0])
            cv_simplex = ConvexHull([atom.coord for atom in atoms[simplex]])
            volume_dict.setdefault(unique_parent, 0)
            volume_dict[unique_parent] += cv_simplex.volume

    return volume_dict

def volume_convex_hull(model):
    """Return dictionary with the volume of each residue in model.

    The volume of each residue is equal to the volume of the convex hull of
    its atoms.

    Parameters
    ----------
    model: Bio.PDB.Model.Model
        Model for the structure of the protein.

    """

    volume_dict = {}

    for residue in model.get_residues():
        conv_residue = ConvexHull([atom.coord for atom in residue])
        volume_dict[label_residue(residue)] = conv_residue.volume

    return volume_dict
