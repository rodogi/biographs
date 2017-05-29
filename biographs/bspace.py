"""Tools to deal with the void around residues in proteins
"""

from collections import defaultdict, deque
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from Bio.PDB import Selection
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
    atoms = [atom for atom in model.get_atoms()]
    DT = Delaunay([a.coord for a in atoms])  # Delaunay tessellation

    def void_residue(residue, *args):
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
        atom_indices = set([atoms.index(atom) for atom in residue])

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
                n_n = indptr[indices[k]:indices[k+1]]  # Neighbors of k

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
                    if ((n, k) in selec_tri) and ((n, j) in selec_tri):
                        simplex = set([v, k, j, n])
                        # Select tuples not containing only points in residue.
                        if len(simplex.intersection(atom_indices)) < 4:
                            cv = ConvexHull([atoms[a].coord for a in simplex])
                            void += cv.volume

        return void

    if (mu or sigma):
        dis_edges = []
        for simplex in DT.simplices:
            for i in range(3):
                for j in range(i+1, 4):
                    dis_edges.append(atoms[simplex[i]] - atoms[simplex[j]])

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


def void_ken_dill(model):
    """Return dict with the void of each residue in `model`.

    Atom radii are taken from [1] except for hidrogen taken from [2].

    Parameters
    ----------
    model: Bio.PDB.Model.Model
    The structure model of the protein

    Notes
    -----
    [1]: D. Flatow et al. (Volumes and surface areas: Geometries and scaling
    relationships between coarse grained and atomic structures).
    [2]: J.C Gaines et al. (Packing in protein cores).

    """
    atom_radii = {
        'C': 1.70,
        'F': 2.00,  # Iron
        'H': 1.02,
        'I': 1.98,
        'M': 2.00,  # Magnesium
        'N': 1.55,
        'O': 1.52,
        'S': 1.80}

    atoms = np.array([atom for atom in model.get_atoms()])
    delaunay_triangulation = Delaunay([atom.coord for atom in atoms])
    simplices = delaunay_triangulation.simplices
    neighbors = delaunay_triangulation.neighbors

    def _empty_triangles(simplex):
        """Return True if simplex is empty, False otherwise.

        """
        empty_triangles = []

        for i in range(2):
            for j in range(i+1, 3):
                for k in range(j+1, 4):
                    atom_i = atoms[simplex[i]]
                    atom_j = atoms[simplex[j]]
                    atom_k = atoms[simplex[k]]

                    try:
                        radius_i = atom_radii[atom_i.id[0]]
                    except KeyError:
                        raise Exception(
                            "Radius of atom {} is not defined".format(atom_i.id))
                    try:
                        radius_j = atom_radii[atom_j.id[0]]
                    except KeyError:
                        raise Exception(
                            "Radius of atom {} is not defined".format(atom_j.id))
                    try:
                        radius_k = atom_radii[atom_k.id[0]]
                    except KeyError:
                        raise Exception(
                            "Radius of atom {} is not defined".format(atom_k.id))

                    if (
                        radius_i + radius_j < atom_i - atom_j and
                            radius_i + radius_k < atom_i - atom_k and
                            radius_j + radius_k < atom_j - atom_k):
                        empty_triangles.append((i, j, k))

        if len(empty_triangles) == 4:
            return True

        return False

    def _depth_first_search(index):
        """Return list of connected component and checked indices.

        """
        checked_indices = []
        stack = deque([index])
        connected_component = []

        if not _empty_triangles(simplices[index]):
            checked_indices.append(index)

            return connected_component, checked_indices

        while stack:
            checked_indices.append(index)
            #empty_triangles = _empty_triangles(simplices[index])
            local_neighbors = neighbors[index]
            # Check if residue at the boundary
            if (-1) in local_neighbors:
                """
                # Where is the boundary
                where = np.where(local_neighbors == -1)[0]
                for _local_neighbor in where:
                    _triangle = tuple(_index for _index in range(4) if
                                      _index != _local_neighbor)
                    if _triangle in empty_triangles:
                        connected_component = []

                        # connected_component is set to zero if it is
                        # connected to a simplex with an empty triangle
                        # on the boundary.
                """
                connected_component = []
                return connected_component, checked_indices

            good_neighbors = []
            for n, neighbor in enumerate(local_neighbors):
                # Check if neighbor shares an empty triangle.
                """
                _triangle = tuple(_index for _index in range(4) if _index != n)
                if _triangle in empty_triangles and neighbor not in checked_indices:
                    good_neighbors.append(neighbor)
                """
                if _empty_triangles(simplices[neighbor]):
                    if neighbor not in checked_indices:
                        good_neighbors.append(neighbor)

            if not good_neighbors:
                connected_component.append(index)
                # Remove `index` of stack.
                stack.pop()
                # If still elements in stack then set index to its last element.
                if stack:
                    index = stack[-1]

            else:
                stack.extend(good_neighbors)
                index = stack[-1]

        return connected_component, checked_indices

    def _bounded_regions():
        """Return the bounded regions of the alpha-complex of the triangulation.

        """
        indices = set(range(len(simplices)))
        checked_indices = []
        connected_components = []

        while indices - set(checked_indices):
            # Take an element at random from the different between indices and
            # checked indices.
            remaining_set = indices - set(checked_indices)
            # Pick one of the remaining set.
            index = next(iter(remaining_set))
            connected_component, _checked_indices = _depth_first_search(index)
            if connected_component:
                connected_components.append(connected_component)
            checked_indices.extend(_checked_indices)

        return connected_components

    def _volume_overlap(conv_simplex, atoms_simplex):
        """Return the volume of the overlap between the sphere and the tetrahedron.

        The volume of the overlap between the tetrahedron and the sphere,
        V_{overlap} is given by the equation:

        V_{overlap} = \frac{2*r^3}{6} * [ -pi + phi_1 + phi_2 + phi_3 ],

        where phi_i is the dihedral angle between the planes intersecting at edge
        `i` of the convex hull. If `n_1` and `n_2` are the normal vectors of the
        planes intersecting at `edge 1`, then

        phi_1 = - [ n_1 \cdot n_2 ].

        Where \cdot represents the dot product of both vectors. Note: The minus
        sign in the formula is added taken into account the direction of the
        normal vectors found in the equation of the planes defined by the facets
        in `scipy.spatial.ConvexHull`.

        Parameters
        ----------
        conv_simplex : scipy.spatial.ConvexHull
        Convex hull of a tetrahedron belonging to the Delaunay tesselation of a set
        of atomic coordinates.

        atoms_simplex : list
        Atoms used for the computation of the convexhull. The order of the
        `atoms_simplex` is the order of `conv_simplex.vertices`.

        """

        # Equations of the planes defined by the facets of the convex hull are of
        # the form [a, b, c, d] where a*x + b*y + c*z + d = 0. Vector
        # n = (a, b, c) is the normal vector of the facet with direction to
        # the `outside` of the convex hull.
        equations = conv_simplex.equations
        # Tehtrahedra
        simplices = conv_simplex.simplices
        vertices = conv_simplex.vertices
        volume_overlaps = []

        for vertex in vertices:
            dihedral_angles = []
            # Adjacent facets of `vertex`
            adjacent_facets = np.where(simplices == vertex)[0]
            for first, second in [[0, 1], [1, 2], [2, 0]]:
                # Normal vectors
                normal_first = equations[adjacent_facets[first]][:-1]
                normal_second = equations[adjacent_facets[second]][:-1]
                # cos(phi) = - [ n_1 \cdot n_2 ]
                cosine_angle = - np.dot(normal_first, normal_second)
                dihedral_angle = np.arccos(cosine_angle)
                dihedral_angles.append(dihedral_angle)

            # Volume overlap
            # For `radius` take only the starting letter of atom.
            radius = atom_radii[atoms_simplex[vertex].id[0]]
            factor = (radius**3) / 3.
            volume_overlap = factor * (-np.pi + sum(dihedral_angles))
            volume_overlaps.append(volume_overlap)

        return volume_overlaps

    connected_components = _bounded_regions()
    void = defaultdict(int)
    # In General expect nothing of this void. Residues are well-packed.
    inner_void = defaultdict(int)

    for connected_component in connected_components:
        for simplex in simplices[connected_component]:
            atoms_simplex = atoms[simplex]
            residues_simplex = Selection.get_unique_parents(atoms_simplex)
            conv_simplex = ConvexHull([atom.coord for atom in atoms_simplex])
            volume_overlap = _volume_overlap(conv_simplex, atoms_simplex)
            simplex_void = conv_simplex.volume - sum(volume_overlap)
            if 1 < len(residues_simplex):
                for residue in residues_simplex:
                    void[label_residue(residue)] += simplex_void
            else:
                residue = label_residue(residues_simplex[0])
                inner_void[residue] += simplex_void

    return void, inner_void


def difference_convex_hulls(protein_name, chain='All',
                            min_seq=None, max_seq=None):
    """
    Return void by computing the difference in volume
    between the amino acid convex hull and a larger convex hull
    defined by at least len(amino_acid) points.
    Parameters
    ----------
    amino_acid: str, amino acid e.g. 'A311'.
    protein_name: str, name of the PDB file
    """
    D = {} #Final dictionary
    def get_set_points_protein(model, chain='All', min_seq=None, max_seq=None):
        """
        Return set of points corresponding to the atomic coordinates of a
        protein part.

        parameters
        ----------
        protein_name: str, name of the PDB file
        chain: str, list of str, name of interested chains
        min_seq: int, minimum sequence number
        max_seq: int, maximum sequence number
        """
        if type(chain) is str and chain!='All':
            chain = list(chain)
        if chain is 'All':
            chain = [c.id for c in model.child_list]
        if min_seq is None:
            min_seq = -float('inf')
        if max_seq is None:
            max_seq = float('inf')

        A = [a for a in Bio.PDB.Selection.unfold_entities(model, 'A') if
             a.parent.parent.id in chain and min_seq<=a.parent.id[1]<=max_seq]

        R = [r for r in Bio.PDB.Selection.unfold_entities(model, 'R') if
             r.parent.id in chain and min_seq<=r.id[1]<=max_seq]

        return [tuple(a.coord) for a in A], R

    Model = bpdb.pdb_model(protein_name)
    protein_points, R = get_set_points_protein(Model, chain=chain,
                                               min_seq=min_seq, max_seq=max_seq)
    for r in R:
        atoms = [tuple(a.coord) for a in r.child_list]
        P = set(protein_points) - set(atoms)
        cv = ConvexHull(atoms)
        simplices = cv.simplices
        vertices = cv.vertices
        eq = cv.equations
        N = [] #Points to extend convex hull
        for j, simplex in enumerate(simplices):
            for i, p_ in enumerate(atoms):
                if i not in simplex:
                    break

            sign = eq[j][0]*p_[0]+eq[j][1]*p_[1]+eq[j][2]*p_[2]+eq[j][3]
            #Separate all points on the right side of the plane:
            half_space = [np.array(p) for p in P if
                          sign*(p[0]*eq[j][0]+p[1]*eq[j][1]+p[2]*eq[j][2]
                                +eq[j][3])<0]

            if len(half_space) is 0:
                break

            u = np.array(atoms[simplex[1]]) - np.array(atoms[simplex[0]])
            v = np.array(atoms[simplex[2]]) - np.array(atoms[simplex[0]])
            n_tri = np.cross(u, v)

            n = eq[j][:3] #Unit normal vector
            in_triangle = []

            for hp in half_space:
                w = hp - np.array(atoms[simplex[0]])
                gamma = np.dot(np.cross(u, w), n_tri)/np.dot(n_tri,n_tri)
                beta = np.dot(np.cross(w, v), n_tri)/np.dot(n_tri,n_tri)
                alpha = 1. - gamma - beta
                #Check if all 3 barycentric coordinates of hp are < 0:
                if 0<=alpha<=1 and 0<=beta<=1 and 0<=gamma<=1:
                    #Select point and annotate its distance to plane:
                    in_triangle.append([hp, np.abs(np.dot(hp, n) + eq[j][3])])

            if len(in_triangle) is not 0:
                Spoint = min(in_triangle, key=lambda x: x[1])
                #Point kept only if distance to plane < 5 angstroms:
                if Spoint[1] <= 5.:
                    N.append(Spoint[0])
                else:
                    #Else, approach point until distance = 5 angstroms:
                    Spoint[0] = Spoint[0] - n*(Spoint[1]-5.)
                    #Make sure is at 5 angstroms from at least one point:
                    _p1 = np.array(atoms[simplex[0]])
                    _p2 = np.array(atoms[simplex[1]])
                    _p3 = np.array(atoms[simplex[2]])
                    _min_dis = min(np.linalg.norm(_p1-Spoint[0]),
                                   np.linalg.norm(_p2-Spoint[0]),
                                   np.linalg.norm(_p3-Spoint[0]))
                    while _min_dis > 5.:
                        Spoint[0] = Spoint[0] - n*(0.1)
                        _min_dis = min(np.linalg.norm(_p1-Spoint[0]),
                                       np.linalg.norm(_p2-Spoint[0]),
                                       np.linalg.norm(_p3-Spoint[0]))
                    N.append(Spoint[0])
        if len(N) is 0:
            print 'ZeroVoid: '+protein_name[-8:-4]+r.parent.id+str(r.id[1]) 
        p_vertices = [np.array(atoms[m]) for m in vertices]
        p_vertices = np.array(p_vertices + N)
        cv1 = ConvexHull(p_vertices)
        D[r.parent.id+str(r.id[1])] = [cv.volume, cv1.volume - cv.volume]

    return D
