"""Tools to deal with the void around residues in proteins
"""

from collections import defaultdict, deque
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from Bio.PDB import Selection
from biographs.lib.bpdb import label_residue


def void_delaunay(model, cutoff=5, mean=0, sigma=0):
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
    atoms = [atom for atom in model.get_atoms()]
    delaunay = Delaunay([a.coord for a in atoms])  # Delaunay tessellation
    indices, indptr = delaunay.vertex_neighbor_vertices

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
        atom_indices = set([atoms.index(atom) for atom in residue])

        # point1 represents a point (atom) of residue
        for point1 in atom_indices:
            edge = []
            # Loop over the neighbors of point1 and select only neighbors
            # point2 such that dist(point1, point2) <= cutoff.
            for point2 in indptr[indices[point1]:indices[point1+1]]:
                if atoms[point1] - atoms[point2] <= cutoff:
                    edge.append(point2)

            triangle = []
            # For each pair of neighbors of point1, (point3, point4), select
            # it if point3 and point4 are neighbors and dist(point3, point4)
            # less or equal than `cutoff`.
            for i_n, point3 in enumerate(edge[:-1]):
                # Neighbors of point3
                n_n = indptr[indices[point3]:indices[point3+1]]

                for point4 in edge[i_n+1:]:
                    if atoms[point3] - atoms[point4] <= cutoff:
                        if point4 in n_n:
                            triangle.append((point3, point4))

            # `triangle` contains all 3-tuple of pair-wise neighbors such that
            # at least one point in the tuple is a point (atom) of residue.
            # For each such tuple, add a point of the triangulation, if that
            # point is itself a neighbor of each point in the tuple and if
            # the resulting 4-tuple does not contain only points of residue.
            for point3, point4 in triangle:
                for point2 in edge:
                    # Check if all four are neighbors.
                    if (point2, point3) in triangle and \
                            (point2, point4) in triangle:
                        simplex = set([point1, point2, point3, point4])
                        # Select tuples not containing only points in residue.
                        if len(simplex.intersection(atom_indices)) < 4:
                            conv_simplex = ConvexHull([atoms[a].coord for a in
                                                       simplex])
                            void += conv_simplex.volume

        return void

    if mean or sigma:
        dis_edges = []
        for simplex in delaunay.simplices:
            for i in range(3):
                for j in range(i+1, 4):
                    dis_edges.append(atoms[simplex[i]] - atoms[simplex[j]])

        cutoff = mean * np.mean(dis_edges) + sigma * np.std(dis_edges)
        del dis_edges

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
    delaunay = Delaunay([atom.coord for atom in atoms])

    for simplex in delaunay.simplices:
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
                            "Radius of atom {} is not defined".format(
                                atom_i.id))
                    try:
                        radius_j = atom_radii[atom_j.id[0]]
                    except KeyError:
                        raise Exception(
                            "Radius of atom {} is not defined".format(
                                atom_j.id))
                    try:
                        radius_k = atom_radii[atom_k.id[0]]
                    except KeyError:
                        raise Exception(
                            "Radius of atom {} is not defined".format(
                                atom_k.id))

                    if (
                            radius_i + radius_j < atom_i - atom_j and
                            radius_i + radius_k < atom_i - atom_k and
                            radius_j + radius_k < atom_j - atom_k):
                        empty_triangles.append((i, j, k))

        # Simplex is open if for all 6 edges in the tetrahedron, the distance
        # between the endpoints is greater than the sum of their Van der Waals
        # radii.
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
            local_neighbors = neighbors[index]
            # Check if residue at the boundary
            if -1 in local_neighbors:
                connected_component = []

                return connected_component, checked_indices

            good_neighbors = []
            for neighbor in local_neighbors:
                # Check if neighbor is an empty triangle and hasn't been
                # checked yet.
                if _empty_triangles(simplices[neighbor]):
                    if neighbor not in checked_indices:
                        good_neighbors.append(neighbor)

            if not good_neighbors:
                connected_component.append(index)
                # Remove `index` of stack as simplex has no `open` neighbors
                # left to check.
                stack.pop()
                # If still simplices in stack then set index to its last
                # simplex.
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

        V_{overlap} = frac{2*r^3}{6} * [ -pi + phi_1 + phi_2 + phi_3 ],

        where phi_i is the dihedral angle between the planes intersecting at edge
        `i` of the convex hull. If `n_1` and `n_2` are the normal vectors of the
        planes intersecting at `edge 1`, then

        phi_1 = - [ n_1 dot n_2 ].

        Where `dot` represents the dot product of both vectors. Note: The minus
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
            if len(residues_simplex) > 1:
                for residue in residues_simplex:
                    void[label_residue(residue)] += simplex_void
            else:
                residue = label_residue(residues_simplex[0])
                inner_void[residue] += simplex_void

    return void, inner_void


def void_convex_hulls(model):
    """
    Return dictionary of with void around each residue using convex hulls.

    For each residue r, the void is defined as the difference in volume between
    the convex hull of r and a larger convex hull noted conv(r_v).

    Parameters
    ----------
    model : Bio.PDB.Model.Model
    The model of the protein structure defined by its atomic coordinates.

    """

    void = {}
    atomic_coordinates = [tuple(atom.coord) for atom in model.get_atoms()]

    for residue in model.get_residues():
        atoms = [tuple(atom.coord) for atom in residue]
        # Set of atomic coordinates NOT in residue
        not_in_residue = set(atomic_coordinates) - set(atoms)
        atoms = [atom.coord for atom in residue]

        conv_residue = ConvexHull(atoms)
        simplices = conv_residue.simplices
        vertices = conv_residue.vertices
        equations = conv_residue.equations

        r_prime = []  # Points `r_v' for extend convex hull
        for index, simplex in enumerate(simplices):
            # Take only points on top of the simplex
            points_on_top = [np.array(point) for point in not_in_residue if
                             equations[index][:-1].dot(point)
                             + equations[index][-1] > 0]

            if not points_on_top:
                break

            vect_u = atoms[simplex[1]] - atoms[simplex[0]]
            vect_v = atoms[simplex[2]] - atoms[simplex[0]]
            normal = np.cross(vect_u, vect_v)

            # Unit normal vector
            unit_normal = equations[index][:3]
            points_in_triangle = []

            for point in points_on_top:
                vect_w = point - atoms[simplex[0]]
                gamma = np.dot(np.cross(vect_u, vect_w),
                               normal) / normal.dot(normal)

                beta = np.dot(np.cross(vect_w, vect_v),
                              normal) / normal.dot(normal)

                alpha = 1 - gamma - beta

                # If the barycentric coordinates of a point are all three
                # within the 0,1 range; then the projection of the point lies
                # inside the triangle.
                if 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1:
                    # Select point and annotate its distance to plane.
                    points_in_triangle.append([point, point.dot(unit_normal)
                                               + equations[index][3]])

            if points_in_triangle:
                # Pick up the point at minimal distance from the triangle.
                selected_point = min(points_in_triangle, key=lambda x: x[1])
                # Point kept only if the distance to plane is less or equal
                # than 5 angstroms.
                if selected_point[1] <= 5:
                    r_prime.append(selected_point[0])
                else:
                    # Approach point until its distance from one vertex of the
                    # triangle is at distance is equal to 5 angstroms.
                    selected_point[0] = selected_point[0] - unit_normal \
                                         * (selected_point[1] - 5)

                    vertex_1 = atoms[simplex[0]]
                    vertex_2 = atoms[simplex[1]]
                    vertex_3 = atoms[simplex[2]]

                    # Select the vertex of the triangle at minimal distance
                    # from `selected_point`. Move `selected_point` closer to
                    # that vertex by 0.1 angstroms until `selected_point` is at
                    # distance less than 5 angstroms from the vertex.

                    min_dis = float("inf")  # Dummy value extremely big

                    point_min_dis = 0
                    for vertex in [vertex_1, vertex_2, vertex_3]:
                        distance = np.linalg.norm(vertex - selected_point[0])
                        if distance <= min_dis:
                            min_dis = distance
                            point_min_dis = vertex

                    while min_dis > 5.:
                        selected_point[0] = selected_point[0] - \
                                            unit_normal * 0.1
                        min_dis = np.linalg.norm(point_min_dis
                                                 - selected_point[0])

                    r_prime.append(selected_point[0])

        atoms_vertices = [atoms[m] for m in vertices]
        larger_convex_hull_vertices = atoms_vertices + r_prime
        larger_convex_hull = ConvexHull(larger_convex_hull_vertices)
        void[label_residue(residue)] = larger_convex_hull.volume \
            - conv_residue.volume

    return void
