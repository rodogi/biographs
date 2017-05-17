from Bio.PDB import *
from collections import Counter

def pdb_model(structure_file, water=False):
    protein_name, file_format = structure_file.rsplit('.', 1)
    accepted_formats = ['cif', 'pdb']
    parsers = [MMCIFParser, PDBParser]

    try:
        parser = parsers[accepted_formats.index(file_format)]
        parser = parser(QUIET=True)
    except ValueError:
        raise ValueError("Accepted structure files are: {}".format(
            accepted_formats))

    structure = parser.get_structure(protein_name, structure_file)
    model = structure[0]

    if not water:
        for chain in model.get_chains():
            for residue in list(chain):
                hetero_flag = residue.id[0].strip()
                # Empty strings evaluate to False.  Therefore hetero_flag
                # returns False if the residue is not a water molecule.
                if hetero_flag:
                    chain.detach_child(residue.id)
            if not list(chain):
                model.detach_child(chain.id)

    return model

def label_residue(residue):
    """ Return a string with the label of the biopython [1] residue object.

    The label of the residue is the following:
        Chain + Position + Type
    Where Type is the 3 letter amino acid format.

    Parameters
    ----------
    residue: Bio.PDB.Residue.Residue
        The residue to be labeled.

    Notes
    -----
    1. http://biopython.org/wiki/Biopython
    """
    name = residue.resname
    position = str(residue.id[1])
    chain = residue.parent.id

    return chain+position+name

def residue_adjacency(structure_file, cutoff=5, weight=True):
    """Return residue adjacency dictionary defined by cutoff distance.

    Parameters
    ----------
    structure_file: string
        File with atomic coordinates.
    cutoff: int or float
        Distance cutoff defining links between atoms.  Two atoms are adjacent
        if their distance is less than the given cutoff.

    See Also
    --------
    pdb_model

    """

    model = pdb_model(structure_file)
    atoms = Selection.unfold_entities(model, 'A')

    ns = NeighborSearch(atoms)

    atomic_adjacency = {atom: set(ns.search(atom.coord, cutoff))
                          - set([atom]) for atom in atoms}
    residue_adjacency = {}

    # Create residue adjacency dictionary with string format, see
    # label_residue.
    for atom, neighbors in atomic_adjacency.items():
        residue = label_residue(atom.get_parent())
        residue_adjacency.setdefault(residue, [])

        # Only different residues are connected by an edge (No loops).
        not_in_residue = []
        for neighbor in atomic_adjacency[atom]:
            n_parent = label_residue(neighbor.get_parent())
            if n_parent is not residue:
                not_in_residue.append(n_parent)

        residue_adjacency[residue].extend(not_in_residue)


    if weight:
        # Make new dictionary mapping each residue to its neighbors taking 
        # into account the weight.
        weighted_adjacency = {}
        for residue in residue_adjacency:
            counter = Counter(residue_adjacency[residue])
            weighted_adjacency[residue] = {
                neighbor: {'weight': counter[neighbor]}
                               for neighbor in counter}
        return weighted_adjacency

    else:
        return residue_adjacency

