from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

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
                if hetero_flag: # empty strings evaluate to False
                                # the above conditional evaluates True
                                # if heter_flag indicates a water residue
                    chain.detach_child(residue.id)
            if not list(chain): # empty lists evaluate to False
                model.detach_child(chain.id)

    return model

