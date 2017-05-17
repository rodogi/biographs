# class to deal with protein structures
# python 2
from bpdb import pdb_model

class Pmolecule(object):

    def __init__(self, structure_file, water=False):
        self.model = pdb_model(structure_file, water=water)
        pass

    pass
