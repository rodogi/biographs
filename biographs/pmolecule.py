# class to deal with protein structures
# python 2
import bpdb
import bgraph

class Pmolecule(object):

    def __init__(self, structure_file, water=False):
        self.model = bpdb.pdb_model(structure_file, water=water)
        self.path_to_file = structure_file
        pass

    def network(self, cutoff=5, weight=True):
        path = self.path_to_file

        return bgraph.network(path, cutoff=cutoff, weight=weight)

    pass
