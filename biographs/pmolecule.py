# class to deal with protein structures
# python 2
import bpdb
import bgraph
import bspace


class Pmolecule(object):

    def __init__(self, structure_file, water=False):
        self.model = bpdb.pdb_model(structure_file, water=water)
        self.path_to_file = structure_file

    def network(self, cutoff=5, weight=True):
        model = self.model

        return bgraph.network(model, cutoff=cutoff, weight=weight)

    def void(self, cutoff=5, mu=0, sigma=0):
        model = self.model

        return bspace.void_delaunay(model, cutoff=cutoff, mu=mu, sigma=sigma)

    def volume_delaunay(self):
        model = self.model

        return bspace.volume_delaunay(model)

    def volume_convex_hull(self):
        model = self.model

        return bspace.volume_convex_hull(model)

    def void_alpha_shape(self):
        model = self.model

        return bspace.void_ken_dill(model)

    pass
