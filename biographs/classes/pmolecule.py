# class to deal with protein structures
# python 2

from __future__ import absolute_import
from biographs.lib.bpdb import pdb_model
from biographs.lib.bgraph import network
from biographs.lib.bspace import (void_delaunay, volume_delaunay,
    volume_convex_hull, void_ken_dill, void_convex_hulls)


class Pmolecule(object):

    def __init__(self, structure_file, water=False):
        self.model = pdb_model(structure_file, water=water)
        self.path_to_file = structure_file

    def network(self, cutoff=5, weight=True):
        model = self.model

        return network(model, cutoff=cutoff, weight=weight)

    def void(self, cutoff=5, mean=0, sigma=0):
        model = self.model

        return void_delaunay(model, cutoff=cutoff, mean=mean,
                                    sigma=sigma)

    def volume_delaunay(self):
        model = self.model

        return volume_delaunay(model)

    def volume_convex_hull(self):
        model = self.model

        return volume_convex_hull(model)

    def void_alpha_shape(self):
        model = self.model

        return void_ken_dill(model)

    def void_convex_hulls(self):
        model = self.model

        return void_convex_hulls(model)
