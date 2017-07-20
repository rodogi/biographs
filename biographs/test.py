#!/usr/local/bin/python
from pmolecule import Pmolecule

mol = Pmolecule('/Users/rdora/Documents/paper/pdbs/1efi.pdb')
print(mol.void_convex_hulls())
