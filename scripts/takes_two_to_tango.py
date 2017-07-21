""" python tools to create the data used in `takes two to tango`

From a text file of PDB (Protein data bank) IDs, return a database.


The database is a comma separated value (csv) file containing for each line a
different residue belonging to one of the 252 proteins used in the article.

For each line we have:
Residue name | Type of residue | Degree (k) | Weight (w) | Atomic number (size)

Where Degree and Weight are network properties of the residue and the atomic
number is its number of atoms as given in the structural file.
"""

import shutil
import sys
import os
from glob import glob
import numpy as np
from Bio.PDB import PDBList
import networkx as nx
sys.path.append(os.path.realpath(__file__).rsplit('/', 2)[0])
import biographs as bg

current_path = os.path.realpath(__file__).rsplit('/', 1)[0]

# Create a temporal directory called `pdb'
pdbs_path = os.path.join(current_path, 'pdb')

if pdbs_path not in glob(os.path.join(current_path, '*')):
    os.mkdir(pdbs_path)

with open(os.path.join(current_path, 'pdbs.txt'), 'r') as f:
    pdbs = [pdb[:-1] for pdb in f]

pdbl = PDBList(obsolete_pdb=True)

pdbs = glob(os.path.join(pdbs_path, '*'))

if not pdbs:
    pdbl.download_pdb_files(pdbs, file_format='pdb', pdir=pdbs_path)

database = [['Residue name', 'Type of residue', 'Degree', 'Weight',
             'Atomic number']]

for pdb in pdbs:
    mol = bg.Pmolecule(pdb)
    net = mol.network()
    for residue in mol.model.get_residues():
        deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
        weight = nx.degree(net, residue.parent.id + str(residue.id[1]),
                           weight='weight')
        restype = residue.resname
        resname = pdb.rsplit('/', 1)[1][:-4] + residue.parent.id \
                + str(residue.id[1])
        size = len(residue)
        database.append([resname, restype, deg, weight, size])

database = np.array(database, dtype='S')
np.savetxt(os.path.join(current_path, 'database.csv'), database, delimiter=',',
           fmt='%s')

# Delete temporal directory `pdb'
#shutil.rmtree(pdbs_path, ignore_errors=True)
