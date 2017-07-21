""" python tools to create the data used in `takes two to tango`

From a text file of PDB (Protein data bank) IDs, return a database.


The database is a comma separated value (csv) file containing for each line a
different residue belonging to one of the 152 proteins used in the article.

For each line we have:
    Residue name | Type of residue | Degree | Weight | Atomic number

Where Degree and Weight are network properties of the residue and the atomic
number is its number of atoms as given in the structural file.
"""

from __future__ import absolute_import

from ..biographs import pmolecule

