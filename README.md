biographs: Amino acid networks in Python
========================================

What is biographs?
------------------

Biographs is a python module converting structural files formats like *pdb* or
*cif* into amino acid networks. It is written over the packages biopython and
networkx. This allows for a flexible and powerful handling of both: protein
networks and protein structures.

Notably biographs uses the great Bio.PDB module to construct the network given
a cutoff distance. Here, the network is a networkx.Graph object, which allows
high customization of the network and tons of tools for networks analysis.

To use biographs just type in a terminal:

		$ git clone https://github.com/rodogi/biographs.git

Here is an example to build the amino acid network of the cholera toxin (PDB
id: 1eei):

```python
import biographs as bg

molecule = bg.Pmolecule('1eei.pdb')
network = molecule.network()
network.nodes()[:10]
['G26', 'G27', 'G24', 'G25', 'G22', 'G23', 'G20', 'G21', 'G28', 'G29']
```

Required packages
-----------------

[biopython](http://biopython.org/wiki/Biopython)
Great to read and work with structure files such as pdb found in the protein
data bank. You can listen to the voices of three of the main contributors of
biopython in [this](https://www.podcastinit.com/biopython-with-peter-cock-wibowo-andrarto-and-tiago-antao-episode-125/>) podcast.

[networkx](https://networkx.github.io)
Has an extensive set of functions to deal with networks.

[scipy](https://www.scipy.org)
Mainstream package for scientific computation, here used mainly for numpy and
geometrical algorithms.
