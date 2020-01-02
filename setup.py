import setuptools

with open("README.org", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biographs",
    version="0.1",
    author="Rodrigo Dorantes-Gilardi",
    author_email="rodgdor@gmail.com",
    description="A package to work with protein structure networks",
    url="https://github.com/rodogi/biographs",
    keywords=["protein contact network", "protein graphs"],
    packages=["biographs", "biographs.classes", "biographs.lib"],
    install_requires=['networkx', 'biopython'],
    long_description=long_description,
    python_requires='>=3.5')
