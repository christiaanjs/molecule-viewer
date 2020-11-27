# Get the crystal structure of the main COVID-19 protease
# (See https://www.rcsb.org/structure/6LU7)
protease = Protein(pdbid="6LU7")
# Get the 3D structure of lopinavir
# (See https://www.ebi.ac.uk/chembl/compound report card/CHEMBL729)
lopinavir = ChemicalMolecule(chemblid="CHEMBL729")
# Interactively showing the molecules
protease.show()
lopinavir.show()