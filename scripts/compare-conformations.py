from molecule_viewer.methods import ConformationMethod

lopinavir = ChemicalMolecule(chemblid="CHEMBL729")

methods = ["MMFF94"]
for method in ConformationMethod:
    lopinavir.show(method=method)