from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

class ChemicalMolecule:
    def __init__(self, chemblid):
        self.chemblid = chemblid
        res = new_client.molecule.filter(chembl_id='CHEMBL729').only(["molecule_chembl_id", "molecule_structures"])
        if not len(res):
            raise ValueError(f"No result for CHEMBL ID {chemblid}")
        else:
            self.molecule = next(res)
        

    def show(self, name=None):
        smiles = self.molecule["molecule_structures"]["canonical_smiles"]
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)
        self.filename = f"{name or self.chemblid}.mol"
        with open(self.filename, "w") as f:
            f.write(Chem.MolToMolBlock(mol))
        cmd.load(self.filename)

