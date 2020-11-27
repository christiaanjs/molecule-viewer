from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

import molecule_viewer
from molecule_viewer.methods import ConformationMethod

class ChemicalMolecule:
    def __init__(self, chemblid):
        self.chemblid = chemblid
        print(f"Querying CHEMBL for molecule{self.chemblid}...")
        res = new_client.molecule.filter(chembl_id='CHEMBL729').only(["molecule_chembl_id", "molecule_structures"])
        if not len(res):
            raise ValueError(f"No result for CHEMBL ID {chemblid}")
        else:
            self.molecule = next(res)
        

    def show(self, name=None, method=ConformationMethod.MMFF94):
        smiles = self.molecule["molecule_structures"]["canonical_smiles"]
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol)
        if method == ConformationMethod.MMFF94:
            print(f"Optimising molecule with {method.name}...")
            AllChem.MMFFOptimizeMolecule(mol)
        elif method == ConformationMethod.UFF:
            print(f"Optimising molecule with {method.name}...")
            AllChem.UFFOptimizeMolecule(mol)
        elif method != ConformationMethod.ETKDG:
            raise ValueError(f"Unknown method {method.name}...")

        mol = Chem.RemoveHs(mol)
        name = name or f"{self.chemblid}-{method.name}"
        filename = molecule_viewer._wd / f"{name}.mol"
        with open(filename, "w") as f:
            f.write(Chem.MolToMolBlock(mol))
        cmd.load(filename, name)

