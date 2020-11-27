from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

import molecule_viewer
import molecule_viewer.methods
from molecule_viewer.methods import ConformationMethod

class ChemicalMolecule:
    """Represents molecules from the ChEMBL database."""

    def __init__(self, chemblid, conformation_method=None):
        """Query ChEMBL for a molecule and produce a 3D conformer.

        Args:
            chemblid (str): The ChEMBL ID for the molecule
            conformation_method (molecule_viewer.methods.ConformationMethod): The method used to generate the 3D conformer

        """
        self.chemblid = chemblid
        self.conformation_method = conformation_method or molecule_viewer.methods._default_conformation_method
        
        print(f"Querying CHEMBL for molecule{self.chemblid}...")
        res = new_client.molecule.filter(chembl_id='CHEMBL729').only(["molecule_chembl_id", "molecule_structures"])
        if not len(res):
            raise ValueError(f"No result for CHEMBL ID {chemblid}")
        else:
            chembl_data = next(res)
        
        smiles = chembl_data["molecule_structures"]["canonical_smiles"]
        self.mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(self.mol)
        if self.conformation_method == ConformationMethod.MMFF94:
            print(f"Optimising molecule with {self.conformation_method.name}...")
            AllChem.MMFFOptimizeMolecule(self.mol)
        elif method == ConformationMethod.UFF:
            print(f"Optimising molecule with {self.conformation_method.name}...")
            AllChem.UFFOptimizeMolecule(self.mol)
        elif method != ConformationMethod.ETKDG:
            raise ValueError(f"Unknown method {self.conformation_method.name}...")

        self.mol = Chem.RemoveHs(self.mol)
        

    def show(self, name=None):
        """Display the 3D conformer of the molecule in PyMOL

        Args:
            name (str): The temporary filename and PyMOL label. Defaults to '{chemblid}-{conformation_method}'

        """
        name = name or f"{self.chemblid}-{self.conformation_method.name}"
        print(f"Showing molecule {name}...")
        filename = molecule_viewer._wd / f"{name}.mol"
        with open(filename, "w") as f:
            f.write(Chem.MolToMolBlock(self.mol))
        cmd.load(filename, name)

    def __repr__(self):
        return f"{self.__class__.__name__}(chemblid={self.chemblid}, conformation_method={self.conformation_method})"

