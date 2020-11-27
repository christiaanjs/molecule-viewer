from pymol import cmd

class Protein:
    """Represents proteins from PDB"""

    def __init__(self, pdbid):
        """
        Args:
            pdbid (str): The PDB ID of the protein
        """
        self.pdbid = pdbid

    def show(self):
        """Display the protein in PyMOL"""
        print(f"Showing protein {self.pdbid}...")
        cmd.fetch(self.pdbid)

    def __repr__(self):
        return f"{self.__class__.__name__}(pdbid={self.pdbid})"
