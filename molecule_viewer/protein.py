from pymol import cmd

class Protein:
    def __init__(self, pdbid):
        self.pdbid = pdbid

    def show(self):
        print(f"Showing protein {self.pdbid}...")
        cmd.fetch(self.pdbid)
