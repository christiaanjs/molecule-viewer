from pymol import cmd

class Protein:
    def __init__(self, pdbid):
        self.pdbid = pdbid

    def show(self):
        cmd.fetch(self.pdbid)
