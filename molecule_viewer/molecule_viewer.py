import pymol
from molecule_viewer.protein import Protein

def init():
    pymol.finish_launching()

def quit():
    pymol.cmd.quit()

def main():
    init()

    protease = Protein(pdbid="6LU7")
    protease.show()