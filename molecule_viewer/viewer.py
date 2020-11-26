import pymol
from molecule_viewer import Protein, ChemicalMolecule

def init():
    pymol.finish_launching()

def quit():
    pymol.cmd.quit()

def main():
    init()
    try:

        protease = Protein(pdbid="6LU7")
        protease.show()

        lopinavir = ChemicalMolecule(chemblid="CHEMBL729")
        lopinavir.show()

        nitroglycerin = ChemicalMolecule(chemblid="CHEMBL730")
        nitroglycerin.show()

        input()

    finally:
        quit()


if __name__ == "__main__":
    main()