import code
import argparse
import pathlib
import molecule_viewer
import pymol

def init():
    print("Initialising PyMOL...")
    pymol.pymol_argv = ["pymol", "-Q"]
    pymol.finish_launching()

def quit():
    print("Quitting PyMOL...")
    pymol.cmd.quit()

def get_imports():
    from molecule_viewer.chemical_molecule import ChemicalMolecule
    from molecule_viewer.protein import Protein
    return dict(ChemicalMolecule=ChemicalMolecule, Protein=Protein)

def viewer(working_directory, scripts):
    if working_directory is not None:
        molecule_viewer._set_wd(working_directory)

    env = get_imports()
    for script in scripts:
        print(f"Executing script {script}...")
        with open(script, "r") as f:
            script_string = f.read()
            try: 
                exec(script_string, env)
            except Exception as ex:
                print(f"Error in script {script}")
                raise ex
    print("Entering interactive mode. Type `quit()` to exit (don't close PyMOL).")
    code.interact(banner="(Molecule Viewer Interactive Console)", exitmsg="Exiting interactive console...", local=env)

def main():
    parser = argparse.ArgumentParser(description="Visualise proteins and chemical molecules using PyMOL.")
    parser.add_argument("--working-directory", "-w", metavar="DIRECTORY",
                        help="Working directory where molecule files are stored.")
    parser.add_argument("scripts", help="Scripts to execute before entering interactive mode.", nargs="*")
    args = parser.parse_args()

    init()
    try:
        viewer(args.working_directory, args.scripts)
    finally:
        quit()