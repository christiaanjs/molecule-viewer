import code
import argparse
import pathlib
import molecule_viewer
import molecule_viewer.methods
import pymol

def init():
    print("Initialising PyMOL...")
    pymol.pymol_argv = ["pymol"]
    pymol.finish_launching()

def quit():
    print("Quitting PyMOL...")
    pymol.cmd.quit()

def get_imports():
    from molecule_viewer.chemical_molecule import ChemicalMolecule
    from molecule_viewer.protein import Protein
    return dict(ChemicalMolecule=ChemicalMolecule, Protein=Protein)

def viewer(working_directory, default_conformation_method, scripts):
    if working_directory is not None:
        molecule_viewer._set_wd(working_directory)

    if default_conformation_method is not None:
        molecule_viewer.methods._set_default_conformation_method(default_conformation_method)

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
    for method in molecule_viewer.methods.ConformationMethod:
        help_text = f"Use the {method.name} method for 3D conformation."
        if method == molecule_viewer.methods._default_conformation_method:
            help_text += " (default)"
        parser.add_argument(f"--{method.name}", action="store_const", const=method, dest="default_conformation_method", help=help_text)

    args = parser.parse_args()
    init()
    try:
        viewer(args.working_directory, args.default_conformation_method, args.scripts)
    finally:
        quit()