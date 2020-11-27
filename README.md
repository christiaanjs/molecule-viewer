# Molecule Viewer

Visualise proteins from [PDB](https://www.rcsb.org/) and molecules from [ChEMBL](https://www.ebi.ac.uk/chembl/) using [PyMOL](https://pymol.org/).

## Dependencies

* Python 3.7
* `pip`
* [PyMOL 2](https://pymol.org/2/)
* [RDKit](https://www.rdkit.org/) for 3D conformer generation
* [`chembl_webresource_client`](https://github.com/chembl/chembl_webresource_client)

## Installation

1. *(optional)* Install dependencies using [`conda`](https://docs.conda.io/en/latest/):
    * `conda env create -f environment.yml`
    * `conda activate molecule-viewer`.
2. `pip install .`

## Usage

If you've installed dependencies with `conda`, you'll need to run `conda activate molecule-viewer` first.

```
molecule-viewer [-h] [--working-directory DIRECTORY] [--MMFF94] [--UFF]
                       [--ETKDG]
                       [scripts [scripts ...]]

positional arguments:
  scripts               Scripts to execute before entering interactive mode.

optional arguments:
  -h, --help            show this help message and exit
  --working-directory DIRECTORY, -w DIRECTORY
                        Working directory where molecule files are stored.
  --MMFF94              Use the MMFF94 method for 3D conformer generation. (default)
  --UFF                 Use the UFF method for 3D conformer generation.
  --ETKDG               Use the ETKDG method for 3D conformer generation.
```

The `molecule-viewer` command will load PyMOL, optionally run some scripts, and enter an interactive Python console.
The [`ChemicalMolecule`](molecule_viewer/chemical_molecule.py) and [`Protein`](molecule_viewer/protein.py) classes, as well as the [`ConformationMethod` enum](molecule_viewer/methods.py) are available in this environment (rather than the PyMOL console).
To exit, the `quit()` at the interactive prompt, rather than exiting the PyMOL window.

See [molecule_viewer/methods.py](molecule_viewer/methods.py) for documentation on 3D conformer generation methods for `ChemicalMolecule`s. This method can also be supplied in the call to `ChemicalMolecule.show`

To run the demo, run:
```
molecule-viewer scripts/demo.py
```
(click 'Zoom' and 'Orient' in PyMOL to see both molecules)