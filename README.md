# Automated Phonon Calculation

This is a backup repo for the script related to the following paper:

[Li J-X, Liu J-X, Baronett S.A., Li R-H, Wang L., Zhu Q. and Chen X-Q (2021). Computation and data driven discovery of topological phononic materials. Nature Communication, 12, 1204](https://www.nature.com/articles/s41467-021-21293-2)

## Setup

The following packages are required
- `ase` : the main script to call manage vasp calculation
- `phonopy`: prepare vasp inputs
- `Pymatgen`: to handle crystal symmetry
- [`vasprun`](https://github.com/qzhu2017/vasprun): to extract information from vasprun.xml 

For `ase`, one needs to set up the vasp script and pseudopotential path following [this link](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html)

## Run
Assuming you have installed all the required packages, 
```
$ python vasp.py -f todo/A/POSCARs-A4-diamond_todo
```

the script will automatically perform geometry optimization and force calculation which are required to obtain the dynamical matrix file for phonopy.
