# Automated Phonon Calculation

This is a backup repo for the script related to the following paper:

[Li J-X, Liu J-X, Baronett S.A., Li R-H, Wang L., Zhu Q. and Chen X-Q (2021). Computation and data driven discovery of topological phononic materials. Nature Communication, 12, 1204](https://www.nature.com/articles/s41467-021-21293-2)

## Setup

The following packages are required
- `ase` : the main script to call manage vasp calculation, ``pip install ase``
- `phonopy`: prepare vasp inputs, ``pip install phonopy``
- `Pymatgen`: to handle crystal symmetry, ``pip install pymatgen``
- [`vasprun`](https://github.com/qzhu2017/vasprun): to extract information from vasprun.xml, ``pip install vasprun-xml``.

In ``vasp.py``, change the followings based on your own setting,

```python
    # Setup VASP here
    VASP_CMD = '/home/x-qiangz/Github/VASP/vasp.5.4.4.pl2/bin/vasp_std > log'
    cmd = 'srun -n $SLURM_NTASKS --mpi=pmi2 ' + VASP_CMD
    os.environ['ASE_VASP_COMMAND'] = cmd                      # The actual command to launch VASP
    os.environ['VASP_PP_PATH'] = '/home/x-qiangz/Github/VASP' # VASP pseudopotential directory path
```


## Run
Assuming you have installed all the required packages, 
```
$ python vasp.py -f libs/POSCARs-A4-diamond
```

In the first run, the script will automatically perform geometry optimization and force calculation which are required to obtain the dynamical matrix file for phonopy.
If the calculation is complete, this will generate a folder according to the chemical formula (e.g., `Si8`) and a file called `FORCE_CONSTANTS` there.

If you call ```python vasp.py -f libs/POSCARs-A4-diamond``` one more time, it will call phonopy-api to compute phonon dispersion and DOS. Finally, it will generate `band.png` and `dos.png` for you to check. 


For your convenience, there is a slurm script called `myrun_vasp`. You must modify the setting based on your own environment.
