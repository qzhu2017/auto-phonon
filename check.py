from pymatgen.core.structure import Structure
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.io import read
import os, time
import warnings
import phonopy
warnings.filterwarnings("ignore")

"""
scripts to perform structure conformation
"""
def pymatgen2ase(struc):
    atoms = Atoms(symbols = struc.atomic_numbers, cell = struc.lattice.matrix)
    atoms.set_scaled_positions(struc.frac_coords)
    return atoms

def ase2pymatgen(struc):
    lattice = struc.get_cell()
    coordinates = struc.get_scaled_positions()
    species = struc.get_chemical_symbols()
    return Structure(lattice, species, coordinates)

def get_subdir(a_dir):
    return sorted([name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))])

cwd = os.getcwd()
count_good = 0
count_imag = 0
good_dirs = []
for dir in ['todo', 'jxli', 'jxliu']:
    os.chdir(dir)
    for dir_type in get_subdir('./'):
        os.chdir(dir_type)
        for dir_type1 in get_subdir('./'):
            if dir_type1.find('mp') > -1:
                os.chdir(dir_type1)
                if os.path.exists('FORCE_CONSTANTS'): 
                    if os.path.exists('band.png'):
                        if os.path.exists('IMAGINARY_FREQUENCY'):
                            count_imag += 1
                        else:
                            count_good += 1
                            good_dirs = dir + '/' + dir_type + '/' + 'dir_type1'
                    else:
                        struc = ase2pymatgen(read('POSCAR-unitcell', format='vasp'))
                        finder = SpacegroupAnalyzer(struc)
                        spg = finder.get_space_group_symbol()[0]
                        formula = struc.composition.get_reduced_formula_and_factor()[0]
                        if spg == 'P':
                            spg = 'auto'
                        with open('phonopy.conf', 'r') as f:
                            lines = f.readlines()
                            for line in lines:
                                if line.find('DIM') > -1:
                                    tmp = line.split('=')[-1]
                                    matrix=[eval(i) for i in tmp.split()]
                                    if len(matrix) == 3:
                                        supercell_matrix = np.diag(matrix)
                                    elif len(matrix) == 9:
                                        supercell_matrix = np.array(matrix).reshape([3,3])
                                    break
                        try:
                            phonon = phonopy.load(unitcell_filename = "POSCAR-unitcell",
                                                  force_constants_filename = "FORCE_CONSTANTS",
                                                  supercell_matrix = supercell_matrix,
                                                  primitive_matrix= spg,
                                                  )
                            phonon.run_mesh(mesh=[10,10,10])
                            phonon.run_total_dos()
                            phonon.plot_total_DOS().savefig('dos.png')
                            phonon.auto_band_structure(plot=True).savefig('band.png')
                            mesh_dict = phonon.get_mesh_dict()
                            frequencies = mesh_dict['frequencies']  
                            freqs = frequencies.ravel()
                            if len(freqs[freqs<-2e-2])>0:
                                print('IMAGINARY FREQUENCY exists in ', dir, dir_type, dir_type1, formula)
                                with open('IMAGINARY_FREQUENCY', 'w') as f:
                                    f.write('IMAGINARY FREQUENCY exists {:s} {:s}\n'.format(dir_type1, formula))
                                count_imag += 1
                            else: 
                                print('Good phonon in ', dir, dir_type, dir_type1, formula)
                                good_dirs = dir + '/' + dir_type + '/' + 'dir_type1'
                                count_good += 1
                        except RuntimeError: 
                            print('RuntimeError: Dynamical matrix has not yet built ', dir, dir_type, dir_type1, formula)
                os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir(cwd)
for dir in good_dirs:
    os.system('zip -r my.zip ' + dir + '/FORCE_CONSTANTS ' + dir + '/POSCAR-unitcell ' + dir + '/phonopy.conf')
print('------good phonons: ', count_good)
print('------imag phonons: ', count_imag)
