from pymatgen.core.structure import Structure
from vasprun.vasprun import vasprun
import numpy as np
from optparse import OptionParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase import Atoms
from ase.io import read
from ase.calculators.vasp import Vasp
import os, time
import warnings
warnings.filterwarnings("ignore")

def Read_POSCARS(filename):
    """
    Read POSCARS in a single file and convert it to Pymatgen structure object
    """
    Struc = []
    ids = []
    with open(filename, 'rb') as f:
        input_content = f.readlines()

    POSCAR_content = []
    N_atom = 0
    for str1 in input_content:
        tmp = str(str1, 'utf-8')
        POSCAR_content.append(tmp)
        if len(POSCAR_content) == 1:
            print(tmp)
            ids.append(tmp.split('|')[0])
        elif len(POSCAR_content) == 7:
            N_atom = sum(int(f) for f in str1.split())
        elif len(POSCAR_content) == 8+N_atom:
            pos_str = ''.join(POSCAR_content)
            try:
                p = Poscar.from_string(pos_str)
                Struc.append(p.structure)
            except:
                print('strucuture is wrong', pos_str)
                raise
            POSCAR_content = []
    return Struc, ids

def parse_ENMAX(filename='POTCAR'):
    enmax = []
    with open(filename, 'r') as f:
        contents = f.readlines()
        for line in contents:
            if line.find('ENMAX') > 0:
                tmp = line.split(';')[0]
                enmax.append(float(tmp.split('=')[-1]))
    return 1.3*max(enmax)

def prepare_vasp(struc, gap=0):
    with open('INCAR', 'w') as f1:
        ENCUT = parse_ENMAX()
        f1.write('PREC = Accurate\n')
        f1.write('ENCUT = {:.2f}\n'.format(ENCUT))
        f1.write('IBRION = 8\n')
        f1.write('EDIFF = 1.0e-08\n')
        f1.write('IALGO = 38\n')
        f1.write('ISMEAR = 0; SIGMA = 0.1\n')
        f1.write('LREAL = .FALSE.\n')
        f1.write('ADDGRID = .TRUE.\n')
        f1.write('LWAVE = .FALSE.\n')
        f1.write('LCHARG = .FALSE.\n')

    with open('KPOINTS', 'w') as f2:
        f2.write('Automatic mesh\n')
        f2.write('0\n')
        f2.write('Gamma\n')
        if gap > 0.001:
            f2.write('1     1     1\n')
        else:
            cell = struc.get_cell()
            a = np.linalg.norm(cell[0,:])
            b = np.linalg.norm(cell[1,:])
            c = np.linalg.norm(cell[2,:])
            kpoints = []
            for x in [a, b, c]:
                if x > 15:
                    kpoints.append(1)
                else:
                    kpoints.append(int(round(24/x)))
            f2.write('{:2d} {:2d} {:2d}\n'.format(kpoints[0], kpoints[1], kpoints[2]))
        f2.write('0.000 0.000 0.000\n')
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

def symmetrize_cell(struc, mode='C'):
    """
    symmetrize structure from pymatgen, and return the struc in conventional/primitive setting
    Args:
    struc: ase type
    mode: output conventional or primitive cell
    """
    P_struc = ase2pymatgen(struc)
    finder = SpacegroupAnalyzer(P_struc,symprec=0.05,angle_tolerance=5)
    if mode == 'C':
        P_struc = finder.get_conventional_standard_structure()
    else:
        P_struc = finder.get_primitive_standard_structure()

    return pymatgen2ase(P_struc)

def good_lattice(struc):
    maxvec = 25.0
    minvec = 2.5
    maxangle = 150
    minangle = 30
    para = struc.get_cell_lengths_and_angles()
    if (max(para[:3])<maxvec) and (max(para[3:])<maxangle) and (min(para[3:])>minangle):
        return True
    else:
        return False


"""
A script to perform multistages vasp calculation
"""

def set_vasp(level=0, encut=None, setup=None):
    default0 = {'xc': 'pbe',
            'npar': 8,
            'kgamma': True,
            'lcharg': False,
            'lwave': False,
            'setups': setup,
            }
    if encut is not None:
        default0['encut'] = encut

    if level==1: #pre-relax

        default1 = {'prec': 'accurate',
                'ibrion': 2,
                'kspacing': 0.15,
                'isif': 3,
                'ediff': 1e-4,
                'nsw': 50,
                }
    elif level==2: #relax before phonon
        default1 = {'prec': 'accurate',
                'ibrion': 2,
                'kspacing': 0.15,
                'isif': 3,
                'ediff': 1e-8,
                'ediffg': -1e-8,
                'nsw': 20,
                'lreal': False,
                'addgrid': True,
                }
    dict_vasp = dict(default0, **default1)
    return Vasp(**dict_vasp)

def read_OUTCAR(path='OUTCAR'):
    """read time and ncores info from OUTCAR"""
    time = 0
    ncore = 0
    for line in open(path, 'r'):
        if line.rfind('running on  ') > -1:
            ncore = int(line.split()[2])
        elif line.rfind('Elapsed time ') > -1:
            time = float(line.split(':')[-1])

    return time, ncore

def single_optimize(struc, level, encut=None, mode='C'):
    """single optmization"""
    struc = symmetrize_cell(struc, mode)
    struc.set_calculator(set_vasp(level, encut))
    energy = struc.get_potential_energy()
    print(energy)
    time, ncore = read_OUTCAR()
    struc = read('CONTCAR',format='vasp')
    symbols = struc.get_chemical_symbols()
    symbols_new = []
    for symbol in symbols:
        if symbol == 'X':
            symbols_new.append('Xe')
        else:
            symbols_new.append(symbol)

    struc.set_chemical_symbols(symbols_new)
    return struc, energy, time

parser = OptionParser()
parser.add_option("-f",  "--file", dest="posfile",
                  help="by filename in POSCAR format")
(options, args) = parser.parse_args()

strucs, ids = Read_POSCARS(options.posfile)
cwd = os.getcwd()
cmd = '/public/software/openmpi-1.8.8/Install/bin/mpirun -np 16 /public/sourcecode/vasp.5.4.4/vasp.5.4.4/bin/vasp_std '

for mid, struc_pmg in zip(ids, strucs):
    formula = struc_pmg.composition.get_reduced_formula_and_factor()[0]
    dir0 = mid
    if not os.path.exists(dir0):
        os.makedirs(dir0)
    os.chdir(dir0)
    print(dir0)
    if os.path.exists('FORCE_CONSTANTS'): # and not os.path.exists('band.png'):
        import phonopy
        struc = read('POSCAR-unitcell', format='vasp')
        finder = SpacegroupAnalyzer(ase2pymatgen(struc))
        spg = finder.get_space_group_symbol()[0]
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
                print('IMAGINARY FREQUENCY exists in ', dir0, formula)
                with open('IMAGINARY_FREQUENCY', 'w') as f:
                    f.write('IMAGINARY FREQUENCY exists {:s} {:s}\n'.format(dir0, formula))
        except RuntimeError:
            print('RuntimeError: Dynamical matrix has not yet built ', dir0, formula)

    else:
        struc_ase = pymatgen2ase(struc_pmg)
        #struc = sort(struc_ase)
        try:
            error = False
            try:
                struc, energy, time = single_optimize(struc_ase, 1)
                results = vasprun('vasprun.xml')
                error = results.error
            except ValueError:
                print('VASP is not converged')
                error = True
            except IndexError:
                print('VASP exits with error, problem in reading OUTCAR')
                error = True

            if not error:
                struc, energy, time = single_optimize(struc, 2)
                os.system('cp INCAR INCAR-Relax')
                os.system('cp OUTCAR OUTCAR-Relax')
                os.system('rm IBZKPT EIGENVAL XDATCAR DOSCAR vasprun.xml')
                #struc.write('POSCAR-unitcell', format='vasp', sort=True, vasp5=True, direct=True)
                struc.write('POSCAR-unitcell', format='vasp', vasp5=True, direct=True)
                finder = SpacegroupAnalyzer(ase2pymatgen(struc))
                spg = finder.get_space_group_symbol()
                spg_number = finder.get_space_group_number()

                cell = struc.get_cell()
                a = np.linalg.norm(cell[0,:])
                b = np.linalg.norm(cell[1,:])
                c = np.linalg.norm(cell[2,:])
                #if np.dot(cell[0,:], cell[1,:])>0:
                if 142 < spg_number < 195:
                    # increase the cell for hexagonal cell
                    a = a/1.5
                    b = b/1.5

                d1 = int(np.ceil(9.5/a))
                d2 = int(np.ceil(9.5/b))
                d3 = int(np.ceil(9.5/c))
                dim = str(d1) + ' ' + str(d2) + ' ' + str(d3)
                os.system('phonopy -d --dim="'+ dim +'" -c POSCAR-unitcell')
                with open('phonopy.conf', 'w') as f:
                    #So far, we consider only F, I, R, P
                    f.write('DIM = {:d} 0 0 0 {:d} 0 0 0 {:d}\n'.format(d1, d2, d3))
                    if spg[0] == 'F':
                        f.write('PRIMITIVE_AXIS = 0 1/2 1/2 1/2 0 1/2 1/2 1/2 0\n')
                    elif spg[0] == 'I':
                        f.write('PRIMITIVE_AXIS = -1/2 1/2 1/2 1/2 -1/2 1/2 1/2 1/2 -1/2\n')
                    elif spg[0] == 'R':
                        f.write('PRIMITIVE_AXIS = 2/3 -1/3 -1/3 1/3 1/3 -2/3 1/3 1/3 1/3\n')
                    elif spg[0] == 'P':
                        f.write('PRIMITIVE_AXIS = 1 0 0 0 1 0 0 0 1\n')
                    else:
                        f.write(spg)

                os.system('cp SPOSCAR POSCAR')
                struc = read('POSCAR', format='vasp')
                prepare_vasp(struc, results.values['gap'])
                os.system(cmd)
                os.system('phonopy --fc vasprun.xml')
                os.system('rm POSCAR CONTCAR POSCAR-0* CHG* DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT REPORT WAVECAR XDATCAR ase-sort.dat')
        except RuntimeError:
            print('No pseudopotential ', formula)
    os.chdir(cwd)
