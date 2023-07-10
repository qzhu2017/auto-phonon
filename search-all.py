from pymatgen import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from optparse import OptionParser
from pymatgen.io.vasp.inputs import Poscar
import pandas as pd
import pymatgen.analysis.structure_matcher as sm
from tabulate import tabulate

mpr = MPRester('fn8WQGgT9rvZAh6H')

def Read_POSCARS(filename):
    """
    Read POSCARS in a single file and convert it to Pymatgen structure object
    """
    # Read INPUT  file
    Struc = []
    with open(filename, 'rb') as f:
        input_content = f.readlines()

    POSCAR_content = []
    N_atom = 0
    for str1 in input_content:
        POSCAR_content.append(str(str1, 'utf-8'))
        if len(POSCAR_content) == 7:
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
    return Struc


parser = OptionParser()
parser.add_option("-f",  "--file", dest="posfile",
                  help="by filename in POSCAR format")
(options, args) = parser.parse_args()

content = []
if options.posfile:
    strucs = Read_POSCARS(options.posfile)
    finder = SpacegroupAnalyzer(strucs[0])
    spg = finder.get_space_group_symbol()
    P_struc = finder.get_primitive_standard_structure()
    nelements = len(P_struc.composition)
    nsites = len(P_struc.sites)

else:
    strucs = []
    spg = 'P6_3mc'
    nelements = 2
    nsites = 4

formula = []
spgs = []
mid = []
e_hull = []
criteria ={
       "e_above_hull": {"$lt":0.02},
       "spacegroup.symbol": spg,
       "nelements": nelements,
       "nsites":nsites,
      }
print(spg, nelements, nsites)
properties = ["structure", "pretty_formula", "material_id", "spacegroup.symbol","e_above_hull"]
entries = mpr.query(criteria=criteria, properties=properties)
print('Total number of returns in ', spg, ':  ', len(entries))
for entry in entries:
    struc = entry['structure']
    include = True
    if len(strucs) > 0:
        if not sm.StructureMatcher().fit_anonymous(struc, strucs[0]):
            include = False
    if include:
        formula.append(entry['pretty_formula'])
        mid.append(entry['material_id'])
        spgs.append(entry['spacegroup.symbol'])
        e_hull.append(entry['e_above_hull'])
        new = True
        for i, s in enumerate(strucs):
            #print(i, s.composition, sm.StructureMatcher().fit(struc, s))
            if sm.StructureMatcher().fit(struc, s):
                #print(entry['pretty_formula'], ' exists')
                new = False
                break
        #print(entry['pretty_formula'], new)
        if new:
            tmp = struc.to(fmt='poscar')
            tmp = entry['material_id'] + '|' + tmp
            print(tmp)
            content.append(tmp)

with open(options.posfile.split('/')[-1]+'_todo', 'w') as f:
    f.writelines(content)
    
col_name = {'Formula': formula,
            'Space group': spgs,
            'material_id': mid,
            'e_above_hull': e_hull,
            }

df = pd.DataFrame(col_name)
print(tabulate(df, headers='keys', tablefmt='psql'))

