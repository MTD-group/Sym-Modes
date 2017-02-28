import pymatgen as mg
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter
import os

if os.path.isfile('files.txt'): # For a list of folder names with POSCARs inside them, will write their symmetrized versions
    with open('files.txt') as f:
        cwd = os.getcwd()
        for line in f.readlines():
            path = cwd + '/' +  line.split()[0]
            os.chdir(path)
            p = Poscar.from_file('POSCAR')
            z = CifWriter(p.structure, symprec=0.001)
            z.write_file(cwd + '/pymat_cifs/' +  '%s_pymat.cif' % line.split()[0])
            
elif os.path.isfile('cif_files.txt'): # For a list of cif filenames, will write their symmetrized versions 
    with open('cif_files.txt') as f:
        for line in f.readlines():
            path = line.split()[0]
            os.mkdir("./pymat_cifs")
            struct = mg.Structure.from_file(path)
            z = CifWriter(struct, symprec=0.001)
            z.write_file('./pymat_cifs/' +  '%s_pymat.cif' % path.split(".")[0])
            
else: # If you just want the symmetrized cif for one POSCAR
    p = Poscar.from_file('POSCAR')
    z = CifWriter(p.structure, symprec=0.001)
    with open('POSCAR') as f:
        your_name = f.readline().split()[0]
    z.write_file(your_name + '.cif')

