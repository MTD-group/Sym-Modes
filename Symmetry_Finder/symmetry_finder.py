from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter
import os

if os.path.isfile('files.txt'):
    with open('files.txt') as f:
        cwd = os.getcwd()
        for line in f.readlines():
            path = cwd + '/' +  line.split()[0]
            os.chdir(path)
            p = Poscar.from_file('POSCAR')
            z = CifWriter(p.structure, symprec=0.001)
            z.write_file(cwd + '/cifs/' +  '%s_pymat.cif' % line.split()[0])
        
else:
    p = Poscar.from_file('POSCAR')
    z = CifWriter(p.structure, symprec=0.001)
    your_name = raw_input('Please provide your cif name without .cif: ')
    z.write_file(your_name + '.cif')

