
# coding: utf-8

# In[1]:

from pymatgen import *
import Tkinter as tk
import tkFileDialog
# help(pymatgen)


# In[14]:

# Read in cif file for structure with strain but zero displacive mode amplitude
root = tk.Tk()
root.withdraw()
print("Please provide parent structure cif with no modes: ")
file_path = tkFileDialog.askopenfilename()
no_modes = Structure.from_file(file_path)

numberomodes = int(raw_input("Please enter number of modes to interpolate (1-3): "))

# Read cif file of distortion of interest with strain at 100 % of its value
print("Please provide cif with 100% of your desired mode: ")
file_path2 = tkFileDialog.askopenfilename()
mode = Structure.from_file(file_path2)

if numberomodes >= 2:
    print("Please provide cif with 100% of your second desired mode: ")
    file_path_mode2 = tkFileDialog.askopenfilename()
    mode2 = Structure.from_file(file_path_mode2)
if numberomodes == 3:
    print("Please provide cif with 100% of your third desired mode: ")
    file_path_mode3 = tkFileDialog.askopenfilename()
    mode3 = Structure.from_file(file_path_mode3)


# all_modes = Structure.from_file("LaVO3(monoclinic)_relaxed_findsym.cif")

# print(no_modes)
# print(mode)


# In[16]:

from pymatgen.analysis.structure_matcher import *
matchy = StructureMatcher(primitive_cell=False)
matchy.fit(no_modes, mode)

trans_mode = matchy.get_s2_like_s1(no_modes, mode) # Arranges atoms in 2nd str to be similar to 1st str
trans_mode.modify_lattice(no_modes.lattice)

if numberomodes >= 2:
    matchy = StructureMatcher(primitive_cell=False, scale=False)
    matchy.fit(no_modes, mode2)
    trans_mode2 = matchy.get_s2_like_s1(no_modes, mode2) # Arranges atoms in 2nd str to be similar to 1st str
    trans_mode2.modify_lattice(no_modes.lattice)

if numberomodes == 3:
    matchy = StructureMatcher(primitive_cell=False, scale=False)
    matchy.fit(no_modes, mode3)
    trans_mode3 = matchy.get_s2_like_s1(no_modes, mode3) # Arranges atoms in 2nd str to be similar to 1st str
    trans_mode3.modify_lattice(no_modes.lattice)

# matchy.fit(no_modes, all_modes)
# trans_all = matchy.get_s2_like_s1(no_modes, all_modes)
# print(all_modes)
#print(mode)

#print(trans_mode)
#print(no_modes)
# help(StructureMatcher)


# In[17]:

import operator, pprint
nm_coords = no_modes.frac_coords
mode_coords = trans_mode.frac_coords
subby = []
for atom_index in range(0, len(nm_coords)):
    subby.append([])
    for coord_index in range(0, 3):
        subby[atom_index].append(mode_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])

if numberomodes >= 2:        
    mode2_coords = trans_mode2.frac_coords
    subby2 = []
    for atom_index in range(0, len(nm_coords)):
        subby2.append([])
        for coord_index in range(0, 3):
            subby2[atom_index].append(mode2_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
        
if numberomodes == 3:        
    mode3_coords = trans_mode3.frac_coords
    subby3 = []
    for atom_index in range(0, len(nm_coords)):
        subby3.append([])
        for coord_index in range(0, 3):
            subby3[atom_index].append(mode3_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
        
        
# no_modes.
pp = pprint.PrettyPrinter()
# print(nm_coords)
#print(mode_coords)
# pp.pprint(subby)


# In[19]:

from pymatgen.io.vasp import Poscar

intervals = 10 
increment = 1.0 / intervals
mode_filename = raw_input("Please enter first mode label: ")

if numberomodes >= 2:
    mode2_filename = raw_input("Please enter second mode label: ")
if numberomodes == 3:
    mode3_filename = raw_input("Please enter third mode label: ")



if numberomodes == 1:
    for interval in range(1, intervals+1):
        copy = Structure.from_file(file_path)
        species = copy.species
        subs = copy.frac_coords

        for atom in range(0, len(subby)):
            for coord in range(0, 3):
                subs[atom][coord] = subby[atom][coord] * interval * increment + nm_coords[atom][coord]

        # Add displacement to zero modes structure
        index = range(0, len(species))
        for i in index:
            copy.replace(i, species[i], coords=subs[i])

        # Export to POSCAR
        struct = Poscar(copy)
        value = increment*100*interval
        struct.write_file("%s_%s%%.vasp" % (mode_filename, value)) #e.g. 'Y2+_100%'
    
if numberomodes == 2:
    for interval2 in range(0, intervals + 1):
        for interval in range(0, intervals+1):
            copy = Structure.from_file(file_path)
            species = copy.species
            subs = copy.frac_coords

            for atom in range(0, len(subby)):
                for coord in range(0, 3):
                    subs[atom][coord] = subby[atom][coord] * interval * increment +                                         subby2[atom][coord] * interval2 * increment + nm_coords[atom][coord]

            # Add displacement to zero modes structure
            index = range(0, len(species))
            for i in index:
                copy.replace(i, species[i], coords=subs[i])

            # Export to POSCAR
            struct = Poscar(copy)
            value = increment*100*interval
            value2 = increment*100*interval2
            struct.write_file("%s_%s%%_%s_%s%%.vasp" % (mode2_filename, value2, mode_filename, value))
    
if numberomodes == 3:
    for interval3 in range(0, intervals +1):
        for interval2 in range(0, intervals + 1):
            for interval in range(0, intervals+1):
                copy = Structure.from_file(file_path)
                species = copy.species
                subs = copy.frac_coords

                for atom in range(0, len(subby)):
                    for coord in range(0, 3):
                        subs[atom][coord] = subby[atom][coord] * interval * increment +                                             subby2[atom][coord] * interval2 * increment +                                             subby3[atom][coord] * interval3 * increment + nm_coords[atom][coord]

                # Add displacement to zero modes structure
                index = range(0, len(species))
                for i in index:
                    copy.replace(i, species[i], coords=subs[i])

                # Export to POSCAR
                struct = Poscar(copy)
                value = increment*100*interval
                value2 = increment*100*interval2
                value3 = increment*100*interval3
                struct.write_file("%s_%s%%_%s_%s%%_%s_%s%%.vasp" % (mode3_filename, value3,                                                             mode2_filename, value2, mode_filename, value))    
    
    
    
    
    
# print(range(0, len(species)))
# copy.replace
# for i in range(0, increments):
    
# print(species[0])

# for i in range(0, len(subby)):
#     copy.__setitem__(s, subby[i])

#print(trans_mode)
#print(copy)
# pp.pprint(subs)

