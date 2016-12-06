
# coding: utf-8

# In[1]:

from pymatgen.core import Structure
import Tkinter as tk
import tkFileDialog
import tkSimpleDialog

while True:
    try:
        guimode = raw_input("GUI mode? (yes/no): ")
    except ValueError:
        print("Sorry, I didn't understand that.")
        continue
    if guimode.lower() in ('yes', 'y', 'no', 'n'):
        break
    else:
        print('Please enter yes or no')

# Read in cif files for structures without and with modes
if guimode.lower() == 'yes' or guimode.lower() == 'y':
    print("Please provide parent structure cif with no modes: ")
    file_path = tkFileDialog.askopenfilename(title="parent structure cif")
    no_modes = Structure.from_file(file_path)
    numberomodes = tkSimpleDialog.askinteger('Number of Modes', 'Please enter number of modes to interpolate (1-3): '                  ,minvalue=1, maxvalue=3)
    file_path_mode1 = tkFileDialog.askopenfilename(title="Cif with 100% of your desired mode")
    mode = Structure.from_file(file_path_mode1)

    if numberomodes >= 2:
        print("Please provide cif with 100% of your second desired mode: ")
        file_path_mode2 = tkFileDialog.askopenfilename(title="Cif with 100% of your second desired mode")
        mode2 = Structure.from_file(file_path_mode2)
    if numberomodes >= 3:
        print("Please provide cif with 100% of your third desired mode: ")
        file_path_mode3 = tkFileDialog.askopenfilename(title="Cif with 100% of your third desired mode")
        mode3 = Structure.from_file(file_path_mode3)
else:
    file_path = raw_input("Please provide parent structure cif with no modes: ")
    no_modes = Structure.from_file(file_path)
    numberomodes = int(raw_input("Please enter number of modes to interpolate (1-3): "))
    
    file_path_mode1 = raw_input("Please cif with 100% of your desired mode: ")
    mode = Structure.from_file(file_path_mode1)
    
    if numberomodes >= 2:
        file_path_mode2 = raw_input("Please cif with 100% of your second desired mode: ")
    if numberomodes >= 3:
        file_path_mode3 = raw_input("Please cif with 100% of your third desired mode: ")


# In[24]:

from pymatgen.analysis.structure_matcher import *
matchy = StructureMatcher(primitive_cell=False, attempt_supercell=True, scale=False)

trans_nomode = matchy.get_s2_like_s1(mode, no_modes) # Changes setting but keeps relative positions constant
trans_nomode.modify_lattice(mode.lattice)        # Adjusts lattice parameters to match that of daughter structures

if numberomodes == 2:
    if len(mode) > len(mode2):
        trans_nomode = matchy.get_s2_like_s1(mode, no_modes)
        trans_mode2 = matchy.get_s2_like_s1(trans_nomode, mode2)
        trans_nomode.modify_lattice(mode.lattice)        
        
    else:
        trans_nomode = matchy.get_s2_like_s1(mode2, no_modes)
        trans_mode1 = matchy.get_s2_like_s1(trans_nomode, mode)
        trans_nomode.modify_lattice(mode2.lattice)
        
if numberomodes == 3:
    if len(mode) > len(mode2) and len(mode) > len(mode3):
        trans_nomode = matchy.get_s2_like_s1(mode, no_modes)
        trans_mode2 = matchy.get_s2_like_s1(trans_nomode, mode2)
        trans_mode3 = matchy.get_s2_like_s1(trans_nomode, mode3)
        trans_nomode.modify_lattice(mode.lattice)
        
    elif len(mode2) > len(mode3):
        trans_nomode = matchy.get_s2_like_s1(mode2, no_modes)
        trans_mode1 = matchy.get_s2_like_s1(trans_nomode, mode)
        trans_mode3 = matchy.get_s2_like_s1(trans_nomode, mode3)
        trans_nomode.modify_lattice(mode2.lattice)
        
    else:
        trans_nomode = matchy.get_s2_like_s1(mode3, no_modes)
        trans_mode1 = matchy.get_s2_like_s1(trans_nomode, mode)
        trans_mode2 = matchy.get_s2_like_s1(trans_nomode, mode2)
        trans_nomode.modify_lattice(mode3.lattice)
        


# In[25]:

nm_coords = trans_nomode.frac_coords
subby = []

if numberomodes == 1:
    mode_coords = mode.frac_coords
    for atom_index in range(0, len(nm_coords)):
        subby.append([])
        for coord_index in range(0, 3):
            subby[atom_index].append(mode_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])

if numberomodes == 2:        
    try: # Use the scaled up version of the structure with less sites in POSCAR
        trans_mode2
    except NameError:
        mode1_coords = trans_mode1.frac_coords
        mode2_coords = mode2.frac_coords
    else:
        mode1_coords = mode.frac_coords
        mode2_coords = trans_mode2.frac_coords        

    for atom_index in range(0, len(nm_coords)):
        subby.append([])
        for coord_index in range(0, 3):
            subby[atom_index].append(mode1_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
    
    subby2 = []
    for atom_index in range(0, len(nm_coords)):
        subby2.append([])
        for coord_index in range(0, 3):
            subby2[atom_index].append(mode2_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
        
if numberomodes == 3:  
    try: # Use the scaled up version of the structure with less sites in POSCAR
        trans_mode3
    except NameError:
        mode1_coords = trans_mode1.frac_coords
        mode2_coords = trans_mode2.frac_coords
        mode3_coords = mode3.frac_coords
    else:
        try:
            trans_mode2
        except NameError:
            mode1_coords = trans_mode1.frac_coords
            mode2_coords = mode2.frac_coords
            mode3_coords = trans_mode3.frac_coords
        else:
            mode1_coords = mode.frac_coords
            mode2_coords = trans_mode2.frac_coords
            mode3_coords = trans_mode3.frac_coords
    
    
    for atom_index in range(0, len(nm_coords)):
        subby.append([])
    for coord_index in range(0, 3):
        subby[atom_index].append(mode1_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
    
    subby2 = []
    for atom_index in range(0, len(nm_coords)):
        subby2.append([])
        for coord_index in range(0, 3):
            subby2[atom_index].append(mode2_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])
    
    mode3_coords = mode3.frac_coords
    subby3 = []
    for atom_index in range(0, len(nm_coords)):
        subby3.append([])
        for coord_index in range(0, 3):
            subby3[atom_index].append(mode3_coords[atom_index][coord_index]-nm_coords[atom_index][coord_index])


# In[26]:

from pymatgen.io.vasp import Poscar
import copy

intervals = 10 
increment = 1.0 / intervals

if guimode.lower() == 'yes' or guimode.lower() == 'y':
    mode_filename=tkSimpleDialog.askstring('Mode Name','''Enter your mode's name''')
    if numberomodes >= 2:
        mode2_filename = tkSimpleDialog.askstring('2nd Mode Name','''Enter your second mode's name''')
    if numberomodes == 3:
        mode3_filename = tkSimpleDialog.askstring('3rd Mode Name','''Enter your third mode's name''')
    
    print('Enter output directory')
    output_path = tkFileDialog.askdirectory(title='Enter output directory')
    
    
else:     
    mode_filename = raw_input("Please enter first mode label: ")

    if numberomodes >= 2:
        mode2_filename = raw_input("Please enter second mode label: ")
    if numberomodes == 3:
        mode3_filename = raw_input("Please enter third mode label: ")

    output_path = raw_input("Please enter output directory: ")


if numberomodes == 1:
    for interval in range(1, intervals+1):
        ccopy = copy.copy(trans_nomode)
        species = ccopy.species
        subs = ccopy.frac_coords

        for atom in range(0, len(subby)):
            for coord in range(0, 3):
                subs[atom][coord] = subby[atom][coord] * interval * increment + nm_coords[atom][coord]

        # Add displacement to zero modes structure
        index = range(0, len(species))
        for i in index:
            ccopy.replace(i, species[i], coords=subs[i])

        # Export to POSCAR
        struct = Poscar(ccopy)
        value = increment*100*interval
        struct.write_file(output_path+"/%s_%s%%.vasp" % (mode_filename, value)) #e.g. 'Y2+_100%'
    
if numberomodes == 2:
    for interval2 in range(0, intervals + 1):
        for interval in range(0, intervals+1):
            ccopy = copy.copy(trans_nomode)
            species = ccopy.species
            subs = ccopy.frac_coords

            for atom in range(0, len(subby)):
                for coord in range(0, 3):
                    subs[atom][coord] = subby[atom][coord] * interval * increment +                                         subby2[atom][coord] * interval2 * increment + nm_coords[atom][coord]

            # Add displacement to zero modes structure
            index = range(0, len(species))
            for i in index:
                ccopy.replace(i, species[i], coords=subs[i])

            # Export to POSCAR
            struct = Poscar(ccopy)
            value = increment*100*interval
            value2 = increment*100*interval2
            struct.write_file(output_path+"/%s_%s%%_%s_%s%%.vasp" % (mode2_filename, value2, mode_filename, value))
    
if numberomodes == 3:
    for interval3 in range(0, intervals +1):
        for interval2 in range(0, intervals + 1):
            for interval in range(0, intervals+1):
                ccopy = copy.copy(trans_nomode)
                species = ccopy.species
                subs = ccopy.frac_coords

                for atom in range(0, len(subby)):
                    for coord in range(0, 3):
                        subs[atom][coord] = subby[atom][coord] * interval * increment +                                             subby2[atom][coord] * interval2 * increment +                                             subby3[atom][coord] * interval3 * increment + nm_coords[atom][coord]

                # Add displacement to zero modes structure
                index = range(0, len(species))
                for i in index:
                    ccopy.replace(i, species[i], coords=subs[i])

                # Export to POSCAR
                struct = Poscar(ccopy)
                value = increment*100*interval
                value2 = increment*100*interval2
                value3 = increment*100*interval3
                struct.write_file(output_path+"/%s_%s%%_%s_%s%%_%s_%s%%.vasp" % (mode3_filename, value3,                                                             mode2_filename, value2, mode_filename, value))

