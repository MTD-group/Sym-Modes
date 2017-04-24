
# coding: utf-8   
    

# In[5]:

# Cite: "Descriptors of oxygen evolution activity for oxides: a statistical evaluation"
# WT Hong, RE Welsch, Y Shao-Horn
# This script queries the Materials Project database for ABO3 perovskites, where B is a transition metal. Exports data on the bond lengths, energies, and Madelung potentials using pymatgen
from pymatgen import MPRester, Composition, Structure
from pymatgen.analysis.bond_valence import calculate_bv_sum
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.util import coord_utils
from math import pi, cos, sin, acos, ceil
import numpy as np
import csv
import pandas as pd

def MObonds(structure,Bsite):
# This function takes a pymatgen structure and perovskite Bsite and returns a list of the bond lengths associated with the Bsite/oxygen bond lengths
    bond_lengths = []
    # determine Bsite and oxygen indexes
    for site in structure.sites:
        if Bsite in str(site):
            bonds = structure.get_neighbors(site, r = 2.7)
            for bond in range(len(bonds)):
                bond_lengths.append(bonds[bond][1])
    return bond_lengths	

def MMbonds(structure,Bsite):
# This function takes a pymatgen structure and perovskite Bsite and returns a list of the lengths associated with the Bsite/Bsite distance
    distances = []
    # determine Bsite indexes
    for site in structure.sites:
        if Bsite in str(site):
            neighbors = structure.get_neighbors(site, r = 5)
            for neighbor in range(len(neighbors)):
                if Bsite in str(neighbors[neighbor][0]):
                    distances.append(neighbors[neighbor][1])
    return distances

def frac2cart(structure, coords):
# Converts fractional coordinates to cartesian coordinates
    a = structure.lattice.a
    b = structure.lattice.b
    c = structure.lattice.c
    alpha = structure.lattice.alpha*pi/180
    beta = structure.lattice.beta*pi/180
    gamma = structure.lattice.gamma*pi/180
    # transformation matrix for fractional coords to cartesian coords
    v = np.sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
    transform = np.array([[a, b*cos(gamma), c*cos(beta)], [0, b*sin(gamma), c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)], [0, 0, c*v/sin(gamma)]])
    abc = np.matrix([[a], [b], [c]])
    return np.dot(transform, coords)
    
def get_angle(A, B, C):
# This function calculates the angle between three points in cartesian space
    AB = B-A
    BC = C-B
    cos_angle = np.dot(AB,BC)/(np.linalg.norm(AB)*np.linalg.norm(BC))
    cos_angle = min(1,max(cos_angle,-1))	# clean to prevent computational errors with floating points
    angle = acos(cos_angle)
    return angle*180/pi
    
def MOMangles(structure,Bsite):
# This function takes a pymatgen structure and perovskite Bsite and returns a list of the MOM angles
    angles = []
    site_angles = []
    # determine Bsite indexes
    for site in structure.sites:
        if Bsite in str(site):
            neighbors = structure.get_neighbors(site, r = 5, include_index=True)
            Oneighbors = structure.get_neighbors(site, r = 4, include_index=True)
            for neighbor in neighbors:
                if Bsite in str(neighbor[0]):
                    for Osite in Oneighbors:
                        if 'O' in str(Osite[0]):
                            # calculate angle between sites
                            angle = coord_utils.get_angle(site.coords-Osite[0].coords, 
                                             neighbor[0].coords-Osite[0].coords, units="degrees")
                            site_angles.append(angle)

                    angles.append(max(site_angles))
                    site_angles = []
    return angles

def oxstates(structure, Bsite):
# This function estimates the oxidation state of the Bsite using bond valence methods
    sites = structure.sites
    Bsite_indexes = []
    Bsite_ox = []
    oxstate = 0
    for site in sites:
        if Bsite in str(site):
            neighbors = structure.get_neighbors(site, r = 2.5)
            oxstate = calculate_bv_sum(site, neighbors)
            Bsite_ox.append(oxstate)
            Bsite_specie = site.specie
        if Bsite not in str(site) and 'O' not in str(site):
            Asite_specie = site.specie

    Bsite_ox = ceil(sum(Bsite_ox)/len(Bsite_ox))
    Asite_ox = 6-Bsite_ox
    # catch unusual oxidation states resulting from calculating on the Bsite oxidation state
    if Asite_ox not in Asite_specie.common_oxidation_states or Bsite_ox not in Bsite_specie.common_oxidation_states:
        if 3 in Asite_specie.common_oxidation_states:
            Asite_ox = 3
        else:
            Asite_ox = max(Asite_specie.common_oxidation_states)
        Bsite_ox = 6-Asite_ox
    final_ox = {str(Asite_specie): Asite_ox, Bsite: Bsite_ox, 'O': -2}
    return final_ox

def EMad(structure, Bsite):
# Determines the site Madelung energies for the Bsite and Osite
    B_indexes = []
    O_indexes = []
    for index in range(len(structure.sites)):
        if Bsite in str(structure.sites[index]):
            B_indexes.append(index)
        elif 'O' in str(structure.sites[index]):
            O_indexes.append(index)

    # calculate Madelung energy
    mad_energy = EwaldSummation(structure)
    site_energies = np.array([])
    site_energies = sum(mad_energy.total_energy_matrix)

    # max V_mad_M and min V_mad_O gives the min CT energy
    V_mad_M = max(site_energies[B_indexes])
    V_mad_O = min(-site_energies[O_indexes])

    return [V_mad_M, V_mad_O]

def CN(site, structure, r):
# This function returns the coordination of a site given a radius defined bond length
    site_indexes = []
    for index in range(len(structure.sites)):
        if site in str(structure.sites[index]):
            site_indexes.append(index)
    coord_num = []
    for index in site_indexes:
        neighbors = structure.get_neighbors(structure.sites[index], r)
        coord_num.append(len(neighbors))
    return sum(coord_num)/len(coord_num)


# In[7]:

# I iterated through structures manually to get feature vectors
struct = Structure.from_file("Structures/Insulators/Sm2O3_LT_40475.cif")
elem_site = 'Sm'

arr = [float('{:0.4f}'.format(np.average(MOMangles(struct, elem_site)))), 
      float('{:0.4f}'.format(np.min(MOMangles(struct, elem_site)))),
      float('{:0.4f}'.format(np.max(MOMangles(struct, elem_site)))),
      float('{:0.4f}'.format(np.average(MObonds(struct, elem_site)))),
      float('{:0.4f}'.format(np.average(MMbonds(struct, elem_site)))),
      float('{:0.4f}'.format(np.min(MObonds(struct, elem_site)))),
      float('{:0.4f}'.format(np.min(MMbonds(struct, elem_site)))),
      float('{:0.4f}'.format(np.max(MObonds(struct, elem_site)))),
      float('{:0.4f}'.format(np.max(MMbonds(struct, elem_site)))),
      float('{:0.4f}'.format(EMad(struct, elem_site)[0])),
      float('{:0.4f}'.format(EMad(struct, elem_site)[1]))]
df = pd.Series(arr)
df.to_csv("./garbage.csv") # Write to temporary csv for copying


