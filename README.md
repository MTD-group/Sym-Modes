# Sym-Modes

This repo is for Python scripts created by the Materials Theory and Design Group to help with symmetry mode analysis.
If you have questions, Nicholas can be contacted at monsieurwagner@gmail.com

Structure_Aligner.py is designed to generate interpolated VASP structures between a high symmetry parent phase and a lower symmetry child phase. A group-subgroup relationship should exist between the two space groups, and both structures should be aligned so that c is the long axis. 

TolFactCalc is a Juptyer notebook for calculating the tolerance factors of ABX3 perovskites. It can be extended in a straightforward manner to other structure types such as RP-phases. It is intended for use with the Bond_valences2016 csv file. 
