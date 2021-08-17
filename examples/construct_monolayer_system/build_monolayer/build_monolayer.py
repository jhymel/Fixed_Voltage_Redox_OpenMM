import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as md

# Starting from a surfactant density of 4.6E-10 mol/cm^2
# found that corresponds with 0.027099 molecules per angstrom^2
density = 0.046 # particles per angstrom^2
box_length = 31.735 # sheet size, length of edge in angstrom

'''
# code for determining placement based on density and number of surfactants
density = 0.027099 # particles per angstrom^2
points = 10 # number of particles/surfacant molecules, needs to be a perfect square
total_area = points/density
box_length = (total_area)**(1/2)
point_area = total_area/points
side_length = (point_area)**(1/2)
'''

# Opens the pdb file containing the electrode and finds the index of the last atom
pdb_template = 'electrode.pdb'
with open(pdb_template, 'r') as f:
	pdb_data = f.readlines()
last_atom = pdb_data[-2].split()[1]
start_atom = int(last_atom) + 1

#input_file_Fc = '%s_only.pdb'%surfactant
input_file_Fc = 'Fc11.pdb'
input_file_CH = 'C11.pdb'

# Pull xyz coords from Fc_only file and center Fc on first lattice point
with open(input_file_CH, 'r') as f:
        data = f.readlines()
with open(input_file_Fc, 'r') as g:
	data_Fc = g.readlines()
x_max = 0.0; x_min = 1000.0; y_max = 0.0; y_min = 1000.0
for line in data:
        if float(line.split()[5]) > x_max:
                x_max = float(line.split()[5])
        if float(line.split()[5]) < x_min:
                x_min = float(line.split()[5])
        if float(line.split()[6]) > y_max:
                y_max = float(line.split()[6])
        if float(line.split()[6]) < y_min:
                y_min = float(line.split()[6])


buff = 0.1
x_min -= buff
y_min -= buff
x_width = x_max-x_min; y_width = y_max-y_min

shifted_Fc = []

total_area = box_length*box_length
points = density*total_area
point_area = total_area/points
side_length = np.sqrt(0.8)*(point_area)**(1/2)
#side_length = (point_area)**(1/2)

#x = side_length/2.0
#y = side_length/2.0

x = 0.0
y = 0.0

lattice = []

# Loop that creates the lattice of points where each sulfur should lie
count = 0
while x < box_length - max(x_width,y_width):
	while y < box_length - max(x_width,y_width):
		lattice.append([x,y])
		count += 1
		y += side_length
	x += side_length
	#y = side_length/2.0
	y = 0.0

print ('num surf',len(lattice))
print ('area',box_length*box_length)
print ('target density',density)
print ('actual density',len(lattice)/(box_length*box_length))
actual_density = len(lattice)/(box_length*box_length)
surfactant_Fc = input_file_Fc.split('.')[0]
surfactant_CH = input_file_CH.split('.')[0]
percent_Fc = 0.25
ratio_Fc = int(1/percent_Fc)
out_file = f"{surfactant_Fc}-{surfactant_CH}-{ratio_Fc}_{round(actual_density,4)}.pdb"

# Loop that writes the final pdb file
# Puts in the electrode. Loops over each lattice point, translating a surfactant molecule into the appropriate place
# on the electrode. Includes an if statement (if surf_count %(ratio_Fc) == 0:) for changing the Fc/C percentage
# of the monolayer
surf_count = 0
count = start_atom
with open(out_file, 'w+') as g:
	for line in pdb_data[0:len(pdb_data)-1]:
		g.write(line)
		line.split()
	for index, shift in enumerate(lattice):
		#print (surf_count,surf_count%2, surf_count%14, points)
		if surf_count %(ratio_Fc) == 0: # and surf_count%14 < 6:
			for linenumber, line in enumerate(data_Fc):
				g.write("%4s  %5s %4s %4s%s%4s    %8.3f%8.3f%8.3f  1.00  0.00      %4s%2s \n" %( line[0:4], count, line[12:16], 'FE_S', line[21:22], index + 1 , float(line[30:38]) - x_min + float(shift[0]), float(line[38:46]) - y_min + float(shift[1]), float(line[46:54]) + 0.7, line[72:76], line[76:78] ))
				count += 1
		else:
			for linenumber, line in enumerate(data):
				g.write("%4s  %5s %4s %4s%s%4s    %8.3f%8.3f%8.3f  1.00  0.00      %4s%2s" %( line[0:4], count, line[12:16], 'C11 ', line[21:22], index + 1 , float(line[30:38]) - x_min + float(shift[0]), float(line[38:46]) - y_min + float(shift[1]), float(line[46:54]) + 0.7, line[72:76], line[76:78] ))
				count += 1
		surf_count += 1
	g.write(pdb_data[-1])


