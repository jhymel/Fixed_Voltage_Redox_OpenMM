# John Hymel : 5/30/21

# This script was created to add electrolyte molecules to the pdb produced by build_monolayer.py
# Need to have created a pdb containing the monolayer+electrode and a pdb containin the electrolyte

# Example run: python monolayer_electrolyte_combiner.py Fc11-C11-4_0.0487.pdb electrolyte.pdb 

import os
import sys

monolayer_pdb = sys.argv[1]
electrolyte_pdb = sys.argv[2]
out_file = 'monolayer.pdb'

with open(monolayer_pdb, 'r') as mono_in:
	with open(electrolyte_pdb, 'r') as electro_in:
		with open(out_file, 'w+') as out:
			for mono_line in mono_in:
				if mono_line != 'END':
					out.write(mono_line)
			for electro_line in electro_in:
				if 'ATOM' in electro_line:
					out.write(electro_line)
			out.write('END')
				
