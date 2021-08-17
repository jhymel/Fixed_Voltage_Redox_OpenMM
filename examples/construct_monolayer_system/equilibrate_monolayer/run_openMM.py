from __future__ import print_function
import sys
import os
#********** Fixed_Voltage_Redox_OpenMM
repo_path = '../../../' # Path to repo head
lib_path = os.path.join(repo_path,'lib')
sys.path.append(lib_path)
ffdir_path = os.path.join(repo_path,'ffdir')
#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from MM_classes_FV_Redox import *
#***************************
import numpy as np
import urllib.request
from FV_shared.lib.add_customnonbond_xml import add_CustomNonbondedForce_SAPTFF_parameters
# other stuff
from sys import stdout
from time import gmtime, strftime
from datetime import datetime

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

# *********************************************************************
#                Fixed-Voltage MD code for electrode simulations
#
#    This code has been generalized to work for flat electrodes (e.g. graphene), and electrodes
#    functionalized with conducting objects such as spherical buckyballs or cylindrical nanotubes
#
#    This includes a module for MC equilibration of the initial system to equilibrate the density.
#**********************************************************************

# a few run control settings...   WARNING:  write_charge = True will generate a lot of data !!!
simulation_time_ns = 10.0 ; freq_charge_update_fs = 200 ; freq_traj_output_ps = 10 ; freq_checkpoint_ps = 10 ; write_charges = False; freq_lambda_derivative_fs = 100*freq_charge_update_fs

flexible_surfactant = True
FV_bool = True # set whether or not the FV routine will be used

analytic_charge_scaling = False # Turn on/off analytic correction in Poisson solver

checkpoint = 'state.chk'
charge_name = 'charges.dat'
if write_charges:
    if os.path.exists(charge_name):
        chargeFile = open(charge_name, "a")
    else:
        chargeFile = open(charge_name, "w")

lambda_value = 0.0
lambda_step_size = 0.1
derivative_name = 'derivative.out'
if os.path.exists(derivative_name):
    derivative_output_file = open(derivative_name, "a")
else:
    derivative_output_file = open(derivative_name, "w")

rigid_body = ('FE','Cr','Hr')

# Issue with thiol head of the surfactant causing cathode/anode asymmetry at 0V, exlusion resolves this issue
surfactant_exceptions=('HS','S','CT1','HT1', 'HT2', 'CT2', 'HT3', 'HT4')
surfactant_exceptions_index=[]


#******************************************************************
#                Choose type of simulation
simulation_type = "MC_equil"  # either "Constant_V" or "MC_equil"
#**********************************************************************

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# set applied voltage in Volts.  This won't be used for "MC_equil" simulation...
Voltage = 0.0 # in Volts, units will be internally converted later...

# electrode names used to exclude intra-electrode non-bonded interactions ...
#cathode_name="cath"; anode_name = "anod"
# use chain indices instead of residue names to identify electrode...
#cathode_index=(3,2,1,0); anode_index=(7,4,5,6) # note chain indices start at 0 ...
#icathode_index=(); anode_index=() # note chain indices start at 0 ...
cathode_index=(0,4,3,2); anode_index=(1,7,6,5) 

# Initialize: Input list of pdb and xml files
home = os.getcwd()
starting_structure = sys.argv[1]
pdb_path = home
residue_xml_list = ['gold_c_residue.xml', 'gold_a_residue.xml', 'gold_residue.xml', 'water_ion_residue.xml', 'Fe_surfactant11_residue.xml'] # Standard residue xml files
ff_xml_list_neutral_flex = ['gold_cathode.xml', 'gold_anode.xml', 'gold.xml', 'gold-water-Fe-flex.xml', 'Fe_surfactant11_flex.xml']
ff_xml_list_neutral_fixed = ['gold_cathode.xml', 'gold_anode.xml', 'gold.xml', 'gold-water.xml', 'Fe_surfactant11_fixed.xml']
ff_xml_list_oxidized_flex = ['gold_cathode.xml', 'gold_anode.xml', 'gold.xml', 'gold-water_FE_oxidized.xml', 'Fe_surfactant11_flex.xml']
ff_xml_list_oxidized_fixed = ['gold_cathode.xml', 'gold_anode.xml', 'gold.xml', 'gold-water_FE_oxidized.xml', 'Fe_surfactant11_fixed.xml']

# Create two MMsys objects, one for each oxidation state
# Different forcefield/xml files used depending on oxidation state and flexiblity of the surfactant
if flexible_surfactant == True:
    MMsys_neutral=MM_FixedVoltage_Redox( pdb_list = [ os.path.join(pdb_path,starting_structure), ] , residue_xml_list = [ os.path.join(ffdir_path,f) for f in residue_xml_list ] , ff_xml_list = [ os.path.join(ffdir_path,f) for f in ff_xml_list_neutral_flex ] , analytic_charge_scaling = analytic_charge_scaling , rigid_body=rigid_body , FV_solver = FV_bool)
    MMsys_oxidized=MM_FixedVoltage_Redox( pdb_list = [ os.path.join(pdb_path,starting_structure), ] , residue_xml_list = [ os.path.join(ffdir_path,f) for f in residue_xml_list ] , ff_xml_list = [ os.path.join(ffdir_path,f) for f in ff_xml_list_oxidized_flex ] , analytic_charge_scaling = analytic_charge_scaling )
else:
    MMsys_neutral=MM_FixedVoltage_Redox( pdb_list = [ os.path.join(pdb_path,starting_structure), ] , residue_xml_list = [ os.path.join(ffdir_path,f) for f in residue_xml_list ] , ff_xml_list = [ os.path.join(ffdir_path,f) for f in ff_xml_list_neutral_fixed ] , analytic_charge_scaling = analytic_charge_scaling , FV_solver = FV_bool)
    MMsys_oxidized=MM_FixedVoltage_Redox( pdb_list = [ os.path.join(pdb_path,starting_structure), ] , residue_xml_list = [ os.path.join(ffdir_path,f) for f in residue_xml_list ] , ff_xml_list = [ os.path.join(ffdir_path,f) for f in ff_xml_list_oxidized_fixed ] , analytic_charge_scaling = analytic_charge_scaling )

constain_atom = "S"
k = 166188.0
with open(starting_structure, 'r') as f:
    sulf_z_total = 0
    sulf_count = 0
    for line in f:
        if ' %s '%constain_atom in line:
            sulf_z_total += float(line.split()[7])
            sulf_count += 1
average_z = round(sulf_z_total/sulf_count,3)/10
MMsys_neutral.set_position_constraint(constain_atom, k, average_z)

# if periodic residue, call this
MMsys_neutral.set_periodic_residue(True)

# Set platform and create simulation object (simmd) inside of MMsys class
MMsys_neutral.set_platform('CUDA')
MMsys_oxidized.set_platform('CPU')

for res in MMsys_neutral.simmd.topology.residues():
    for atom in res._atoms:
        if atom.name in surfactant_exceptions:
            surfactant_exceptions_index.append( atom.index )

# dictionary containng key:value pairs of resid:[atom_indices]
resid_dict = filter_redox_molecules(atom_types=('FE','Cr','Hr'), MMsys_redox0 = MMsys_neutral, MMsys_redox1 = MMsys_oxidized)

# redox_system : list of redox_molecule objects
# active_redox_molecule : redox_molecule object that actually has it's charges scaled by lambda in the simulation, selected via a resid
redox_system = []
active_redox_molecule_resid = 1
for counter, resid in enumerate(resid_dict):
    redox_system.append(redox_molecule(resid_dict[resid], nbondedForce0 = MMsys_neutral.nbondedForce, nbondedForce1 = MMsys_oxidized.nbondedForce))
    if resid == active_redox_molecule_resid:
        active_redox_molecule = redox_system[counter]

set_redox_system_oxidized = False
if set_redox_system_oxidized:
    for mol in redox_system:
        mol.scale_charges_by_lambda(context = MMsys_neutral.simmd.context, lambda_scale = 1.0)

#***********  Initialze OpenMM API's, this method creates simulation object

# initialize Virtual Electrodes, these are electrode `sheets' that solve electrostatics for constant Voltage ...
# can choose electrodes by residue name (this is default)
# can also choose electrodes by chain name (set chain=True)
# can input tuple exclude_element with elements to exclude from virtual electrode, such as dummy Hydrogen atoms ...i
MMsys_neutral.initialize_electrodes( Voltage, cathode_identifier = cathode_index , anode_identifier = anode_index , chain=True , exclude_element=("H",) )  # chain indices instead of residue names
MMsys_oxidized.initialize_electrodes( Voltage, cathode_identifier = cathode_index , anode_identifier = anode_index , chain=True , exclude_element=("H",) )  # chain indices instead of residue names
# initialize atoms indices of electrolyte, we need this for analytic charge correction.  Currently we electrode residue > 100 atoms, electrolyte residue < 100 atoms ... this should be fine?
MMsys_neutral.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...
MMsys_oxidized.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...

# IMPORTANT: generate exclusions for SAPT-FF.  If flag_SAPT_FF_exclusions = True , then will assume SAPT-FF force field and put in appropriate exclusions.
# set flag_SAPT_FF_exclusions = False if not using SAPT-FF force field


#MMsys_neutral.generate_exclusions( water_name = 'HOH' , flag_hybrid_water_model = True , flag_SAPT_FF_exclusions = False , surfactant_exclude_index = surfactant_exceptions_index )
# Need to call generate_exclusions on MMsys_oxidized since charges are scaled for whole surfactant
print ('going to try zeroing charges...')
MMsys_neutral.generate_exclusions( water_name = 'HOH' , flag_hybrid_water_model = True , flag_SAPT_FF_exclusions = False , surfactant_zero_charge_index = surfactant_exceptions_index )
MMsys_oxidized.generate_exclusions( water_name = 'HOH' , flag_hybrid_water_model = True , flag_SAPT_FF_exclusions = False , surfactant_zero_charge_index = surfactant_exceptions_index )
print ('charges zeroed succesfully!')

# load checkpoint for restarting if it is here ...
if os.path.exists(checkpoint):
    print( "restarting simulation from checkpoint ", checkpoint )
    MMsys_neutral.simmd.loadCheckpoint(checkpoint)

# write initial pdb with Drudes, and setug trajectory output
state = MMsys_neutral.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
positions = state.getPositions()
PDBFile.writeFile(MMsys_neutral.simmd.topology, positions, open('start_drudes.pdb', 'w'))

if simulation_type == "MC_equil":
    # Monte Carlo equilibration, initialize parameters ...  currently set to move Anode, keep Cathode fixed...
    celldim = MMsys_neutral.simmd.topology.getUnitCellDimensions()
    MMsys_neutral.MC = MC_parameters(  MMsys_neutral.temperature , celldim , electrode_move="Anode" , pressure = 1000.0*bar , barofreq = 100 , shiftscale = 0.2 )
    trajectory_file_name = 'equil_MC.dcd'
else :
    trajectory_file_name = 'FV_NVT.dcd'

append_trajectory = False
MMsys_neutral.set_trajectory_output( trajectory_file_name , freq_traj_output_ps * 1000 , append_trajectory , checkpoint , freq_checkpoint_ps * 1000 )

if lambda_value == 0.0:
    print ('start minimizing!')
    MMsys_neutral.simmd.minimizeEnergy(maxIterations=100)
    print ('done minimizing!')

#for lam in np.arange(0,1+lambda_step_size,lambda_step_size):
for lam in [lambda_value]:

    print ('setting lambda = %s' % lam)
    derivative_output_file.write('setting lambda = %s\n' % lam)
    active_redox_molecule.scale_charges_by_lambda(context = MMsys_neutral.simmd.context, lambda_scale = lam)

    cur_FE_charges = []
    for res in MMsys_neutral.simmd.topology.residues():
        if int(res.id) in resid_dict:
            for atom in res._atoms:
                if atom.name[0:2] == 'FE':
                    (q_i, sig, eps) = MMsys_neutral.nbondedForce.getParticleParameters(atom.index)
                    cur_FE_charges.append(q_i._value)
    print ("initial charges on FE's", cur_FE_charges)

    for i in range( int(simulation_time_ns * 1000 / freq_traj_output_ps ) ):
        state = MMsys_neutral.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=False,getPositions=True)
        print(str(state.getKineticEnergy()), flush=True)
        print(str(state.getPotentialEnergy()))
        for j in range(MMsys_neutral.system.getNumForces()):
            f = MMsys_neutral.system.getForce(j)
            print(type(f), str(MMsys_neutral.simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))
        print(" charge on Cathode/Anode " , MMsys_neutral.Cathode.get_total_charge() , MMsys_neutral.Anode.get_total_charge() )
        print( "analytic charges ", MMsys_neutral.Cathode.Q_analytic , MMsys_neutral.Anode.Q_analytic )
 
        #**********  Monte Carlo Simulation ********
        if simulation_type == "MC_equil":
            for j in range( int(freq_traj_output_ps * 1000 / MMsys_neutral.MC.barofreq) ):
                MMsys_neutral.MC_Barostat_step()

        #**********  Constant Voltage Simulation ****
        elif simulation_type == "Constant_V":
            for j in range( int(freq_traj_output_ps * 1000 / freq_charge_update_fs) ):
                
                # Fixed Voltage Electrostatics ..
                MMsys_neutral.Poisson_solver_fixed_voltage()
                MMsys_neutral.simmd.step( freq_charge_update_fs )
                
                if (j * freq_charge_update_fs) % freq_lambda_derivative_fs == 0:
                    dE_dlambda = dE_dLambda(active_redox_molecule = active_redox_molecule, MMsys = MMsys_neutral, lambda_scale = lam, dlambda=0.02)
                    derivative_output_file.write(str(dE_dlambda) + '\n')
                    derivative_output_file.flush() # flush buffer

            if write_charges :
                # write charges...
                MMsys_neutral.write_electrode_charges( chargeFile )
            

        else:
            print('simulation type not recognized ...')
            sys.exit()

print('done!')
sys.exit()





