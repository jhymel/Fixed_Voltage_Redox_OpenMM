from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******************************************************

# Thermodynamic integration module

# Defines the redox_molecule class
#  - redox_molecule has charges scaled linearly between two oxidation states
# Computes change in energy w.r.t. change in scaling factor lambda

# redox_state variable: If == 0, system initializes as all redox0.
#                       If == 1, system initializes as all redox1.
#******************************************************

def dE_dLambda(active_redox_molecule, MMsys, lambda_scale, dlambda=0.02):
     
    # Optimize electrode charges before energy calculation
    MMsys.Poisson_solver_fixed_voltage()
    state = MMsys.simmd.context.getState(getEnergy=True,getForces=False,getVelocities=False,getPositions=False)
    energy_ref = state.getPotentialEnergy()
    lambda_1 = lambda_scale + dlambda
    active_redox_molecule.scale_charges_by_lambda(context = MMsys.simmd.context, lambda_scale = lambda_1)

    # Optimize electrode charges before energy calculation
    MMsys.Poisson_solver_fixed_voltage()
    state = MMsys.simmd.context.getState(getEnergy=True,getForces=False,getVelocities=False,getPositions=False)
    energy_1 = state.getPotentialEnergy()

    # Set charges back to their original value
    active_redox_molecule.scale_charges_by_lambda(context = MMsys.simmd.context, lambda_scale = lambda_scale)
    MMsys.Poisson_solver_fixed_voltage()

    # Output numerical derivative
    return (energy_1 - energy_ref) / dlambda

class redox_molecule(object):
    def __init__(self, atom_charges_list, nbondedForce):

        self.atom_charges_list = atom_charges_list
        self.nbondedForce = nbondedForce

    def scale_charges_by_lambda(self, context, lambda_scale):

        for atom in self.atom_charges_list:
            scaled_charge = (1-lambda_scale)*atom[1] + lambda_scale*atom[2]
            (q_i, sig, eps) = self.nbondedForce.getParticleParameters(atom[0])
            self.nbondedForce.setParticleParameters(atom[0], scaled_charge, sig, eps)
        self.nbondedForce.updateParametersInContext(context)

