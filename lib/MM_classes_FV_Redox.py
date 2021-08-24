from sys import *
#******** import parent class
from FV_shared.lib.MM_classes_FV import *
#******** import thermodynamic integration classes/methods
from thermodynamic_integration import *

#*************************** README  **************************************
# This is a child MM system class for Fixed-Voltage simulations involving 
# oxidation-reduction reactions. Inherits from parent MM_fixedVoltage
# class in MM_classes_FV.py
#
#**************************************************************************
# note that using super().__init__() in all child classes ensures that
# MM_base.__init__() isn't doubly called.
#**************************************************************************
class MM_FixedVoltage_Redox(MM_FixedVoltage):
    # required input: 1) list of pdb files, 2) list of residue xml files, 3) list of force field xml files.
    def __init__( self , pdb_list , residue_xml_list , ff_xml_list , **kwargs  ):

        # constructor runs through MM_QMMM and MM_FixedVoltage ...
        super().__init__( pdb_list , residue_xml_list , ff_xml_list , **kwargs )

        # reading inputs from **kwargs

    # only restrains atoms in the z-direction but can be easily generalized to x and y
    def set_position_constraint( self, atom_name, kz, z0 ):
        ZForce = CustomExternalForce("0.5*kz*periodicdistance(x,y,z,x,y,z0)^2")
        for res in self.modeller.topology.residues():
            for atom in res._atoms:
                if atom.name in atom_name:
                    ZForce.addParticle(atom.index)
        ZForce.addGlobalParameter('z0', z0)
        ZForce.addGlobalParameter('kz', kz)
        self.system.addForce(ZForce)

    #**************************************
    # trivial scaling algorithm for updating charges during iterative solver
    #**************************************
    def update_charge_iteration( self , q_old , q_new , lambda_scale = 1.0 ):
        q_update = (1.0 - lambda_scale) * q_old + lambda_scale * q_new
        return q_update
 

    #***************************************
    # this generates exclusions for intra-electrode interactions,
    #
    # if flag_SAPT_FF_exclusions=True, then will also set exclusions for SAPT-FF force field...
    #***************************************
    def generate_exclusions(self, water_name = 'HOH', flag_hybrid_water_model = False ,  flag_SAPT_FF_exclusions = True , **kwargs ):

        
        # if Cathode/Anode are present, setup exclusions ...
        if self.Cathode is not None and self.Anode is not None :

            if 'surfactant_exclude_index' in kwargs :
                surfactant_exclude_index = kwargs['surfactant_exclude_index']
                self.Cathode.electrode_extra_exclusions.append( surfactant_exclude_index )


        if 'surfactant_zero_charge_index' in kwargs :
            surfactant_zero_charge_index = kwargs['surfactant_zero_charge_index']
            #zeroing charges on specified atoms
            charge_sum = 0.0
            for i in surfactant_zero_charge_index:
                (q_i, sig, eps) = self.nbondedForce.getParticleParameters(i)
                self.nbondedForce.setParticleParameters(i, 0.0, sig , eps)
                charge_sum += q_i._value
            self.nbondedForce.updateParametersInContext(self.simmd.context)

            for i in surfactant_zero_charge_index:
                (q_i, sig, eps) = self.nbondedForce.getParticleParameters(i)
            if abs(charge_sum) > 10e-6:
                print ('zeroing charges has changed neutrality of molecule, this isnt allowed!')
                print ('charge_sum', charge_sum)
                sys.exit()

        super().generate_exclusions(water_name = water_name, flag_hybrid_water_model = flag_hybrid_water_model ,  flag_SAPT_FF_exclusions = flag_SAPT_FF_exclusions) # , **kwargs )

        # now reinitialize to make sure changes are stored in context
        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(positions)


    #*******************************************
    #   this method performs MC/MD steps for moving electrode sheets to
    #   equilibrate the density of the electrolyte
    #
    #         !!! assumes object self.MC exists !!!!
    #
    #   before calling this method, make sure to initialize MC parameters by
    #   construcing self.MC object of class MC_parameters(object):
    #*******************************************
    def MC_Barostat_step( self ):

        # inner functions ...
        def metropolis(pecomp):
            if pecomp < 0.0 * self.MC.RT :
                return True
            elif (random.uniform(0.0,1.0) < numpy.exp(-pecomp/self.MC.RT)):
                return True

        def intra_molecular_vectors( residue_object , pos_ref, positions_array ):
            intra_vec = []
            # loop over atoms in residue
            for atom in residue_object._atoms:
                pos_res_i = positions_array[atom.index]
                vec_i = pos_res_i - pos_ref
                intra_vec.append(numpy.asarray(vec_i))
            return numpy.asarray(intra_vec)


        self.MC.ntrials += 1

        # ************** normal MD steps ****************
        self.simmd.step(self.MC.barofreq)
        print (self.simmd.context.getParameter('z0'))

        # get final positions
        state = self.simmd.context.getState(getEnergy=True, getPositions=True)
        positions = state.getPositions().value_in_unit(nanometer)

        # energy before move
        oldE = state.getPotentialEnergy()
        # store positions
        oldpos = numpy.asarray(positions)
        newpos = numpy.asarray(positions)

        # now generate trial move
        # randomly choose move distance ... (-1 , 1) Angstrom for now...
        deltalen = self.MC.shiftscale*(random.uniform(0, 1) * 2 - 1)

        # need a reference point for scaling relative positions.  Choose the stationary electrode for this.
        reference_atom_index = -1

        # Currently, can only move Anode, because we might have other conductors on Cathode (Buckyballs, nanotubes)
        # could easily generalize this to move other conductors as well, but haven't yet...
        self.MC.electrode_move = 'Cathode'
        if self.MC.electrode_move == "Anode" :
            reference_atom_index = self.Cathode.electrode_atoms[0].atom_index # since we are moving Anode, choose stationary cathode as reference...
            # move Anode
            for atom in self.Anode.electrode_atoms:
                newpos[atom.atom_index,2] += deltalen
            # now see if we need to move any other chains in Anode
            if len( self.Anode.electrode_extra_exclusions ) > 0 :
                for extra_anode_sheet in self.Anode.electrode_extra_exclusions :
                    for index in extra_anode_sheet:
                        newpos[index,2] += deltalen
        else:
            reference_atom_index = self.Anode.electrode_atoms[0].atom_index # since we are moving Anode, choose stationary cathode as reference...
            # move Cathode
            for atom in self.Cathode.electrode_atoms:
                newpos[atom.atom_index,2] -= deltalen
            # now see if we need to move any other chains in Anode
            if len( self.Cathode.electrode_extra_exclusions ) > 0 :
                for extra_cathode_sheet in self.Cathode.electrode_extra_exclusions :
                    for index in extra_cathode_sheet:
                        newpos[index,2] -= deltalen

        z0 = self.simmd.context.getParameter('z0')
        z0 -= deltalen
        self.simmd.context.setParameter('z0', z0)

        Lcell_old = self.Lcell
        Lcell_new = Lcell_old + deltalen


        N_electrolyte_mol=0
        # now loop over electrolyte molecules and move their COM
        for res in self.electrolyte_residues:
            N_electrolyte_mol += 1
            # use first atom in residue as reference...
            for atom in res._atoms:
                pos_ref = newpos[atom.index]
                break
            # get relative coordinates
            intra_vec = intra_molecular_vectors( res , pos_ref, newpos )

            # convert to coordinate system relative to stationary electrode ...
            pos_ref[2] = pos_ref[2] - newpos[reference_atom_index,2]
            # now scale this reference position by ratio of change in electrochemical cell volume
            pos_ref[2] = pos_ref[2] * Lcell_new / Lcell_old
            # now back to global coordinate system ...
            pos_ref[2] = pos_ref[2] + newpos[reference_atom_index,2]
            # now update positions of all atoms in this molecule
            index_in_molecule = 0
            for atom in res._atoms:
                newpos[atom.index] = pos_ref + intra_vec[index_in_molecule]
                index_in_molecule+=1


        #  Energy of trial move
        self.simmd.context.setPositions(newpos)
        statenew = self.simmd.context.getState(getEnergy=True,getPositions=True)
        newE = statenew.getPotentialEnergy()

        w = newE-oldE + self.MC.pressure*(deltalen * nanometer) - N_electrolyte_mol * self.MC.RT * numpy.log(Lcell_new/Lcell_old)
        if metropolis(w):
            self.MC.naccept += 1
            # move is accepted, update electrochemical cell parameters...
            boxVecs = self.simmd.topology.getPeriodicBoxVectors()
            positions = statenew.getPositions()
            self.set_electrochemical_cell_parameters( positions, boxVecs )
        else:
            # move is rejected, revert to old positions ...
            self.simmd.context.setPositions(oldpos)
            z0 += deltalen
            self.simmd.context.setParameter('z0', z0)

        if self.MC.ntrials > 50 :
            #print(" After 50 more MC steps ...")
            #print("dE, exp(-dE/RT) ", w, numpy.exp(-w/ self.MC.RT))
            #print("Accept ratio for last 50 MC moves", self.MC.naccept / self.MC.ntrials)
            if (self.MC.naccept < 0.25*self.MC.ntrials) :
                self.MC.shiftscale /= 1.1
            elif self.MC.naccept > 0.75*self.MC.ntrials :
                self.MC.shiftscale *= 1.1
            # reset ...
            self.MC.ntrials = 0
            self.MC.naccept = 0



