import os
import sys
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import lxml.etree as ET

# Given an XML file, creates dictionary mapping atom_classes -> atom_charges
def store_XML_charges(xml_file):
    class_charge_dict = {}
    tree_base = ET.parse(xml_file)
    root = tree_base.getroot()
    for child in root:
        if child.tag == 'NonbondedForce':
            for subchild in child:
                if subchild.tag == 'Atom':
                    atom_charge = subchild.attrib['charge']
                    try:
                        atom_class = subchild.attrib['class']
                    except:
                        atom_class = subchild.attrib['type']
                    try:
                        class_charge_dict[atom_class] = atom_charge
                    except:
                        print ('Something weird happened with the xml.')
                        #pass
    return (class_charge_dict)

# Get initial data for all atoms (name, index, residue.name, residue.index)
# Needs to be passed: topology.atoms()
def store_atom_data(topology_atoms):
    #atom_data = [[i.name, i.index, i.residue.name, i.residue.index] for i in MMsys_neutral.simmd.topology.atoms()]
    atom_data = [[i.name, i.index, i.residue.name, i.residue.index] for i in topology_atoms]
    return (atom_data)


# Given the forcefield._templates object, creates dictionary mapping atom_names -> atom_types
def store_name_type_dict(forcefield_templates):
    name_type_dict = {}
    for template in forcefield_templates:
        for atom in forcefield_templates[template].atoms:
            name_type_dict[atom.name] = atom.type
    return (name_type_dict)

# Combines all other functions within this submodule.
# Returns a list of lists containing all important atom data for the system.
# Data stored as [atom_name, atom_index, residue_name, residue_index, atom_class, atom_charge_neutral, atom_charge_oxidized]
def combiner(topology_atoms, forcefield_templates, neutral_xml_file, oxidized_xml_file):
    atom_data = store_atom_data(topology_atoms)
    name_type_dict = store_name_type_dict(forcefield_templates)
    neutral_class_charge_dict = store_XML_charges(neutral_xml_file)
    oxidized_class_charge_dict = store_XML_charges(oxidized_xml_file)
    for atom in atom_data:
        atom.append(name_type_dict[atom[0]])
        try:
            atom.append(float(neutral_class_charge_dict[name_type_dict[atom[0]]]))
        except:
            pass
            #print ('Could not append neutral charges for atom: %s:%s to atom_data' % (atom[0],atom[1]) )
        try:
            atom.append(float(oxidized_class_charge_dict[name_type_dict[atom[0]]]))
        except:
            pass
            #print ('Could not append oxidized charges for atom: %s:%s to atom_data' % (atom[0],atom[1]) )
    return (atom_data)
      
def store_formatted_dict(topology_atoms, forcefield_templates, neutral_xml_file, oxidized_xml_file):
    atom_data = combiner(topology_atoms, forcefield_templates, neutral_xml_file, oxidized_xml_file)
    formatted_dict = {}
    for atom in atom_data:
        if (len(atom) == 7) and (atom[-1] != atom[-2]):
            if atom[3] not in formatted_dict:
                formatted_dict[atom[3]] = [(atom[1],atom[-2],atom[-1])]
            elif atom[3] in formatted_dict:
                formatted_dict[atom[3]].append((atom[1],atom[-2],atom[-1]))
    return (formatted_dict)






