# coding: utf-8
# Distributed under the terms of the MIT License.
from __future__ import print_function
from pymatgen.io.xyz import XYZ
from pymatgen.core import structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.core import Site
from pymatgen import Lattice, Structure, Molecule
import yaml
from pymatgen.io.vasp import Poscar, sets
from string import digits


def generate_atom_list(items):
	generate_atom_list = []
	for each in items:
		generate_atom_list.append(str(each))
	return generate_atom_list	

def generate_coordinate_list(molecule):
	generate_coordinate_list = [None]*len(molecule.sites)

	gensites = molecule.sites

	for each in range(0,len(molecule.sites)):
		every = molecule.sites[each].coords
		generate_coordinate_list[each] = [every[0], every[1], every[2]]
	return generate_coordinate_list	 

def import_Structure(type, config):
	if type == 'xyz':
		molecule = Molecule.from_file(config['filename'])
		lattice = Lattice([config["lattice_parameters"]["vectors"]['i'],config["lattice_parameters"]["vectors"]["j"],config["lattice_parameters"]["vectors"]["k"]])
		return Structure(lattice, generate_atom_list(molecule.species), generate_coordinate_list(molecule))
	else:
		return None 	

def known_atom_types_Connectors():
	"""
	Returns types_Required_Connectors: required connector types for a type
	"""
	types_Required_Connectors = []
	for each in config['types']:
		types_Required_Connectors.append(config['types'][str(each)]["REQUIRED_CONNECTORS"])
	return types_Required_Connectors

def known_atom_types():
	"""
	Returns known atom types in a list
	"""
	types = []
	for each in config['types']:
		types.append(str(each)) 
	return types

def known_atom_types_general():
	"""
	Returns known atom types in a list
	"""
	types = known_atom_types()
	known_atom_types_general = []
	for every in types:
		charges = config['types'][every]['CHARGES']
		lj_parameters = (config['types'][every]['LJ_PARAMETERS'][0],config['types'][every]['LJ_PARAMETERS'][1])
		mass = config['types'][every]['MASS']
		required_connectors = config['types'][every]['REQUIRED_CONNECTORS']

		known_atom_types_general.append( (every,charges,lj_parameters,mass,required_connectors) )
	return known_atom_types_general

def known_bond_types():
	"""
	Returns known_bond_types: list of known atom bonds(list of sorted tuples)
	"""
	types = known_atom_types()
	list_bond_types = []
	for every in range(0,len(types)):
		i = config['bonds'][every]['i']
		j = config['bonds'][every]['j']
		energy = config['bonds'][every]['ENERGY']
		length = config['bonds'][every]['LENGTH']
		list_bond_types.append( (i,j,energy,length) )
	return list_bond_types

def known_angle_types():
	"""
	Returns known_angle_types: list of known angles(list of sorted tuples)
	"""
	types = []
	for each in config['angles']:
		types.append(str(each)) 

	list_angle_types = []
	for every in range(0,len(types)):
		i = config['angles'][every]['i']
		j = config['angles'][every]['j']
		k = config['angles'][every]['k']
		energy = config['angles'][every]['ENERGY']
		theta = config['angles'][every]['THETA']
		list_angle_types.append( (i,j,k,energy,theta) )
	return list_angle_types

def known_torsion_types():
	"""
	Returns known_torsion_types: list of known torsions(list of sorted tuples)
	"""
	types = []
	for each in config['propertorsions']:
		types.append(str(each)) 

	known_torsion_types = []
	for every in range(0,len(types)):
		i = config['propertorsions'][every]['i']
		j = config['propertorsions'][every]['j']
		k = config['propertorsions'][every]['k']
		l = config['propertorsions'][every]['l']
		energy = config['propertorsions'][every]['ENERGY']
		angle = config['propertorsions'][every]['ANGLE']
		multiplicity = config['propertorsions'][every]['MULTIPLICITY']
		known_torsion_types.append( (i,j,k,l,energy,angle,multiplicity) )
	return known_torsion_types

def max_bond_length():
	max_bond_length = 0
	for each in known_bond_types():
		bond_length = each[3]
		if bond_length > max_bond_length:
				max_bond_length = bond_length
	return float(max_bond_length + max_bond_length*config['bond_length_tolerance_factor'])


def is_slice_in_list(s,l):
    len_s = len(s) #so we don't recompute length of s on every iteration
    return any(s == l[i:len_s+i] for i in xrange(len(l) - len_s+1))

def hasNumbers(inputLIST):
	for each in inputLIST:
		if any(char.isdigit() for char in each):
			return any(char.isdigit() for char in each)
	return False

def type_assignment_execution_order():
	"""
		Returns an array for the order of execution of the type assignment 
		eg. C1 is assigned before C2
		Order:
		None types assigned first
		Those with just "X" are assigned 2nd, X is assumed to be equal to X1
		Those with "C2" assigned second, etc.
	"""
	def position_order(types, order):
		position_order = []
		for each in order:
			for i in range(0,len(types)):
				if types[i] == each:
					position_order.append(i)
		return position_order

	#intialize
	atom_types_general = known_atom_types_general()
	i, types, dependencies, type_execution_order = 0, [None]*len(atom_types_general), [None]*len(atom_types_general), []
	
	for each in atom_types_general:
		types[i] = each[0]
		dependencies[i] = each[4].values()
		i = i + 1

	# #Intial Sorting
	# types = sorted(types, key=str.lower)

	#Change to order in types at some point, for readability keep variable # for now.
	#1 Check for items without dependencies.
	switch = 0
	while len(type_execution_order) < len(types):
		if switch == 0:
			for i in range(0,len(types)):
				if 'None' in dependencies[i] and types[i] not in type_execution_order:
					type_execution_order.append(types[i])

		if switch == 1:
			for i in range(0,len(types)):
				if 'None' in dependencies[i]:
					continue
				if is_slice_in_list(dependencies[i],type_execution_order) and types[i] not in type_execution_order:
					type_execution_order.append(types[i])
		if switch == 2:
			for i in range(0,len(types)):
				if 'None' in dependencies[i]:
					continue
				if is_slice_in_list(dependencies[i],type_execution_order) and types[i] not in type_execution_order:
					continue
				if not hasNumbers(dependencies[i]) and types[i] not in type_execution_order:
					type_execution_order.append(types[i])

		if switch == 3:
			for i in range(0,len(types)):
				if types[i] not in type_execution_order:
					type_execution_order.append(types[i])

		switch = switch + 1
		if switch > 10:
			print("Arguments for input not properly ordered")
			break
	return position_order(types, type_execution_order)

def atom_type_NUMBER(type_of_atom):
	i = 0
	for each in known_atom_types_general():
		if type_of_atom == each[0]:
			return i
		i = i + 1
	return None

def site_bonded(site1,site2, bond_length):
	"""
	Requirements for a bond:
	
	Must have a bond length = required bond length +/- varriance 
	Must be the correct types.
	"""
	site1type = atom_sites[site1].type
	site2type = atom_sites[site2].type
	upper_range = bond_length + bond_length*float(config['bond_length_tolerance_factor'])
	lower_range = bond_length - bond_length*float(config['bond_length_tolerance_factor'])
	count = 0
	for each in known_bond_types(): # i=0,j=1,length=3
		if each[0].translate(None, digits) == site1type and each[1].translate(None, digits) == site2type: 
			if  bond_length >= lower_range and bond_length <= upper_range:
				return (True,count)
		if each[0].translate(None, digits) == site2type and each[1].translate(None, digits) == site1type:
			if  bond_length >= lower_range and bond_length <= upper_range:
				return (True,count)
		count = count + 1
	return (False, None)


def redefine_site_type():
	"""
	Once the bonds have been updated, we check again to see if species are correctly identified.
	For each site:
		Identfy type.
			Identify types bonded to.
			See if these correspond to a the required connectors.
	"""
	atom_types_general = known_atom_types_general()
	i,types, dependency = 0,[None]*len(atom_types_general),[None]*len(atom_types_general)
	for each in atom_types_general:
		types[i] = each[0]
		dependencies[i] = each[4].values()
		i = i+ 1

	for i in range(len(structure.sites)):
		type_atom = atom_sites[i].type
		bonds_current = atom_sites[i].bonds
		bond_types_current = atom_sites[i].bond_types
		nn_types, nn_types_new= [],[]

		for each in bonds_current:
			nn_types.append(str(sites[each].species_string))
			nn_types_new.append(str(atom_sites[each].type))
		nn_types = sorted(nn_types)
		nn_types_new = sorted(nn_types_new)

		for ii in range(len(atom_types_general)):
			if None not in dependencies[ii]:
				if sorted(dependencies[ii]) == nn_types:
					atom_sites[i].type = types[ii]
				if sorted(dependencies[ii]) == nn_types_new:
					atom_sites[i].type = types[ii]
			
def count_BONDS():
	bonds = 0
	for each in atom_sites:
		bonds = bonds + len(each.bonds)
	return bonds/2

def save_BLANK_LINES(number_blank_lines, f):
	for i in range(0,number_blank_lines):
		print('', file=f)

def generate_DATA_FILE():

"""
	See: http://lammps.sandia.gov/doc/2001/data_format.html
	for an example of data file. 
"""

	fil = config['filename'].split(".", 1)
	f = open(str(fil[0] + ".data"),'w')

	######### LAMMPS Description           (1st line of file)
	print('LAMMPS DATA FILE FOR ' + fil[0], file=f)
	
	save_BLANK_LINES(1, f)

	print(str(len(atom_sites))+ " atoms", file=f)
	print(str(count_BONDS()) + " bonds", file=f)
	#print(len(atom_sites)+ " angles",f)
	#print(len(atom_sites)+ " dihedrals",f)
	#print(len(atom_sites)+ " impropers",f)

	save_BLANK_LINES(1, f)

	print(str(len(known_atom_types()))+ " atom types", file=f)
	print(str(len(known_bond_types()))+ " bond types", file=f)
	print(str(len(known_angle_types()))+ " angle types", file=f)
	print(str(len(known_torsion_types()))+ " dihedral types", file=f)
	print("0 improper types", file=f)

	save_BLANK_LINES(1, f)

	#print(Lattice.lengths_and_angles(structure))
	#print(len(atom_sites)+ " xlo xhi",f) #-0.5 0.5 xlo xhi       #(for periodic systems this is box size,
	#print(len(atom_sites)+ " ylo yhi",f) #-0.5 0.5 ylo yhi       # for non-periodic it is min/max extent of atoms)
	#print(len(atom_sites)+ " zlo zhi",f) #-0.5 0.5 zlo zhi       #(do not include this line for 2-d simulations)
      
	save_BLANK_LINES(1, f)
	print("Masses", file=f)
	gen,i = known_atom_types_general(),0
	for each in gen:
		print(str(i) + " " + str(each[3]), file=f)
		i = i + 1

	# save_BLANK_LINES(1, f)
	# print("Nonbond Coeffs",f)
	#Add this feature if you require it.

	#Bond Coeffs, ENERGY [kcal/mol/A^2], LENGTH [A]
	save_BLANK_LINES(1, f)
	print("Bond Coeffs", file=f) 
	bonds,i = known_bond_types(),0
	for each in bonds:
		print(str(i)+" "+str(each[2])+" "+str(each[3]), file=f)	
		i = i + 1

	#ANGLES (BENDING), ENERGY [kcal/mol/rad^2], THETA [deg]
	save_BLANK_LINES(1, f)
	print("Angle Coeffs", file=f) 
	angles,i = known_angle_types(),0
	for each in angles:
		print(str(i)+" "+str(each[3])+" "+str(each[4]), file=f)	
		i = i + 1

	#PROPER TORSIONS, ENERGY [kcal/mol], ANGLE [deg]
	save_BLANK_LINES(1, f)
	print("Dihedral Coeffs", file=f) 
	torsions,i = known_torsion_types(),0
	for each in torsions:
		print(str(i)+" "+str(each[4])+" "+str(each[5])+" "+str(each[6]), file=f)	
		i = i + 1

	#Add additional properties as needed.

	save_BLANK_LINES(1, f) 
	print("Atoms", file=f)
	i = 0
	katg = known_atom_types_general()
	molecule_number = 1
	sites = structure.sites
	for each in atom_sites:
		# type_of_atom = atom_type_NUMBER(each.type)
		type_of_atom = 1
		print("WARNING: ENABLE ABOVE BEFORE RUNNING")
		coords = sites[i].coords
		katg_each = katg[type_of_atom]
		#q is atomic charge
		print(str(i)+" "+str(molecule_number)+" "+str(type_of_atom)+" "+str(katg_each[1])+" "+str(coords[0])+" "+str(coords[1])+" "+str(coords[2])+" 0 0 0", file=f)
		i = i + 1

	#Add a velocity property to the site class if you need to add this feature.
	save_BLANK_LINES(1, f) 
	print("Velocities", file=f)
	i = 0
	for each in atom_sites:
		print(str(i)+" 0 0 0", file=f)
		i = i + 1

	#Bonds
	save_BLANK_LINES(1, f)
	print("Bonds", file=f) 
	i = 0
	for each in atom_sites:
		print(str(i)+" "+str(each[4])+" "+str(each[5])+" "+str(each[6]), file=f)	
		i = i + 1

		list_bond_types.append( (i,j,energy,length) )

class site_in_Structure: 
	# bonds is a list of atoms that the atom is bonded to
	# eg. (site)
	def __init__(self, type):
		self.type = type
		self.bonds = [] #list(set(bonds))
		self.bond_types = []
		self.bond_lengths = []
		self.angles = []

	def add_bond(self, bond,type,bond_length):
		if bond not in self.bonds:
			self.bonds.append(bond)
			self.bond_types.append(type)
			self.bonds = self.bonds
			self.bond_types = self.bond_types
			self.bond_lengths.append(bond_length)
	def bonded(self, atom2):
		for each in self.bonds:
			if atom2 in each:
				return True 
	print("Add angle is under construction")			
	def add_angle(self, atom1, atom2, atom3):
		if [atom1,atom2,atom3] not in self.angles:
			if atom_sites[atom1].bonded(atom2):
				if atom_sites[atom2].bonded(atom3):
					#Find angle type.
					self.angle_type = angle_type
					self.angles.append()






## STEPS TO RUN ##

#1.
########### Import Config File ###################
global config
with open("config", 'r') as ymlfile:
    config = yaml.load(ymlfile)


#-------------- Load in Structure ------------- 
global structure
structure_pos = Poscar.from_file("CHA_ZIF8.vasp")#import_Structure(config['structure_type'], config)
structure = structure_pos.structure
sites = structure.sites
position_order_assignment = type_assignment_execution_order()
global atom_sites
atom_sites = [None]*len(sites)
for i in range(len(sites)):
	#This generates the class from above for site_in_Structure

	atom_sites[i] = site_in_Structure(str(sites[i].species_string))

#Intiallize Lists
#print "Distance of 0 and 39 index: " + str(structure.get_distance(0,39))
atom_types_general = known_atom_types_general()
i,types, dependencies = 0,[None]*len(atom_types_general),[None]*len(atom_types_general)

for each in atom_types_general:
	types[i] = each[0]
	dependencies[i] = each[4].values()
	# if 'None' in dependencies[i]:
	# 	print "None"
	i = i + 1




#List of all sites
nn_sites = (structure.get_neighbors(sites[0], max_bond_length(), include_index=True),structure.get_neighbors(sites[1], max_bond_length(), include_index=True))
# items = structure.get_all_neighbors(config["bond_length_tolerance_factor"],include_index=True)

"""
This is to debug use other version for all nn geneations.
Need to iterate through "level" by level. 

None-type are already assigned with proper site types.
Next levels can be interdependant and need to be iterated through one at a time.

"""
for order_assign_number in position_order_assignment:
	siteval = 0	
	for each_site in nn_sites:
		site_species = sites[siteval].species_string
		for each_nn in each_site:	
			bond_length = each_nn[1]
			(checked_if_bonds, type_bond) = site_bonded(siteval,each_nn[2], bond_length)
			if checked_if_bonds:#each_site, site_species, siteval):
				atom_sites[siteval].add_bond(int(each_nn[2]), type_bond, bond_length)
		siteval = siteval + 1
	redefine_site_type()





#Check input structure	
#structure.to(filename="POSCAR")


########### Generate Output Files ###################
generate_DATA_FILE()

########### Generate Facts ###################
