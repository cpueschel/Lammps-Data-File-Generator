# coding: utf-8
# Distributed under the terms of the MIT License.
# If you make an update that you feel others would need, feel free to make a merge request to the main branch with your update (https://github.com/cpueschel/Lammps-Data-File-Generator.git). 

from __future__ import print_function
import numpy as np
from numpy import linalg as LA
import math
from pymatgen.io.xyz import XYZ
from pymatgen.core import structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.core import Site
from pymatgen import Lattice, Structure, Molecule
import yaml
from pymatgen.io.vasp import Poscar, sets
from string import digits

#Returns known atom types in a list
def known_atom_types():
	types = []
	for each in config['types']:
		types.append(str(each)) 
	return types
#Returns known atom information in a list
def known_atom_types_general():
	types = known_atom_types()
	known_atom_types_general = []
	for every in types:
		charges = config['types'][every]['CHARGES']
		lj_parameters = (config['types'][every]['LJ_PARAMETERS'][0],config['types'][every]['LJ_PARAMETERS'][1])
		mass = config['types'][every]['MASS']
		required_connectors = config['types'][every]['REQUIRED_CONNECTORS']
		known_atom_types_general.append( (every,charges,lj_parameters,mass,required_connectors) )
	return known_atom_types_general
	
#Returns list_bond_types: list of known atom bonds(list of sorted tuples)
def known_bond_types():
	list_bond_types = []
	for every in range(0,len(config['bonds'])):
		i = config['bonds'][every]['i']
		j = config['bonds'][every]['j']
		energy = config['bonds'][every]['ENERGY']
		length = config['bonds'][every]['LENGTH']
		list_bond_types.append( (i,j,energy,length) )
	return list_bond_types

#Returns list_angle_types: list of known angles(list of sorted tuples)
def known_angle_types():
	list_angle_types = []
	for every in range(0,len(config['angles'])):
		i = config['angles'][every]['i']
		j = config['angles'][every]['j']
		k = config['angles'][every]['k']
		energy = config['angles'][every]['ENERGY']
		theta = config['angles'][every]['THETA']
		list_angle_types.append( (i,j,k,energy,theta) )
	return list_angle_types

#Returns known_torsion_types: list of known torsions(list of sorted tuples)
def known_torsion_types():
	known_torsion_types = []
	for every in range(0,len(config['propertorsions'])):
		i = config['propertorsions'][every]['i']
		j = config['propertorsions'][every]['j']
		k = config['propertorsions'][every]['k']
		l = config['propertorsions'][every]['l']
		energy = config['propertorsions'][every]['ENERGY']
		angle = config['propertorsions'][every]['ANGLE']
		multiplicity = config['propertorsions'][every]['MULTIPLICITY']
		known_torsion_types.append( (i,j,k,l,energy,angle,multiplicity) )
	return known_torsion_types

# Max possible bond length for the system with tolerance factor
def max_bond_length():
	max_bond_length = 0
	for each in known_bond_types():
		bond_length = each[3]
		if bond_length > max_bond_length:
				max_bond_length = bond_length
	return float(max_bond_length + max_bond_length*config['bond_length_tolerance_factor'])

# Checks if smaller list is in another list
def is_slice_in_list(s,l):
    len_s = len(s)
    return any(s == l[i:len_s+i] for i in xrange(len(l) - len_s+1))

# Does this list have digits
def hasNumbers(inputLIST):
	for each in inputLIST:
		if any(char.isdigit() for char in each):
			return any(char.isdigit() for char in each)
	return False

# Assigns the order for execution of atom types.
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
		i += 1

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

	# Allow user to verify order		
	print_types = []
	for each in (position_order(types, type_execution_order)):
		print_types.append(types[each])
	print("Order of execution for type assignment:", print_types)
	return position_order(types, type_execution_order)

# Given a atom type (eg. C3) returns number 
def atom_type_NUMBER(type_of_atom):
	i = 0
	for each in known_atom_types_general():
		if type_of_atom == each[0]:
			return i
		i += 1
	return None

# Determines if two sites are bonded together
def site_bonded(site1,site2, bond_length):
	"""
	Requirements for a bond:
		Must have a bond length = required bond length +/- varriance 
		Must be the correct atomic species_string.
		Then it is identified as a "None" type of bond. Since specific atomic types need to be determined based on bonds.
	"""
	site1type = atom_sites[site1].type
	site2type = atom_sites[site2].type
	upper_range = bond_length + bond_length*float(config['bond_length_tolerance_factor'])
	lower_range = bond_length - bond_length*float(config['bond_length_tolerance_factor'])
	bond_types_general = known_bond_types()
	for each in bond_types_general:
		if each[0].translate(None, digits) == site1type and each[1].translate(None, digits) == site2type: 
				if  bond_length >= lower_range and bond_length <= upper_range:
					return True
		if each[0].translate(None, digits) == site2type and each[1].translate(None, digits) == site1type:
				if  bond_length >= lower_range and bond_length <= upper_range:
					return True
	return False				

#Checks for digits within an input string
def contains_digits(s):
    return any(char.isdigit() for char in s)

# Assigns bond types
def bond_type_assignment(siteval):
	atom_number = 0
	for bonded_atom in atom_sites[siteval].bonds:
		bond_type_counter = 0
		for bond_types in known_bond_types():
			if (bond_types[0] == atom_sites[bonded_atom].type and bond_types[1] == atom_sites[siteval].type) or (bond_types[1] == atom_sites[bonded_atom].type and bond_types[0] == atom_sites[siteval].type):
				atom_sites[siteval].bond_types[atom_number] = bond_type_counter
				break
			bond_type_counter += 1
		atom_number += 1

# Removes bonds that are left untyped, we do not have information for this type of bonding behavior.
# This can happen if the bond tolerances are raised too high. 
def bond_untype_removal(i,errors_list):
	indexes = []

	for ii in range(0,len(atom_sites[i].bonds)):
		if atom_sites[i].bond_types[ii] == None:
			#Remove This Bond
			indexes.append(ii)

			if not sorted([atom_sites[i].type,atom_sites[atom_sites[i].bonds[ii]].type]) in errors_list:
				errors_list.append(sorted([atom_sites[i].type,atom_sites[atom_sites[i].bonds[ii]].type]))
				
	if not indexes == []:
		for index in sorted(indexes, reverse=True):		
			del atom_sites[i].bond_types[index]
			del atom_sites[i].bonds[index]
			del atom_sites[i].bond_lengths[index]


	return errors_list
#Assigns atom types			
def type_assignment(siteval,level, order_assign_number):
	atom_types_general,bond_types_general = known_atom_types_general(),known_bond_types()
	i,types, dependencies = 0,[None]*len(atom_types_general),[None]*len(atom_types_general)
	for each in atom_types_general:
		types[i] = each[0]
		dependencies[i] = each[4].values()
		i += 1

	#Check if atomic species is correct for the order asign number.
	if types[order_assign_number].translate(None, digits) == sites[siteval].species_string:
		bonded_atoms,nn_types = atom_sites[siteval].bonds,[]

		# Level 0: (None) Assign higher level types based on basic species, requires no bonds.
		if level == 0:
			atom_sites[siteval].type = types[order_assign_number]

		for each in bonded_atoms:
			#Level 1: Remove type from atom to get species
			if level == 1:
				# Level 1: Assign based on atomic species.
				nn_types.append(atom_sites[each].type.translate(None, digits))
			# Level 2: Assigned based on type attached.
			if level == 2:
				nn_types.append(atom_sites[each].type)
				
		#Case 1: Lists are Identical, assign type
		if sorted(nn_types) == sorted(dependencies[order_assign_number]):
			#Assign new type
			atom_sites[siteval].type = types[order_assign_number]		
		
		#Case 2: Lists are possible, determine if larger than prior diescription
		if sorted(dependencies[order_assign_number]) < sorted(nn_types):
			if contains_digits(atom_sites[siteval].type):
				# Assign highest order assign number
				if len(dependencies[order_assign_number]) >= len(dependencies[atom_type_NUMBER(atom_sites[siteval].type)]):
					atom_sites[siteval].type = types[order_assign_number]
			else:
				atom_sites[siteval].type = types[order_assign_number]
# Counts number of bonds
def count_BONDS():
	bonds,site_number,bonds_added = 0,0,[]
	for each in atom_sites:
		for every in each.bonds:
			if sorted([site_number,every]) not in bonds_added:
				bonds_added.append(sorted([site_number,every]))	
				bonds += 1
		site_number += 1
	return bonds

# Counts number of angles
def count_ANGLES():
	angles = 0
	for each in atom_sites:
		angles += len(each.angles)
	return angles

#Counts number of dihedrals
def count_dihedrals():
	dihedrals = 0
	for each in atom_sites:
		dihedrals += len(each.dihedrals)
	return dihedrals

# Adds a specified number of blank lines
def save_BLANK_LINES(number_blank_lines, f):
	for i in range(0,number_blank_lines):
		print('', file=f)

# Generates the data file
def generate_DATA_FILE():

	"""
	See: http://lammps.sandia.gov/doc/2001/data_format.html
	for an example of data file. 
	"""

	fil = config['filename'].split(".", 1)
	f = open(str(fil[0] + ".data"),'w')

	# LAMMPS Description           (1st line of file)
	print('LAMMPS DATA FILE FOR ' + fil[0], file=f)
	
	save_BLANK_LINES(1, f)

	print(str(len(atom_sites))+ " atoms", file=f)
	print(str(count_BONDS()) + " bonds", file=f)
	print(str(count_ANGLES())+ " angles", file=f)
	num_dihedrals = count_dihedrals()
	if num_dihedrals > 0:
		print(str(num_dihedrals)+ " dihedrals", file=f)
	#print(len(atom_sites)+ " impropers", file=f)

	save_BLANK_LINES(1, f)

	print(str(len(known_atom_types()))+ " atom types", file=f)
	print(str(len(known_bond_types()))+ " bond types", file=f)
	print(str(len(known_angle_types()))+ " angle types", file=f)
	print(str(len(known_torsion_types()))+ " dihedral types", file=f)
	print("0 improper types", file=f)

	save_BLANK_LINES(1, f)

	X = config["lattice_parameters"]["vectors"]['i']
	Y = config["lattice_parameters"]["vectors"]['j']
	Z = config["lattice_parameters"]["vectors"]['k']
	alpha = math.radians(config["lattice_parameters"]["angles"]['alpha'])
	beta = math.radians(config["lattice_parameters"]["angles"]['beta'])
	gamma = math.radians(config["lattice_parameters"]["angles"]['gamma'])
	xlo = 0
	ylo = 0
	zlo = 0
	xhi = LA.norm(X)
	xy = LA.norm(Y)*math.cos(gamma)
	yhi = LA.norm(Y)*math.sin(gamma)
	xz = LA.norm(Z)*math.cos(beta)
	yz = (np.dot(Y,Z) - xy*xz)/yhi #Was yz changed to xz, veryify this was a typo. -Pueschel
	zhi = math.sqrt(math.pow(LA.norm(Z),2) - math.pow(xz,2) - math.pow(yz,2))
	print(str(len(atom_sites))+" "+ str(xlo)+" " + str(xhi) + " xlo xhi", file=f) #-0.5 0.5 xlo xhi       #(for periodic systems this is box size,
	print(str(len(atom_sites))+" "+ str(ylo)+" " + str(yhi) + " ylo yhi", file=f) #-0.5 0.5 ylo yhi       # for non-periodic it is min/max extent of atoms)
	print(str(len(atom_sites))+" "+ str(zlo)+" " + str(zhi) + " zlo zhi", file=f) #-0.5 0.5 zlo zhi       #(do not include this line for 2-d simulations)
      
	save_BLANK_LINES(1, f)

	print("Masses", file=f)
	gen,i = known_atom_types_general(),1
	for each in gen:
		print(str(i) + " " + str(each[3]), file=f)
		i += 1

	# save_BLANK_LINES(1, f)
	# print("Nonbond Coeffs",f)
	#Add this feature if you require it.

	#Bond Coeffs, ENERGY [kcal/mol/A^2], LENGTH [A]
	save_BLANK_LINES(1, f)
	print("Bond Coeffs", file=f) 
	bonds,i = known_bond_types(),1
	for each in bonds:
		print(str(i)+" "+str(each[2])+" "+str(each[3]), file=f)	
		i += 1

	#ANGLES (BENDING), ENERGY [kcal/mol/rad^2], THETA [deg]
	save_BLANK_LINES(1, f)
	print("Angle Coeffs", file=f) 
	angles,i = known_angle_types(),1
	for each in angles:
		print(str(i)+" "+str(each[3])+" "+str(each[4]), file=f)	
		i += 1

	#PROPER TORSIONS, ENERGY [kcal/mol], ANGLE [deg]
	save_BLANK_LINES(1, f)
	print("Dihedral Coeffs", file=f) 
	torsions,i = known_torsion_types(),1
	for each in torsions:
		print(str(i)+" "+str(each[4])+" "+str(each[5])+" "+str(each[6]), file=f)	
		i += 1

	# ATOMS
	save_BLANK_LINES(1, f) 
	print("Atoms", file=f)
	i = 1
	katg = known_atom_types_general()
	molecule_number = 1
	sites = structure.sites
	for each in atom_sites:
		type_of_atom = atom_type_NUMBER(each.type)
		#type_of_atom = 1
		#print("WARNING: ENABLE ABOVE BEFORE RUNNING")
		coords = sites[i-1].coords
		katg_each = katg[type_of_atom]
		#q is atomic charge
		print(str(i)+" "+str(molecule_number)+" "+str(type_of_atom+1)+" "+str(katg_each[1])+" "+str(coords[0])+" "+str(coords[1])+" "+str(coords[2])+" 0 0 0", file=f)
		i += 1

	#Velocities, ADD IF NEEDED
	save_BLANK_LINES(1, f) 
	print("Velocities", file=f)
	i = 1
	for each in atom_sites:
		print(str(i)+" 0 0 0", file=f)
		i += 1

	#Bonds
	save_BLANK_LINES(1, f)
	print("Bonds", file=f) 
	i = 1
	site_number = 1
	bonds_added = []
	for each in atom_sites:
		for every in each.bonds:
			if sorted([site_number,every+1]) not in bonds_added:
				print(str(i)+" "+str(each.bond_type(every)+1)+" "+str(site_number)+" "+str(every+1), file=f)
				bonds_added.append(sorted([site_number,every+1]))	
				i += 1
		site_number = site_number + 1

	#Angles
	save_BLANK_LINES(1, f)
	print("Angles", file=f) 
	i = 1
	site_number = 1
	for each in atom_sites:
		ii = 0
		for every in each.angles:
			print(str(i)+" "+str(each.angle_types[ii]+1)+" "+str(site_number)+" "+str(every[1]+1)+" "+str(every[2]+1), file=f)	
			ii = ii + 1
			i += 1
		site_number = site_number + 1

	#Dihedrals
	save_BLANK_LINES(1, f)
	print("Dihedrals", file=f) 
	i = 1
	site_number = 1
	for each in atom_sites:
		ii = 0
		for every in each.dihedrals:
			print(str(i)+" "+str(each.dihedral_types[ii]+1)+" "+str(site_number)+" "+str(every[1]+1)+" "+str(every[2]+1)+" "+str(every[3]+1), file=f)	
			ii = ii + 1
			i += 1
		site_number = site_number + 1
	
	#Impropers
	save_BLANK_LINES(1, f)
	print("Impropers", file=f) 
	#If you need this functionality, feel free to add this feature and then request a merge to the main branch.

def generate_VDW_DATA_FILE():
	atom_types_general = known_atom_types_general()
	f = open(str("VDW_LAMMPS.data"),'w')
	for i in range(len(atom_types_general)):
		lj_i = atom_types_general[i][2]
		for j in range(len(atom_types_general)):
			lj_j = atom_types_general[j][2]
			mixed_epsilon = np.sqrt(float(lj_i[1])*float(lj_j[1]))
			mixed_sigma = (float(lj_i[0])+float(lj_j[0]))/float(2)
			f.write('pair_coeff '+str(i+1)+' '+str(j+1)+' lj/cut '+str(mixed_epsilon)+' '+str(mixed_sigma)+'\n')		
	f.close()

class StructureSite: 
	# bonds is a list of atoms that the atom is bonded to
	def __init__(self, type):
		self.type = type
		self.bonds = [] #list(set(bonds))
		self.bond_types = []
		self.bond_lengths = []
		self.angles = []
		self.angle_types = []
		self.dihedrals = []
		self.dihedral_types = []

	def add_bond(self, bond,type,bond_length):
		if bond not in self.bonds:
			self.bonds.append(bond)
			self.bond_types.append(type)
			self.bonds = self.bonds
			self.bond_types = self.bond_types
			self.bond_lengths.append(bond_length)

	def bond_type(self, atom2):
		i = 0
		bond_types = self.bond_types
		for each in self.bonds:
			if atom2 == each:
				return bond_types[i] 
			i += 1
		return None

	def bonded(self, atom2):
		for each in self.bonds:
			if atom2 == each:
				return True 
		
	def add_angle(self, atom1, atom2, atom3, angle_type):
		if [atom1,atom2,atom3] not in self.angles:
			if atom_sites[atom1].bonded(atom2):
				if atom_sites[atom2].bonded(atom3):
					if angle_type != None:
						self.angle_types.append(angle_type)
						self.angles.append([atom1, atom2, atom3])

	def check_if_angle(self, atom1, atom2, atom3):
		atoms_list = (atom_sites[atom1].type, atom_sites[atom2].type, atom_sites[atom3].type)
		if atom_sites[atom1].bonded(atom2):
			if atom_sites[atom2].bonded(atom3):
				#Check if angle type is defined in our list of angle types.
				known_angle_type = known_angle_types()
				angle_type_list = []
				for each in known_angle_type:
					angle_type_list.append((each[0],each[1],each[2]))

				count = 0
				for each in angle_type_list:
					if atoms_list == each:
						return (True, count)
					count += 1	
		return (False, None)

	def add_dihedral(self, atom1, atom2, atom3, atom4, dihedral_type):
		#Assuming Only Proper Torsions for now
		if [atom1,atom2,atom3,atom4] not in self.dihedrals:
			if atom_sites[atom1].bonded(atom2):
				if atom_sites[atom2].bonded(atom3):
					if atom_sites[atom3].bonded(atom4):
						self.dihedral_types.append(dihedral_type)
						self.dihedrals.append([atom1, atom2, atom3, atom4])

	def check_if_dihedral(self, atom1, atom2, atom3, atom4):
		#Assuming Only Proper Torsions for now
		atoms_sorted_list = (atom_sites[atom1].type, atom_sites[atom2].type, atom_sites[atom3].type,atom_sites[atom4].type)
		if atom_sites[atom1].bonded(atom2):
			if atom_sites[atom2].bonded(atom3):
				if atom_sites[atom3].bonded(atom4):
					#Check if dihedral type is defined in our list of dihedral types.
					# i,j,k,l for proper torsions
					known_torsion_type = known_torsion_types()
					dihedral_type_list = []
					for each in known_torsion_type:
						dihedral_type_list.append((each[0],each[1],each[2],each[3]))

					count = 0
					for each in dihedral_type_list:
						if atoms_sorted_list == each:
							return (True, count)
						count += 1	
		return (False, None)

#-------------- Import Config File -----------
global config
with open("config", 'r') as ymlfile:
    config = yaml.load(ymlfile)

#-------------- Load in Structure ------------- 
global structure
structure_pos = Poscar.from_file(config['filename'])

structure = structure_pos.structure
sites = structure.sites
position_order_assignment = type_assignment_execution_order()
global atom_sites
atom_sites = [None]*len(sites)
for i in range(len(sites)):
	#This generates the class from above for site_in_Structure
	atom_sites[i] = StructureSite(str(sites[i].species_string))

#-------------- Intiallize Lists ------------- 
atom_types_general = known_atom_types_general()
i,types, dependencies = 0,[None]*len(atom_types_general),[None]*len(atom_types_general)

for each in atom_types_general:
	types[i] = each[0]
	dependencies[i] = each[4].values()
	i += 1
	
#--------- Find Nearest Neighbors -------------
nn_sites = structure.get_all_neighbors(max_bond_length(),include_index=True)
#------------- Bonded Atoms Assignment ----------------
#Iterates through each position assignment for Bond Assignment.
siteval = 0	
for each_site in nn_sites:
	for each in each_site:
		if site_bonded(siteval,each[2], each[1]):
			atom_sites[siteval].add_bond(each[2], None, each[1])
	siteval += 1

#------------- Atom Type Assignment ----------------
level,siteval = 0,0
for order_assign_number in position_order_assignment:
	if dependencies[order_assign_number] == 'None':
		level = 0
	if hasNumbers(dependencies[order_assign_number]):
		level = 2
	else:
		level = 1

	siteval = 0	

	for each_site in nn_sites:
		type_assignment(siteval,level, order_assign_number)
		siteval += 1

#------------- Bond Type Assignment ----------------
for i in range(0,len(nn_sites)):
	bond_type_assignment(i)

#------------- Remove Un-Typed Bonds ---------------
errors = []
for i in range(0,len(nn_sites)):
	errors = bond_untype_removal(i,errors)
for each in errors:		
	print("WARNING: I hope you know what you are doing, there may be undefined bonds within structure:",each[0],"with",each[1])
if not errors == []:
	print("This can happen if the bond tolerances are raised too high.")
	print("Check your tolerances, located in the config file.")
#------------- Dihedral Assignment -----------------
#Iterates through each position assignement for Angle and Dihedrals
siteval = 0	
for each_site in nn_sites:
	for each_nn in each_site:
		for each_second_nn in atom_sites[int(each_nn[2])].bonds: 

			#Checks for Angles (and adds Angles)
			(checked_if_angles, type_angle) = atom_sites[siteval].check_if_angle(siteval,int(each_nn[2]), each_second_nn)
			if checked_if_angles:				
				atom_sites[siteval].add_angle(siteval, int(each_nn[2]), each_second_nn, type_angle)

			#Checks for Diherals (and adds Diherals) 
			for each_third_nn in atom_sites[int(each_second_nn)].bonds: 	
				#Checks for Diherals (and adds Diherals) 

				(checked_if_dihedrals, type_dihedral) = atom_sites[siteval].check_if_dihedral(siteval,int(each_nn[2]),each_second_nn,each_third_nn)
				if checked_if_dihedrals:				
					atom_sites[siteval].add_dihedral(siteval, int(each_nn[2]), each_second_nn, each_third_nn, type_dihedral)				
	siteval += 1
#-------------- Generate Output Files -----------
generate_DATA_FILE()
generate_VDW_DATA_FILE()

#-------------- Generate Output Notes -----------
print("Sucessfully generated data files.")
print(str(len(atom_sites))+ " atoms")
print(str(count_BONDS()) + " bonds")
print(str(count_ANGLES())+ " angles")
num_dihedrals = count_dihedrals()
if num_dihedrals > 0:
	print(str(num_dihedrals)+ " dihedrals")

