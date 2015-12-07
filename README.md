#LAMMPS Data File Generator
Generate data file for lammps using force field and a POSCAR (VASP) like input file.

## Note: Currently adding functionality for angles and proper torsions (dihedrals)

##1. What does it do?
Quickly generate a LAMMPS data file to define a force field for a periodic system and Van Der Waals interactional parameters data file. 
This includes:
* Bonds
* Angles
* Proper Torsions (dihedrals)

##2. Required Dependencies
Make sure you have installed pymatgen. This can be done:
    `easy_install pymatgen`
or 
    `pip install pymatgen`
Alternatively, go to their [website](http://pymatgen.org/) for installation instructions.

##3. Quick-start Guide
Obtain all of the necessary input files; 
1. Config containing proper atom, bond, angle, dihedral, and filename.
2. Molecule/System defined in a POSCAR (VASP). Additional support for just xyz is being implemented.
3. Ensure you have pymatgen and python 2.7+

To run simply have the input files in the same folder;
`python LDFG.py`
