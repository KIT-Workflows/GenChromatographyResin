import LammpsData as LD
import numpy as np
import copy

def addmol(label, ld, pos=[0.0, 0.0, 0.0]):
	ldn=LD.LammpsData()
	ldn.readFile("data.{}".format(label))
	ldo=LD.LammpsData()

	#copy ref class to output, then add the new info
	#ldo.atoms=ld.atoms[:]
	#ldo.bonds=ld.bonds[:]
	#ldo.angles=ld.angles[:]
	#ldo.dihedrals=ld.dihedrals[:]
	#ldo.impropers=ld.impropers[:]

	#ldo.masses=ld.masses[:]
	#ldo.box=ld.box[:]

	#ldo
	ldo=ld

	offset=[len(ldo.atoms), len(ldo.bonds), len(ldo.angles), len(ldo.dihedrals), len(ldo.impropers)]# if atomIDs are properly set (no omitted indicies)

	offsetT=[ld.natomtypes, ld.nbondtypes, ld.nangletypes, ld.ndihedraltypes, ld.nimpropertypes]

	molmax=np.max([atom[0] for atom in ld.atoms])
	#print("prevmol", molmax)
	for atom in ldn.atoms:
		#mol, type, q, x, y, z
		atom[0]+=molmax
		atom[1]+=offsetT[0]
		atom[3]+=pos[0]
		atom[4]+=pos[1]
		atom[5]+=pos[2]
		ldo.atoms.append(atom)
	for bond in ldn.bonds:
		bond[0]+=offsetT[1]
		bond[1]+=offset[0]
		bond[2]+=offset[0]
		ldo.bonds.append(bond)
	for angle in ldn.angles:
		angle[0]+=offsetT[2]
		angle[1]+=offset[0]
		angle[2]+=offset[0]
		angle[3]+=offset[0]
		ldo.angles.append(angle)
	for dihedral in ldn.dihedrals:
		dihedral[0]+=offsetT[3]
		dihedral[1]+=offset[0]
		dihedral[2]+=offset[0]
		dihedral[3]+=offset[0]
		dihedral[4]+=offset[0]
		ldo.dihedrals.append(dihedral)
	for dihedral in ldn.impropers:
		dihedral[0]+=offsetT[3]
		dihedral[1]+=offset[0]
		dihedral[2]+=offset[0]
		dihedral[3]+=offset[0]
		dihedral[4]+=offset[0]
		ldo.impropers.append(dihedral)
	return ldo, molmax+1, offsetT[0]+1