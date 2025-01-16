import LammpsDataFF as LDff
import os
import sys

lSpacer=int(sys.argv[1])

fname="fTOPO_EQ"

ldff=LDff.LammpsDataFF()
ldff.readFile("FFparams.combined")

if not lSpacer:
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    
    ldff.angles.append(['74.000', '108.900'])
else:
    ldff.masses.append(12.011)
    ldff.masses.append(1.008)
    ldff.masses.append(12.011)
    ldff.masses.append(12.011)
    ldff.masses.append(12.011)    
    ldff.atoms.append(['0.066', '3.500']) #generic carbon
    ldff.atoms.append(['0.030', '2.500']) #generic hydrogen
    ldff.atoms.append(['0.066', '3.500']) #second carbon to the left
    ldff.atoms.append(['0.066', '3.500']) #first carbon to the left
    ldff.atoms.append(['0.066', '3.500']) #right-most carbon

    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    ldff.bonds.append(['268.0000', '1.5290'])
    
    ldff.bonds.append(['340.0000', '1.0900'])
    ldff.bonds.append(['268.0000', '1.5290'])
    
    ldff.angles.append(['37.5000', '110.7'])
    ldff.angles.append(['33.0000', '107.8'])
    ldff.angles.append(['58.3500', '112.7'])
    ldff.angles.append(['74.000', '108.900'])

ldff.writeFile("FFparams.update")

with open("FFparams.update", 'r') as fr:
        fl=fr.readlines()

with open(fname, 'a') as f:
    f.write("\n")
    for line in fl:
        f.write(" ".join(line.split())+"\n")
        
os.system('mv {} fPoly_FF'.format(fname))