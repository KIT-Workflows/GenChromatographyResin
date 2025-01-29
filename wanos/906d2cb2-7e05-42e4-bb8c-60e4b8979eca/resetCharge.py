import sys
import LammpsData as LD
from totcharge import totcharge as CH

def readcomponents(fname="poly.config"):
    with open(fname, 'r') as fh:
        fl=fh.readlines()
        numMonomer, numXLinker, numLigand, lSpacer=fl[0].strip().split('are mixed')[0].split(',')
        #print(numMonomer, numXLinker, numLigand)
        numMonomer, typMonomer=numMonomer.split()
        numXLinker, typXLinker=numXLinker.split()
        numLigand, typLigand=numLigand.split('and')[1].split()
        lSpacer = lSpacer.split("with spacer size ")[1]

        numMonomer=int(numMonomer)
        numXLinker=int(numXLinker)
        numLigand=int(numLigand)
        lSpacer=int(lSpacer)

        return numMonomer, typMonomer, numXLinker, typXLinker, numLigand, typLigand, lSpacer

components=readcomponents(sys.argv[1])
if components[5]=="DEAE":
    explicitcharge=components[4]*1.0
elif components[5]=="SO3":
    explicitcharge=components[4]*-1.0
else:
    explicitcharge=0.0

ld=LD.LammpsData()
ld.readFile("fTOPO")

natoms=ld.natoms

netCharge=CH("fTOPO")

for idx, atom in enumerate(ld.atoms):
    ld.atoms[idx][2]-=(netCharge-explicitcharge)/natoms

ld.writeFile("fTOPO_EQ")
print((netCharge-explicitcharge)/natoms)