import LammpsDataFF as LD
import sys

ff1=LD.LammpsDataFF()
monomertyp=sys.argv[1]
ff1.readFile("FF.{}2".format(monomertyp))

ff2=LD.LammpsDataFF()
xlinktyp=sys.argv[2]
ff2.readFile("FF.{}half".format(xlinktyp))

ff3=LD.LammpsDataFF()
ligandtyp=sys.argv[3]
ff3.readFile("FF.{}".format(ligandtyp))

ffo=LD.LammpsDataFF()
ATOMS1=[]
for item in ff1.atoms:
    ATOMS1.append(" ".join(item))
BONDS1=[]
for item in ff1.bonds:
    BONDS1.append(" ".join(item))
ANGLES1=[]
for item in ff1.angles:
    ANGLES1.append(" ".join(item))
DIHEDRALS1=[]
for item in ff1.dihedrals:
    DIHEDRALS1.append(" ".join(item))
IMPROPERS1=[]
for item in ff1.impropers:
    IMPROPERS1.append(" ".join(item))
ATOMS2=[]
for item in ff2.atoms:
    ATOMS2.append(" ".join(item))
BONDS2=[]
for item in ff2.bonds:
    BONDS2.append(" ".join(item))
ANGLES2=[]
for item in ff2.angles:
    ANGLES2.append(" ".join(item))
DIHEDRALS2=[]
for item in ff2.dihedrals:
    DIHEDRALS2.append(" ".join(item))
IMPROPERS2=[]
for item in ff2.impropers:
    IMPROPERS2.append(" ".join(item))
ATOMS3=[]
for item in ff3.atoms:
    ATOMS3.append(" ".join(item))
BONDS3=[]
for item in ff3.bonds:
    BONDS3.append(" ".join(item))
ANGLES3=[]
for item in ff3.angles:
    ANGLES3.append(" ".join(item))
DIHEDRALS3=[]
for item in ff3.dihedrals:
    DIHEDRALS3.append(" ".join(item))
IMPROPERS3=[]
for item in ff3.impropers:
    IMPROPERS3.append(" ".join(item))

for item in ff1.atoms:
    ffo.atoms.append(" ".join(item))
for item in ff1.bonds:
    ffo.bonds.append(" ".join(item))
for item in ff1.angles:
    ffo.angles.append(" ".join(item))
for item in ff1.dihedrals:
    ffo.dihedrals.append(" ".join(item))
for item in ff1.impropers:
    ffo.impropers.append(" ".join(item))

for item in ff2.atoms:
    ffo.atoms.append(" ".join(item))
for item in ff2.bonds:
    ffo.bonds.append(" ".join(item))
for item in ff2.angles:
    ffo.angles.append(" ".join(item))
for item in ff2.dihedrals:
    ffo.dihedrals.append(" ".join(item))
for item in ff2.impropers:
    ffo.impropers.append(" ".join(item))

for item in ff3.atoms:
    ffo.atoms.append(" ".join(item))
for item in ff3.bonds:
    ffo.bonds.append(" ".join(item))
for item in ff3.angles:
    ffo.angles.append(" ".join(item))
for item in ff3.dihedrals:
    ffo.dihedrals.append(" ".join(item))
for item in ff3.impropers:
    ffo.impropers.append(" ".join(item))

#print(ffo.atoms)
ffo.writeFile("FFparams.combined")
