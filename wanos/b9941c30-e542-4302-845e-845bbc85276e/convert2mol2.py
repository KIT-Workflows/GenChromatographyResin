import numpy as np
import LammpsDataDict as LD
import sys

infile=sys.argv[1]

ldW = LD.LammpsData()
ldW.readFile(infile)

def findelement(mass):
    refmass={'C':12.0107, 'H': 1.00784, 'N': 14.0067, 'O': 15.999, 'S': 32.065}
    found=False
    for element in refmass:
        if np.around(mass, decimals=2)==np.around(refmass[element], decimals=2):
            found=True
            return element
    if not found: #use carbon as a place holder
        return "C"

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

def wrap(x, y, z, ld):
    boxx=ld.box[1][0]-ld.box[0][0]
    boxy=ld.box[1][1]-ld.box[0][1]
    boxz=ld.box[1][2]-ld.box[0][2]
    xlo=ld.box[0][0]
    ylo=ld.box[0][1]
    zlo=ld.box[0][2]

    x=x#+int(-x/boxx+0.5)*boxx-xlo
    y=y#+int(-y/boxy+0.5)*boxy-ylo
    z=z#+int(-z/boxz+0.5)*boxz-zlo

    return x, y, z

def chkatyp(atype, types, masses):
    ele=findelement(masses[atype-1])
    #print(ele)

    if ele=="C":
        if types==["HEMA", "EGDMA", "TRP"]:
            if atype in [7, 13, 44, 72]: #7, 13:in HEMA2:as C=0, 44:in EDMA as C=O, 72:in TRP as C=O
                return "C.2"
            elif atype in [62, 63, 66, 67, 68, 69, 70, 71]:
                return "C.ar"
            else:
                return "C.3"
        elif types==["DHPMA", "EGDMA", "TRP"]:
            if atype in [7, 13, 52, 80]: #7, 13:in HEMA2:as C=0, 44:in EDMA as C=O, 72:in TRP as C=O
                return "C.2"
            elif atype in [70, 71, 74, 75, 76, 77, 78, 79]:
                return "C.ar"
            else:
                return "C.3"
        elif types==["HEMA", "EGDMA", "DEAE"]:
            if atype in [7, 13, 44]: #7, 13:in HEMA2:as C=0, 44:in EDMA as C=O
                return "C.2"
            else:
                return "C.3"
        elif types==["DHPMA", "EGDMA", "DEAE"]:
            if atype in [7, 13, 52]:
                return "C.2"
            else:
                return "C.3"
        elif types==["HEMA", "EGDMA", "SO3"]:
            if atype in [7, 13, 44]: #7, 13:in HEMA2:as C=0, 44:in EDMA as C=O
                return "C.2"
            else:
                return "C.3"                
        elif types==["DHPMA", "EGDMA", "SO3"]:
            if atype in [7, 13, 52]: 
                return "C.2"
            else:
                return "C.3"

    elif ele=="N":
        if types==["HEMA", "EGDMA", "TRP"]:
            if atype in [64]:
                return "N.ar"
            elif atype in [59]:
                return "N.am"
        elif types==["DHPMA", "EGDMA", "TRP"]:
            if atype in [72]:
                return "N.ar"
            elif atype in [67]:
                return "N.am"
        elif types==["HEMA", "EGDMA", "DEAE"]:
            if atype in [69]:
                return "N.4"
        elif types==["DHPMA", "EGDMA", "DEAE"]:
            if atype in [69]:
                return "N.4"


    elif ele=="O":
        if types==["HEMA", "EGDMA", "TRP"]:
            if atype in [8, 14, 45, 73]: #7, 14:in HEMA2:as C=O, 45:in EDMA as C=O, 73:in TRP as C=O
                return "O.co2"
            else:
                return "O.3"
        elif types==["DHPMA", "EGDMA", "TRP"]:
            if atype in [8, 53, 81]: #7, 14:in HEMA2:as C=O, 53:in EDMA as C=O, 73:in TRP as C=O
                return "O.co2"
            else:
                return "O.3"
        elif types==["HEMA", "EGDMA", "DEAE"]:
            if atype in [8, 14, 45]: #7, 14:in HEMA2:as C=O, 45:in EDMA as C=O, 73:in TRP as C=O
                return "O.co2"
            else:
                return "O.3"
        elif types==["DHPMA", "EGDMA", "DEAE"]:
            if atype in [8, 53]: #8 in DHPMA as C=O 53:in EDMA as C=O
                return "O.co2"
            else:
                return "O.3"
        elif types==["HEMA", "EGDMA", "SO3"]:
            if atype in [8, 14, 45, 60, 61, 62]: #8 in DHPMA as C=O 53:in EDMA as C=O
                return "O.co2"
            else:
                return "O.3"
        elif types==["DHPMA", "EGDMA", "SO3"]:
            if atype in [8, 53, 68, 69, 70]: #8 in DHPMA as C=O 53:in EDMA as C=O
                return "O.co2"
            else:
                return "O.3"

    elif ele=="H":
        return "H"
    elif ele=="S":
        return "S.2"

car=[]
oxy=[]
with open(sys.argv[2], 'r') as fr:
    fl=fr.readlines()
    for line in fl:
        carc, oxyc=line.split()
        car.append(int(carc))
        oxy.append(int(oxyc))

numMonomer, typMonomer, numXLinker, typXLinker, numLigand, typLigand, lSpacer=readcomponents(sys.argv[3])
types=[typMonomer, typXLinker, typLigand]

with open(infile+"_conv.mol2", 'w') as f:
    
    elementIDs=[]

    for mass in ldW.masses:
        """
         1:12.011000:# C
         2:1.008000:# H
         3:14.007000:# N
         4:15.999000:# O
        """
        elementIDs.append(findelement(mass))
    count=1
    atompos={}
    nmol=0
    for atom in ldW.atoms:
        #mol, typ, q, x, y, z=atom
        mol, typ, q, x, y, z, ix, iy, iz=ldW.atoms[atom]
        x, y, z=wrap(x, y, z, ldW)
        atompos[atom]=[mol, typ, q, x, y, z]
        count+=1
        if mol>nmol: nmol=mol
    #print(len(ldW.atoms), count)

    #0:0:for num of features and sets

    ATOMSECTION=[]    

    for atomid in range(1, count, 1):
        atom=atompos[atomid]
        mol, typ, q, x, y, z=atom 
        """
@<TRIPOS>ATOM
      1:  C        37.3010:   6.6870:  34.2710:C.3:    1: IG1:       -0.0483
      2:  C        37.1700:   5.9810:  35.6780:C.3:    1: IG1:        0.0813
      3:  C        36.6860:   7.1370:  36.6650:C.3:    1: IG1:       -0.0209
        """
        atyp=chkatyp(typ, types, ldW.masses)
        if atomid in car:
            atyp="C.2"
        elif atomid in oxy:
            atyp="O.co2"

        #atomtype in "C.3:for sp3, O.3:for sp3, O.2, H, O.co2, C.ar, N. pl3, N.ar, N.3, N.4, N.2"
        line="{cID:>7}{ele:>4}     {fx:10.4f}{fy:10.4f}{fz:10.4f} {atyp:<5} {mID:>3} {molStr:>4}     {charge:6.4f}\n".format(cID=atomid, ele=elementIDs[typ-1], mID=mol, molStr="IG"+str(mol), fx=x,fy=y,fz=z, charge=q, atyp=atyp)
        ATOMSECTION.append(line)


    bdict={}
    bcount=0
    #print("bondscount_vanilla: {}".format(len(ldW.bonds)))

    for bond in ldW.bonds:
        #typ, a1, a2=bond

        typ, a1, a2=ldW.bonds[bond]
        atype1=ldW.atoms[a1][1]
        atype2=ldW.atoms[a2][1]
        if ".ar" in chkatyp(atype1, types, ldW.masses) and ".ar" in chkatyp(atype2, types, ldW.masses):
            btyp="ar"
        elif "O.co2" in chkatyp(atype1, types, ldW.masses) or "O.co2" in chkatyp(atype2, types, ldW.masses):
            btyp="2"
        elif a1 in oxy or a2 in oxy:#special bond for C=O for HEMA-TRP connections
            btyp="2"
        else:
            btyp="1"
        #print(a1, a2)
        pos1=np.array(atompos[a1][3:6])
        pos2=np.array(atompos[a2][3:6])
        if (np.linalg.norm(pos1-pos2)>5.0):
            continue
        else:
            if a2>a1:
                bdict[bcount+1]=[a1, a2, btyp]
            else:
                bdict[bcount+1]=[a2, a1, btyp]
            bcount+=1

    header="""@<TRIPOS>MOLECULE
3DprintWorkflow
{} {} {} 0 0
SMALL
GASTEIGER

""".format(ldW.natoms, bcount, nmol)
    f.write(header)

    f.write("@<TRIPOS>ATOM\n")
    for line in ATOMSECTION:
        f.write(line)

    f.write("@<TRIPOS>BOND\n")
    for bidx in sorted(bdict):
        a1, a2, btyp=bdict[bidx]
        
        f.write("{bidx:>6} {a1:>5} {a2:>5} {btyp:>5}\n".format(bidx=bidx, a1=a1, a2=a2, btyp=btyp))

#print("bondscount: {}".format(bcount))
if len(ldW.atoms)>99999:
    print("Warning! Number of atoms exceeds limit {} > 99999".format(len(ldW.atoms)))
