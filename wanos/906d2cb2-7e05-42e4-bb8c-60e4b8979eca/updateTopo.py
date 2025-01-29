import LammpsData as LD
import numpy as np
from addmol import addmol

ldref=LD.LammpsData()
ldref.readFile("fRef_pos")

ldall=LD.LammpsData()
ldall.readFile("fCOORD")

masses=ldall.masses[:]

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

#if is critical atom and there is a bond, delete atom
#count number of atoms in the deleted list
#update mapping using a dictionary as update[atomid]=newid
#update bonds, angles, dihedrals, impropers

#atomid in ref=molid in allatom-1
delatoms=[]
BONDS=[]
ANGLES=[]

numMonomer, typMonomer, numXLinker, typXLinker, numLigand, typLigand, lSpacer=readcomponents()
numTypes={'DHPMA':[48, 47, 84, 115, 26], 'HEMA':[40, 39, 70, 91, 22], 'EGDMA':[18, 17, 28, 28, 9], 'TRP':[27, 28, 46, 69, 15], 'Q':[20, 19, 36, 45, 12], 'S':[5, 3, 3, 0, 0]}

def findmolecule(molID, typ, LDclass):
    found=-1
    for idx, atom in enumerate(LDclass.atoms):
        #print(idx, molID, atom[0])
        if atom[0]==molID:
            if typ==atom[1]:
                found=idx
                break
            elif typ==-1:
                print(atom)            
    return found

#dimerization 1(0)-2*(0) -> 4(1)-3(1) dimer btyp 1

#poly1 2*(0)-3(1) -> 4*(1)-5(2) extension btyp 2
#poly2 4*(1)-3(1) -> 5(2)-7*(2) extension btyp 2
#poly3 3(1)-7*(2) -> 7*(2)-6(3) branching of PMMA as of 10.1002/pola.23154? btyp 4
#poly4 1(0)-4*(1) -> 3(1)-7*(2) extension btyp 2
#poly5 2*(0)-4*(1) -> 3(1)-5(2) extension btyp 2
#poly6 1(0)-7*(2) -> 4*(1)-6(3) branching btyp 4
#poly7 2(0)-7*(2) -> 3(1)-6(3) branching btyp 4

#cross1 4*(1)-6(3) -> 7*(2)-8(4) crosslinking btyp 3
#cross2 2*(0)-6(3) -> 4*(1)-8(4) crosslinking btyp 3
#cross3 7*(2)-(5) -> 6(3)-6(3) GDMA involved? radical? crosslinking btyp 6

#end 3(1)-11(0) -> 5(2)-12(1) btyp 5 corrected

#make use of original typing: type 1 bond for dimerization, type 2 for polymerization, type 3 for crosslinking, type 4 for PMMA-END
"""
bond types=
1 HEMA0-HEMA0
2 HEMA0-HEMA1
3 EGDMA0-EGDMA0
4 EGDMA0-EGDMA1
5 HEMA0-EGDMA0
6 HEMA0-EGDMA1
7 HEMA1-EGDMA0
8 HEMA1-EGDMA1
9 HEMA2-TRP
"""

def candDHPMA(a, ldall, aoffset):
    Aidx1=findmolecule(a, 26+aoffset, ldall) #R H
    Aidx2=findmolecule(a, 31+aoffset, ldall) #L H
    Aidx3=findmolecule(a, 2+aoffset, ldall) #R C
    Aidx4=findmolecule(a, 5+aoffset, ldall) #L C
    AidxCOH1C=findmolecule(a, 13+aoffset, ldall) #CH2OH1 C
    AidxCOH1HR=findmolecule(a, 41+aoffset, ldall) #COH1 OH H
    AidxCOH1H1=findmolecule(a, 40+aoffset, ldall) #COH1 COH1 H
    AidxCOH1H2=findmolecule(a, 39+aoffset, ldall) #COH1 COH2 H
    AidxCOH1O=findmolecule(a, 14+aoffset, ldall) #COH1 OH O
        
    AidxCOH2C=findmolecule(a, 21+aoffset, ldall) #CH2OH2 C
    AidxCOH2HR=findmolecule(a, 48+aoffset, ldall) #COH2 OH H
    AidxCOH2H1=findmolecule(a, 47+aoffset, ldall) #COH2 COH1 H
    AidxCOH2H2=findmolecule(a, 46+aoffset, ldall) #COH2 COH2 H
    AidxCOH2O=findmolecule(a, 22+aoffset, ldall) #COH2 OH O

    AidxCO1C=findmolecule(a, 7+aoffset, ldall) #C=O1 C
    AidxCO1O=findmolecule(a, 8+aoffset, ldall) #C=O1 O    
    AidxCO2C=findmolecule(a, 15+aoffset, ldall) #C=O2 C
    AidxCO2O=findmolecule(a, 16+aoffset, ldall) #C=O2 O 

    #remove H (idx1) and connect C (idx2)
    return Aidx1, Aidx2, Aidx3, Aidx4, AidxCOH1C, AidxCOH1HR, AidxCOH1H1, AidxCOH1H2, AidxCOH1O, AidxCOH2C, AidxCOH2HR, AidxCOH2H1, AidxCOH2H2, AidxCOH2O, AidxCO1C, AidxCO1O, AidxCO2C, AidxCO2O

def candHEMA(a, ldall, aoffset):
    Aidx1=findmolecule(a, 22+aoffset, ldall) #R H
    Aidx2=findmolecule(a, 27+aoffset, ldall) #L H
    Aidx3=findmolecule(a, 2+aoffset, ldall) #R C
    Aidx4=findmolecule(a, 5+aoffset, ldall) #L C
    AidxCOH1C=findmolecule(a, 17+aoffset, ldall) #COH C
    AidxCOH1HR=findmolecule(a, 40+aoffset, ldall) #COH OH H
    AidxCOH1H1=findmolecule(a, 38+aoffset, ldall) #COH COH1 H
    AidxCOH1H2=findmolecule(a, 39+aoffset, ldall) #COH COH2 H
    AidxCOH1O=findmolecule(a, 18+aoffset, ldall) #COH OH O
    
    AidxCOH2C=findmolecule(a, 11+aoffset, ldall) #COH C
    AidxCOH2HR=findmolecule(a, 35+aoffset, ldall) #COH OH H
    AidxCOH2H1=findmolecule(a, 33+aoffset, ldall) #COH COH1 H
    AidxCOH2H2=findmolecule(a, 34+aoffset, ldall) #COH COH2 H
    AidxCOH2O=findmolecule(a, 12+aoffset, ldall) #COH OH O
    
    AidxCO1C=findmolecule(a, 7+aoffset, ldall) #C=O1 C
    AidxCO1O=findmolecule(a, 8+aoffset, ldall) #C=O1 O    
    AidxCO2C=findmolecule(a, 13+aoffset, ldall) #C=O2 C
    AidxCO2O=findmolecule(a, 14+aoffset, ldall) #C=O2 O 

    #remove H (idx1) and connect C (idx2)
    return Aidx1, Aidx2, Aidx3, Aidx4, AidxCOH1C, AidxCOH1HR, AidxCOH1H1, AidxCOH1H2, AidxCOH1O, AidxCOH2C, AidxCOH2HR, AidxCOH2H1, AidxCOH2H2, AidxCOH2O, AidxCO1C, AidxCO1O, AidxCO2C, AidxCO2O

def candEGDMA(a, ldall, aoffset):
    #atom types are offset due to different types in monomer, otherwise they cannot be found properly
    Aidx1=findmolecule(a, 2+aoffset, ldall) #RR C
    Aidx2=findmolecule(a, 1+aoffset, ldall) #RL C
    Aidx3=findmolecule(a, 17+aoffset, ldall) #C H
    Aidx4=findmolecule(a, 8+aoffset, ldall) #C C
    #RR, RL, H, C
    return Aidx1, Aidx2, Aidx3, Aidx4
    
def candTRP(a, ldall, aoffset):
    Aidx1=findmolecule(a, 18+aoffset, ldall) #R H
    Aidx2=findmolecule(a, 1+aoffset, ldall) #R N
    #COOH end
    ACOOH1=findmolecule(a, 14+aoffset, ldall)
    ACOOH2=findmolecule(a, 15+aoffset, ldall)
    ACOOH3=findmolecule(a, 16+aoffset, ldall)
    ACOOH4=findmolecule(a, 27+aoffset, ldall)
    #spacer
    ALeftL=findmolecule(a, ldall.natomtypes-1, ldall)
    ARightL=findmolecule(a, ldall.natomtypes, ldall)
    return Aidx1, Aidx2, ACOOH1, ACOOH2, ACOOH3, ACOOH4, ALeftL, ARightL

def candQ(a, ldall, aoffset):
    Aidx1=findmolecule(a, 7+aoffset, ldall) #R H

    Aidx2=findmolecule(a, 3+aoffset, ldall) #R N
    #spacer
    ALeftL=findmolecule(a, ldall.natomtypes-1, ldall)
    ARightL=findmolecule(a, ldall.natomtypes, ldall)
    
    Aidxtail1=findmolecule(a, 1+aoffset, ldall) #R C
    Aidxtail2=findmolecule(a, 2+aoffset, ldall) #R C
    Aidxtail3=findmolecule(a, 7+aoffset, ldall) #R H
    Aidxtail4=findmolecule(a, 8+aoffset, ldall) #R H
    Aidxtail5=findmolecule(a, 9+aoffset, ldall) #R H
    Aidxtail6=findmolecule(a, 10+aoffset, ldall) #R H    
    Aidxtail7=findmolecule(a, 11+aoffset, ldall) #R H
    return Aidx1, Aidx2, Aidxtail1, Aidxtail2, Aidxtail3, Aidxtail4, Aidxtail5, Aidxtail6, Aidxtail7, ALeftL, ARightL
    
def candS(a, ldall, aoffset):
    Aidx1=findmolecule(a, 1+aoffset, ldall)
    Aidx2=findmolecule(a, 2+aoffset, ldall) #O in S=O
    Aidx3=findmolecule(a, 3+aoffset, ldall) #O in S=O
    Aidx4=findmolecule(a, 4+aoffset, ldall) #assign as single S-O bond 
    AidxDummy=findmolecule(a, 5+aoffset, ldall) #assign as single S-O bond 
    ALeftL=findmolecule(a, ldall.natomtypes-1, ldall)
    ARightL=findmolecule(a, ldall.natomtypes, ldall)
    return AidxDummy, Aidx1, Aidx2, Aidx3, Aidx4, ALeftL, ARightL

#remove the cap on HEMA, and use O on C-OH as C=O
#need to record the indices of the C=O bond during the conversion to .mol2 in the bond definition
specialMarkers=[] 

bondOccu=[]#recording bonds used in EGDMA: since no H caps indicate whether or not the C=C bond is already opened, a dedicated list is used for tracking
for bid, bond in enumerate(ldref.bonds):
    typ, a, b=bond
    if ldref.atoms[a-1][1] in [1, 2, 3, 4, 5, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]:
        typA="Monomer"
    elif ldref.atoms[a-1][1] in [6, 7, 8, 9, 10, 12]:
        typA="Crosslinker"
    elif ldref.atoms[a-1][1] in [14, 15]:
        typA="Ligand"
    if ldref.atoms[b-1][1] in [1, 2, 3, 4, 5, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]:
        typB="Monomer"
    elif ldref.atoms[b-1][1] in [6, 7, 8, 9, 10, 12]:
        typB="Crosslinker"
    elif ldref.atoms[b-1][1] in [14, 15]:
        typB="Ligand"
    if typA=="Monomer":
        if typMonomer=="HEMA":
            candA=candHEMA(a, ldall, aoffset=0)
        elif typMonomer=="DHPMA":
            candA=candDHPMA(a, ldall, aoffset=0)
    elif typA=="Crosslinker":
        candA=candEGDMA(a, ldall, aoffset=numTypes[typMonomer][0])
    elif typA=="Ligand":
        if typLigand=="TRP":
            candA=candTRP(a, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])
        elif typLigand=="Q":
            candA=candQ(a, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])
        elif typLigand=="S":
            candA=candS(a, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])

    if typB=="Monomer":
        if typMonomer=="HEMA":
            candB=candHEMA(b, ldall, aoffset=0)
        elif typMonomer=="DHPMA":
            candB=candDHPMA(b, ldall, aoffset=0)
    elif typB=="Crosslinker":
        candB=candEGDMA(b, ldall, aoffset=numTypes[typMonomer][0])
    elif typB=="Ligand":
        if typLigand=="TRP":
            candB=candTRP(b, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])
        elif typLigand=="Q":
            candB=candQ(b, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])
        elif typLigand=="S":
            candB=candS(b, ldall, aoffset=numTypes[typMonomer][0]+numTypes[typXLinker][0])
            
    #if a in [52]:
    #    print(typ, a, b)
    #    print(numTypes[typMonomer][0]+numTypes[typXLinker][0])
    #    print(typA, typB)
    #    print(candA, candB)
    #    print(typLigand, typMonomer, typXLinker)
    if typ in [1, 2]:
        if typA=="Monomer" and typB=="Monomer":#identical protocols for HEMA and DHPMA
            if not candA[2] in bondOccu:#ldref.atoms[a-1][3]<ldref.atoms[b-1][3]: #A on the left, connect right end to B
                #DHPMA2-DHPMA2 dimers
                delatoms.append(candA[0]+1)
                idx1=candA[2]
                #torsional definitions are not updated
            elif not candA[3] in bondOccu:
                delatoms.append(candA[1]+1)
                idx1=candA[3]
            else:
                raise("full")

            if not candB[3] in bondOccu:
                #DHPMA2-DHPMA2 dimers
                delatoms.append(candB[1]+1)
                idx2=candB[3]
            elif not candB[2] in bondOccu:
                delatoms.append(candB[0]+1)
                idx2=candB[2]
            else:
                raise("full")
            BONDS.append([ldall.nbondtypes+1, idx1+1, idx2+1])
            bondOccu.append(idx1)
            bondOccu.append(idx2)

            #torsional definitions are not updated

        elif typA=="Crosslinker" and typB=="Crosslinker": #default to bond type 1 during create_box
            #connect crosslinking site
            delatoms.append(candA[2]+1)
            delatoms.append(candB[2]+1)
            BONDS.append([ldall.nbondtypes+2, candA[3]+1, candB[3]+1])
                

    elif typ in [3, 4]: #EGDMA-EGDMA
        
        #open C=C
        if not candA[1] in bondOccu:
            idx1=candA[1]
        elif not candA[0] in bondOccu:
            idx1=candA[0]
        else:
            raise("no free C=C bonds")
        if not candB[0] in bondOccu:
            idx2=candB[0]
        elif not candB[1] in bondOccu:
            idx2=candB[1]
        else:
            raise("no free C=C bonds")

        BONDS.append([ldall.nbondtypes+3, idx1+1, idx2+1])
        bondOccu.append(idx1)
        bondOccu.append(idx2)
                
    elif typ in [5, 6, 7, 8]:
        #DHPMA-EDGMA
        if typA=="Monomer" and typMonomer=='DHPMA':
            #debug
            #print(candA)
            #print(candB)
            if not candA[2] in bondOccu:#ldref.atoms[a-1][3]<ldref.atoms[b-1][3]: #HEMA on left
                delatoms.append(candA[0]+1)
                idx1=candA[2]
            elif not candA[3] in bondOccu:
                delatoms.append(candA[1]+1)
                idx1=candA[3]
            else:
                raise("no free C=C on DHPMA")
            if not candB[0] in bondOccu:
                idx2=candB[0]
            elif not candB[1] in bondOccu:
                idx2=candB[1]
            else:
                raise("no free C=C on DHPMA")
            BONDS.append([ldall.nbondtypes+4, idx1+1, idx2+1])
            bondOccu.append(idx1)
            bondOccu.append(idx2)
        elif typA=="Monomer" and typMonomer=='HEMA':
            if not candA[2] in bondOccu:#ldref.atoms[a-1][3]<ldref.atoms[b-1][3]: #HEMA on left
                delatoms.append(candA[0]+1)
                idx1=candA[2]
            elif not candA[3] in bondOccu:
                delatoms.append(candA[1]+1)
                idx1=candA[3]
            else:
                raise("full")
            if not candB[0] in bondOccu:
                idx2=candB[0]
            elif not candB[1] in bondOccu:
                idx2=candB[1]
            else:
                raise("full")
            BONDS.append([ldall.nbondtypes+4, idx1+1, idx2+1])
            bondOccu.append(idx1)
            bondOccu.append(idx2)

        else: # typ "Crosslinker"
            if not candA[0] in bondOccu:
                idx1=candA[0]
            elif not candA[1] in bondOccu:
                idx1=candA[1]
            else:
                raise("no free C=C on crosslinker")
            if not candB[2] in bondOccu:
                delatoms.append(candB[0]+1)
                idx2=candB[2]
            elif not candB[3] in bondOccu:
                delatoms.append(candB[1]+1)
                idx2=candB[3]
            else:
                raise("no free C=C on crosslinker")
            BONDS.append([ldall.nbondtypes+4, idx1+1, idx2+1])
            bondOccu.append(idx1)
            bondOccu.append(idx2)
            #torsional definitions are not updated

    elif typ in [9]: 
        #Monomer-Ligand, remove 2Hs on C and 1H on O in CH2OH, always remove H cap (dummy cap) as 1st atom on Ligand
        if typA=="Monomer" and typMonomer=='DHPMA':
            if typLigand=="TRP": #one CH2OH branch turns into C=O and C-N bond is formed
                delatoms.append(candA[5]+1)
                delatoms.append(candA[6]+1)
                delatoms.append(candA[7]+1)
                delatoms.append(candB[0]+1)
                ldall.atoms[candA[8]][2]=-0.450 #O for OH turns into O.co2
                specialMarkers.append([candA[8], 'O.co2'])
                specialMarkers.append([candA[4], 'C.2'])

            elif typLigand=="S":
                delatoms.append(candA[5]+1)
                #delatoms.append(candA[6]+1)
                #delatoms.append(candA[7]+1)
                delatoms.append(candA[8]+1)
                delatoms.append(candB[0]+1)
            elif typLigand=="Q":
                delatoms.append(candA[5]+1)
                #delatoms.append(candA[6]+1)
                #delatoms.append(candA[7]+1)
                delatoms.append(candA[8]+1)
                delatoms.append(candB[0]+1)
                delatoms.append(candB[2]+1)
                delatoms.append(candB[3]+1)
                delatoms.append(candB[4]+1)
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)
            
            if typLigand=="S" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[4]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])

                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="S":
                ANGLES.append([ldall.nangletypes+1, candA[4]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candA[4]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candA[4]+1, candB[1]+1, candB[4]+1])
            elif lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[4]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            else:
                BONDS.append([ldall.nbondtypes+5, candA[4]+1, candB[1]+1])

        elif typA=="Monomer" and typMonomer=='HEMA':
            if typLigand=="TRP":
                delatoms.append(candA[5]+1)
                delatoms.append(candA[6]+1)
                delatoms.append(candA[7]+1)
                delatoms.append(candB[0]+1)
                ldall.atoms[candA[8]][2]=-0.450
                specialMarkers.append([candA[8], 'O.co2'])
                specialMarkers.append([candA[4], 'C.2'])
            elif typLigand=="S":
                delatoms.append(candA[5]+1)
                #delatoms.append(candA[6]+1)
                #delatoms.append(candA[7]+1)
                delatoms.append(candB[0]+1)
                
            elif typLigand=="Q":
                delatoms.append(candA[5]+1)
                #delatoms.append(candA[6]+1)
                #delatoms.append(candA[7]+1)
                delatoms.append(candB[0]+1)
                delatoms.append(candB[2]+1)
                delatoms.append(candB[3]+1)
                delatoms.append(candB[4]+1)
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)

            if typLigand=="S" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[8]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])

                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="S":
                ANGLES.append([ldall.nangletypes+1, candA[8]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candA[8]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candA[8]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="Q" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[8]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            elif typLigand=="Q":
                BONDS.append([ldall.nbondtypes+5, candA[8]+1, candB[1]+1])
            elif typLigand=="TRP" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[4]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            elif typLigand=="TRP":
                BONDS.append([ldall.nbondtypes+5, candA[4]+1, candB[1]+1])

        elif typA=="Ligand" and typLigand=="TRP":
            if typMonomer=="HEMA":
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candA[0]+1)
            elif typMonomer=="DHPMA":
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candA[0]+1)

            if lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[4]+1])
            else:
                BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[4]+1])

            ldall.atoms[candB[8]][2]=-0.450
            specialMarkers.append([candB[8], 'O.co2'])
            specialMarkers.append([candB[4], 'C.2'])
            #torsional definitions are not updated
        elif typA=="Ligand" and typLigand=="Q":
            if typMonomer=="DHPMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[5]+1)
                #delatoms.append(candB[6]+1)
                #delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)
                
                delatoms.append(candA[2]+1)
                delatoms.append(candA[3]+1)
                delatoms.append(candA[4]+1)
                delatoms.append(candA[5]+1)
                delatoms.append(candA[6]+1)
                delatoms.append(candA[7]+1)
                delatoms.append(candA[8]+1)
                
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[4]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[4]+1])

            elif typMonomer=="HEMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[5]+1)
                #delatoms.append(candB[6]+1)
                #delatoms.append(candB[7]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[8]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[8]+1])
                
        elif typA=="Ligand" and typLigand=="S":
            if typMonomer=="DHPMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[5]+1)
                #delatoms.append(candB[6]+1)
                #delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[4]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[4]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[4]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[4]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[4]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[4]+1, candA[1]+1, candA[4]+1])

            elif typMonomer=="HEMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[5]+1)
                #delatoms.append(candB[6]+1)
                #delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[8]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[4]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[8]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[8]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[8]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[8]+1, candA[1]+1, candA[4]+1])

    elif typ in [10]: #type 10 bonds are for connecting the 2nd Ligand on a dimer, i.e. the "left" branch
        #DHPMA-TRP, remove 2Hs on C and 1H on O in CH2OH 
        if typA=="Monomer" and typMonomer=='DHPMA':
            if typLigand=="TRP": #one CH2OH branch turns into C=O and C-N bond is formed
                delatoms.append(candA[10]+1)
                delatoms.append(candA[11]+1)
                delatoms.append(candA[12]+1)
                delatoms.append(candB[0]+1)
                ldall.atoms[candA[13]][2]=-0.450 #O for OH turns into O.co2
                specialMarkers.append([candA[13], 'O.co2'])
                specialMarkers.append([candA[9], 'C.2'])

            elif typLigand=="S":
                delatoms.append(candA[10]+1)
                #delatoms.append(candA[11]+1)
                #delatoms.append(candA[12]+1)
                delatoms.append(candA[13]+1)
                delatoms.append(candB[0]+1)
                
            elif typLigand=="Q":
                delatoms.append(candA[10]+1)
                #delatoms.append(candA[11]+1)
                #delatoms.append(candA[12]+1)
                delatoms.append(candA[13]+1)
                delatoms.append(candB[0]+1)
                delatoms.append(candB[2]+1)
                delatoms.append(candB[3]+1)
                delatoms.append(candB[4]+1)
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)
            
            if typLigand=="S" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[9]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])

                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="S":
                ANGLES.append([ldall.nangletypes+1, candA[9]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candA[9]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candA[9]+1, candB[1]+1, candB[4]+1])
            elif lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[9]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            else:
                BONDS.append([ldall.nbondtypes+5, candA[9]+1, candB[1]+1])

        elif typA=="Monomer" and typMonomer=='HEMA':
            if typLigand=="TRP":
                delatoms.append(candA[10]+1)
                delatoms.append(candA[11]+1)
                delatoms.append(candA[12]+1)
                delatoms.append(candB[0]+1)
                ldall.atoms[candA[13]][2]=-0.450
                specialMarkers.append([candA[13], 'O.co2'])
                specialMarkers.append([candA[9], 'C.2'])
            elif typLigand=="S":
                delatoms.append(candA[10]+1)
                #delatoms.append(candA[11]+1)
                #delatoms.append(candA[12]+1)
                delatoms.append(candB[0]+1)
            elif typLigand=="Q":
                delatoms.append(candA[10]+1)
                #delatoms.append(candA[11]+1)
                #delatoms.append(candA[12]+1)

                delatoms.append(candB[0]+1)
                delatoms.append(candB[2]+1)
                delatoms.append(candB[3]+1)
                delatoms.append(candB[4]+1)
                delatoms.append(candB[5]+1)
                delatoms.append(candB[6]+1)
                delatoms.append(candB[7]+1)
                delatoms.append(candB[8]+1)

            if typLigand=="S" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[13]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])

                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candB[-2]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="S":
                ANGLES.append([ldall.nangletypes+1, candA[13]+1, candB[1]+1, candB[2]+1])
                ANGLES.append([ldall.nangletypes+1, candA[13]+1, candB[1]+1, candB[3]+1])
                ANGLES.append([ldall.nangletypes+1, candA[13]+1, candB[1]+1, candB[4]+1])
            elif typLigand=="Q" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[13]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            elif typLigand=="Q":
                BONDS.append([ldall.nbondtypes+5, candA[13]+1, candB[1]+1])
            elif typLigand=="TRP" and lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[9]+1, candB[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candB[-1]+1, candB[1]+1])
            elif typLigand=="TRP":
                BONDS.append([ldall.nbondtypes+5, candA[9]+1, candB[1]+1])

        elif typA=="Ligand" and typLigand=="TRP":
            if typMonomer=="HEMA":
                delatoms.append(candB[10]+1)
                delatoms.append(candB[11]+1)
                delatoms.append(candB[12]+1)
                delatoms.append(candA[0]+1)
            elif typMonomer=="DHPMA":
                delatoms.append(candB[10]+1)
                delatoms.append(candB[11]+1)
                delatoms.append(candB[12]+1)
                delatoms.append(candA[0]+1)

            if lSpacer:
                BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[9]+1])
            else:
                BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[9]+1])

            ldall.atoms[candB[13]][2]=-0.450
            specialMarkers.append([candB[13], 'O.co2'])
            specialMarkers.append([candB[9], 'C.2'])
            #torsional definitions are not updated
        elif typA=="Ligand" and typLigand=="Q":
            if typMonomer=="DHPMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candA[2]+1)
                delatoms.append(candA[3]+1)
                delatoms.append(candA[4]+1)
                delatoms.append(candA[5]+1)
                delatoms.append(candA[6]+1)
                delatoms.append(candA[7]+1)
                delatoms.append(candA[8]+1)
                
                delatoms.append(candB[10]+1)
                #delatoms.append(candB[11]+1)
                #delatoms.append(candB[12]+1)
                delatoms.append(candB[13]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[9]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[9]+1])

            elif typMonomer=="HEMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[10]+1)
                #delatoms.append(candB[11]+1)
                #delatoms.append(candB[12]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[13]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[13]+1])
                
        elif typA=="Ligand" and typLigand=="S":
            if typMonomer=="DHPMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[10]+1)
                #delatoms.append(candB[11]+1)
                #delatoms.append(candB[12]+1)
                delatoms.append(candB[13]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[9]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[4]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[9]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[9]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[9]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[9]+1, candA[1]+1, candA[4]+1])

            elif typMonomer=="HEMA":
                delatoms.append(candA[0]+1)
                delatoms.append(candB[10]+1)
                #delatoms.append(candB[11]+1)
                #delatoms.append(candB[12]+1)
                if lSpacer:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candA[-2]+1])
                    BONDS.append([ldall.nbondtypes+5, candA[-1]+1, candB[13]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candA[-2]+1, candA[1]+1, candA[4]+1])
                else:
                    BONDS.append([ldall.nbondtypes+5, candA[1]+1, candB[8]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[13]+1, candA[1]+1, candA[2]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[13]+1, candA[1]+1, candA[3]+1])
                    ANGLES.append([ldall.nangletypes+1, candB[13]+1, candA[1]+1, candA[4]+1])

AATOMS=[]
delatoms=list(set(delatoms)) #unique atom ids in the remove list
#print(delatoms)
#print(BONDS)

typs=[0,0,0,0,0]
for idx,atom in enumerate(ldall.atoms):
    if idx+1 in delatoms:
        continue
    else:
        AATOMS.append(atom[:])
        if typs[0]<atom[1]:
            typs[0]=atom[1]

ABONDS=[]
for bond in ldall.bonds:
    if bond[1] in delatoms or bond[2] in delatoms:
        continue
    else:
        typ, a, b=bond
        ABONDS.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms)])
        if typs[1]<typ:
            typs[1]=typ

AANGLES=[]
for bond in ldall.angles:
    if bond[1] in delatoms or bond[2] in delatoms or bond[3] in delatoms:
        continue
    else:
        typ, a, b, c=bond
        AANGLES.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms), c - sum(i < c for i in delatoms)])
        if typs[2]<typ:
            typs[2]=typ

ADIHEDRALS=[]
for bond in ldall.dihedrals:
    if bond[1] in delatoms or bond[2] in delatoms or bond[3] in delatoms or bond[4] in delatoms:
        continue
    else:
        typ, a, b, c, d=bond
        ADIHEDRALS.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms), c - sum(i < c for i in delatoms), d - sum(i < d for i in delatoms)])
        if typs[3]<typ:
            typs[3]=typ

AIMPROPERS=[]
for bond in ldall.impropers:
    if bond[1] in delatoms or bond[2] in delatoms or bond[3] in delatoms or bond[4] in delatoms:
        continue
    else:
        typ, a, b, c, d=bond
        AIMPROPERS.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms), c - sum(i < c for i in delatoms), d - sum(i < d for i in delatoms)])
        if typs[4]<typ:
            typs[4]=typ

for bond in BONDS:
    typ, a, b=bond
    ABONDS.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms)])
    if typs[1]<typ:
        typs[1]=typ
        
for angle in ANGLES:
    typ, a, b, c=angle
    AANGLES.append([typ, a - sum(i < a for i in delatoms), b - sum(i < b for i in delatoms), c - sum(i < c for i in delatoms)])
    if typs[2]<typ:
        typs[2]=typ

ldall.atoms=AATOMS
ldall.bonds=ABONDS
ldall.angles=AANGLES
ldall.dihedrals=ADIHEDRALS
ldall.impropers=AIMPROPERS

ldall.natoms=len(ldall.atoms)
ldall.nbonds=len(ldall.bonds)
ldall.nangles=len(ldall.angles)
ldall.ndihedrals=len(ldall.dihedrals)
ldall.nimpropers=len(ldall.impropers)

if lSpacer:
    ldall.natomtypes=numTypes[typMonomer][0]+numTypes[typXLinker][0]+numTypes[typLigand][0]+5
    ldall.nbondtypes=numTypes[typMonomer][1]+numTypes[typXLinker][1]+numTypes[typLigand][1]+2+5
    ldall.nangletypes=numTypes[typMonomer][2]+numTypes[typXLinker][2]+numTypes[typLigand][2]+1+3
else:
    ldall.natomtypes=numTypes[typMonomer][0]+numTypes[typXLinker][0]+numTypes[typLigand][0]
    ldall.nbondtypes=numTypes[typMonomer][1]+numTypes[typXLinker][1]+numTypes[typLigand][1]+5
    ldall.nangletypes=numTypes[typMonomer][2]+numTypes[typXLinker][2]+numTypes[typLigand][2]+1
    
ldall.ndihedraltypes=numTypes[typMonomer][3]+numTypes[typXLinker][3]+numTypes[typLigand][3]
ldall.nimpropertypes=numTypes[typMonomer][4]+numTypes[typXLinker][4]+numTypes[typLigand][4]

ldall.masses=masses

#print(ldall.natoms, ldall.nbonds, ldall.nangles, ldall.ndihedrals, ldall.nimpropers)
#ldall.writeFile("fTOPO_NOCOV")
ldall.writeFile("fTOPO")

with open("COrecorder.dat", 'w') as f:
    for idx in range(len(specialMarkers)):
        aid, typ=specialMarkers[idx]
        atomreID=aid-sum(i<aid for i in delatoms)+1
        f.write("{} {}\n".format(atomreID, typ))
