import LammpsData as LD
import numpy as np

import sys

ldref=LD.LammpsData()
ldref.readFile("fRef_pos")

ldref1=LD.LammpsData()

nspacer=int(sys.argv[1])

ldref1.readFile("Monomer.lmp")

ldref2=LD.LammpsData()
ldref2.readFile("Crosslinkerhalf.lmp")

ldref3=LD.LammpsData()
ldref3.readFile("Ligand.lmp")

ldout=LD.LammpsData()
boxx=ldref.box[0][1]-ldref.box[0][0]
boxy=ldref.box[1][1]-ldref.box[1][0]
boxz=ldref.box[2][1]-ldref.box[2][0]
ldout.box=[[0.0, boxx*11], [0.0, boxy*11], [0.0, boxz*11]]
curmol=0
curnatoms, curnbonds, curnangles, curndihedrals, curnimpropers=(0, 0, 0, 0, 0)
curnatomtypes, curnbondtypes, curnangletypes, curndihedraltypes, curnimpropertypes=(0, 0, 0, 0, 0)

for idx, refatom in enumerate(ldref.atoms):

    if refatom[1] in [3, 4, 5, 11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]: #ignore individual (unbound) monomers (type 1 and 2)
        for patom in ldref1.atoms:
            mol=idx+1
            typ=patom[1]+curnatomtypes
            q=patom[2]
            x=patom[3]+refatom[3]*11
            y=patom[4]+refatom[4]*11
            z=patom[5]+refatom[5]*11
            ldout.atoms.append([mol, typ, q, x, y, z])

        for patom in ldref1.bonds:
            typ=patom[0]+curnbondtypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            ldout.bonds.append([typ, a, b])
        for patom in ldref1.angles:
            typ=patom[0]+curnangletypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            ldout.angles.append([typ, a, b, c])
        for patom in ldref1.dihedrals:
            typ=patom[0]+curndihedraltypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.dihedrals.append([typ, a, b, c, d])
        for patom in ldref1.impropers:
            typ=patom[0]+curnimpropertypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.impropers.append([typ, a, b, c, d])
        curnatoms=len(ldout.atoms)

curnatomtypes=np.max([atom[1] for atom in ldout.atoms])
curnbondtypes=np.max([atom[0] for atom in ldout.bonds])
curnangletypes=np.max([atom[0] for atom in ldout.angles])
curndihedraltypes=np.max([atom[0] for atom in ldout.dihedrals])
curnimpropertypes=np.max([atom[0] for atom in ldout.impropers])

for idx, refatom in enumerate(ldref.atoms):
    if refatom[1] in [6, 7, 8, 9, 10, 12]: #should be adding half of the EGDMA to the all-atom representation
    #treat dimer as one or add half
        for patom in ldref2.atoms:
            mol=idx+1
            typ=patom[1]+curnatomtypes
            q=patom[2]
            x=patom[3]+refatom[3]*11
            y=patom[4]+refatom[4]*11
            z=patom[5]+refatom[5]*11
            ldout.atoms.append([mol, typ, q, x, y, z])
        for patom in ldref2.bonds:
            typ=patom[0]+curnbondtypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            ldout.bonds.append([typ, a, b])
        for patom in ldref2.angles:
            typ=patom[0]+curnangletypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            ldout.angles.append([typ, a, b, c])
        for patom in ldref2.dihedrals:
            typ=patom[0]+curndihedraltypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.dihedrals.append([typ, a, b, c, d])
        for patom in ldref2.impropers:
            typ=patom[0]+curnimpropertypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.impropers.append([typ, a, b, c, d])
        curnatoms=len(ldout.atoms)

curnatomtypes=np.max([atom[1] for atom in ldout.atoms])
curnbondtypes=np.max([atom[0] for atom in ldout.bonds])
curnangletypes=np.max([atom[0] for atom in ldout.angles])
curndihedraltypes=np.max([atom[0] for atom in ldout.dihedrals])
curnimpropertypes=np.max([atom[0] for atom in ldout.impropers])

templateATOMS=[[0.0, 0.0, 0.0], [0.0, 1.09, 0.0], [0.0, -1.09, 0.0]]
templateBONDS=[[curnbondtypes+ldref3.nbondtypes, 1, 2], [curnbondtypes+ldref3.nbondtypes, 1, 3]]
templateANGLES=[[curnangletypes+ldref3.nangletypes,2,1,3]]

for idx, refatom in enumerate(ldref.atoms):
    if refatom[1] in [15]: #ignore individual (unbound) ligands (type 14)
        curnatoms=len(ldout.atoms)
        for patom in ldref3.atoms:
            mol=idx+1
            typ=patom[1]+curnatomtypes
            q=patom[2]
            x=patom[3]+refatom[3]*11
            y=patom[4]+refatom[4]*11
            z=patom[5]+refatom[5]*11
            ldout.atoms.append([mol, typ, q, x, y, z])
        for patom in ldref3.bonds:
            typ=patom[0]+curnbondtypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            ldout.bonds.append([typ, a, b])
        for patom in ldref3.angles:
            typ=patom[0]+curnangletypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            ldout.angles.append([typ, a, b, c])
        for patom in ldref3.dihedrals:
            typ=patom[0]+curndihedraltypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.dihedrals.append([typ, a, b, c, d])
        for patom in ldref3.impropers:
            typ=patom[0]+curnimpropertypes
            a=patom[1]+curnatoms
            b=patom[2]+curnatoms
            c=patom[3]+curnatoms
            d=patom[4]+curnatoms
            ldout.impropers.append([typ, a, b, c, d])

        curnatoms=len(ldout.atoms)
            
        spacerATOMS=[]#np.empty((0,3), dtype=np.float64)
        spacerBONDS=[]#np.empty((0,3), dtype=int)
        spacerANGLES=[]#np.empty((0,4), dtype=int)
        #spacerDIHEDRALS=[]
        
        x0, y0, z0=x, y, z
            
        first=True
        second=False
        curnatoms=len(ldout.atoms)
        for n in range(nspacer):

            offsetcoord=np.add([n*2.91, -4.0, -4.0], [x0, y0, z0]) #offset for the -CH2- short chain
            curATOMS=np.add(offsetcoord, templateATOMS)
            curBONDS=np.add([[1, n*3+curnatoms, n*3+curnatoms], [1, n*3+curnatoms, n*3+curnatoms]], templateBONDS)
            curANGLES=np.add([[1, n*3+curnatoms, n*3+curnatoms, n*3+curnatoms]], templateANGLES) #H-C-H angle
            for idx, patom in enumerate(curATOMS):
                x, y, z=patom
                if idx % 3 == 0: #carbon
                    
                    if n==0: #C in left-most CH2 group
                        if first:
                            ldout.atoms.append([mol, curnatomtypes+ldref3.natomtypes+4, -0.1804, x, y, z])
                            first=False
                            second=True
                        elif second:
                            ldout.atoms.append([mol, curnatomtypes+ldref3.natomtypes+3, -0.1804, x, y, z])
                                            
                    elif n==nspacer-1: #right end
                        ldout.atoms.append([mol, curnatomtypes+ldref3.natomtypes+5, -0.1804, x, y, z])
                    else:
                        ldout.atoms.append([mol, curnatomtypes+ldref3.natomtypes+1, -0.1804, x, y, z])                
                else: #hydrogen
                    ldout.atoms.append([mol, curnatomtypes+ldref3.natomtypes+2, 0.0902, x, y, z])

            for pbond in curBONDS:
                typ, a, b=pbond
                ldout.bonds.append([typ, a, b])
            for pangle in curANGLES:
                typ, a, b, c=pangle
                ldout.angles.append([typ, a, b, c])
            if n>0:
                ldout.bonds.append([curnbondtypes+ldref3.nbondtypes+2, (n-1)*3+1+curnatoms, n*3+1+curnatoms])#C-C bonds
                ldout.angles.append([curnangletypes+ldref3.nangletypes+2, (n-1)*3+2+curnatoms, (n-1)*3+1+curnatoms, n*3+1+curnatoms]) #C-C-H bond #1
                ldout.angles.append([curnangletypes+ldref3.nangletypes+2, (n-1)*3+3+curnatoms, (n-1)*3+1+curnatoms, n*3+1+curnatoms]) #C-C-H bond #2
                ldout.angles.append([curnangletypes+ldref3.nangletypes+2, (n-1)*3+1+curnatoms, n*3+1+curnatoms, n*3+2+curnatoms]) #C-C-H bond #3
                ldout.angles.append([curnangletypes+ldref3.nangletypes+2, (n-1)*3+1+curnatoms, n*3+1+curnatoms, n*3+3+curnatoms]) #C-C-H bond #4
                
            if n>1:
                ldout.angles.append([curnangletypes+ldref3.nangletypes+3, (n-2)*3+1+curnatoms, (n-1)*3+1+curnatoms, n*3+1+curnatoms]) #C-C-C angle                     

ldout.natomtypes=np.max([atom[1] for atom in ldout.atoms])
ldout.nbondtypes=np.max([atom[0] for atom in ldout.bonds])
ldout.nangletypes=np.max([atom[0] for atom in ldout.angles])
ldout.ndihedraltypes=np.max([atom[0] for atom in ldout.dihedrals])
ldout.nimpropertypes=np.max([atom[0] for atom in ldout.impropers])

if nspacer:
    ldout.masses=ldref1.masses+ldref2.masses+ldref3.masses+[12.0110, 1.0080, 12.0110, 12.0110, 12.0110]
else:
    ldout.masses=ldref1.masses+ldref2.masses+ldref3.masses

ldout.natoms=len(ldout.atoms)
ldout.nbonds=len(ldout.bonds)
ldout.nangles=len(ldout.angles)
ldout.ndihedrals=len(ldout.dihedrals)
ldout.nimpropers=len(ldout.impropers)

ldout.writeFile("fCOORD")
