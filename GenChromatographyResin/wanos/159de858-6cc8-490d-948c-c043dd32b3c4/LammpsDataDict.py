##
## Python class containing LAMMPS data
##

## lmodan edited on apr 2021. omitting velocities and coeffs written by vmd or lammps

import math
import numpy as np

class DataFileError(Exception):
    pass

class LammpsData:
    def __init__(self):
        self.comment = "Created with LammpsData"
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.nimpropers = 0
        self.natomtypes = 0
        self.nbondtypes = 0
        self.nangletypes = 0
        self.ndihedraltypes = 0
        self.nimpropertypes = 0
        self.box = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        self.masses = []
        self.atoms = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.impropers = {}
        self.reorder = False

    def removeComments(self, line):
        line = line.split('#', 1)[0]
        return line.strip()

    def analyseLine(self, ln, line):
        line = line.strip()
        line = self.removeComments(line)
        elements = line.split()
        if not line:
            return

        if ln == 0:
            self.comment = line
            return

        if len(elements) == 2:
            if elements[1] == "atoms":
                self.natoms = int(elements[0])
                return
    
            if elements[1] == "bonds":
                self.nbonds = int(elements[0])
                return

            if elements[1] == "angles":
                self.nangles = int(elements[0])
                return

            if elements[1] == "dihedrals":
                self.ndihedrals = int(elements[0])
                return

            if elements[1] == "impropers":
                self.nimpropers = int(elements[0])
                return

        if len(elements) == 3:
            if ' '.join(elements[1:]) == "atom types":
                self.natomtypes = int(elements[0])
                return

            if ' '.join(elements[1:]) == "bond types":
                self.nbondtypes = int(elements[0])
                return

            if ' '.join(elements[1:]) == "angle types":
                self.nangletypes = int(elements[0])
                return

            if ' '.join(elements[1:]) == "dihedral types":
                self.ndihedraltypes = int(elements[0])
                return

            if ' '.join(elements[1:]) == "improper types":
                self.nimpropertypes = int(elements[0])
                return

        if len(elements) == 4:
            if ' '.join(elements[2:]) == "xlo xhi":
                self.box[0][0] = float(elements[0])
                self.box[1][0] = float(elements[1])
                return

            if ' '.join(elements[2:]) == "ylo yhi":
                self.box[0][1] = float(elements[0])
                self.box[1][1] = float(elements[1])
                return

            if ' '.join(elements[2:]) == "zlo zhi":
                self.box[0][2] = float(elements[0])
                self.box[1][2] = float(elements[1])
                return

        if len(elements) == 6 and ' '.join(elements[3:]) == "xy xz yz":
            self.box[2][0] = float(elements[0])
            self.box[2][1] = float(elements[1])
            self.box[2][2] = float(elements[2])
            return

        if len(elements) == 1:
            if elements[0] == 'Masses':
                self.currentsection = "Masses"
                return

            if elements[0] == 'Atoms':
                self.currentsection = "Atoms"
                return

            if elements[0] == 'Bonds':
                self.currentsection = "Bonds"
                return

            if elements[0] == 'Angles':
                self.currentsection = "Angles"
                return

            if elements[0] == 'Dihedrals':
                self.currentsection = "Dihedrals"
                return

            if elements[0] == 'Impropers':
                self.currentsection = "Impropers"
                return

            if elements[0] == 'Velocities':
                self.currentsection = "NA"
                return

        elif len(elements) == 2 and elements[1] == "Coeffs":
            self.currentsection = "NA"
            return
        
        if self.currentsection == "Masses":
            if len(elements) == 2:
                if int(elements[0]) != len(self.masses)+1:
                    raise DataFileError("Masses need to be ordered.")
                
                self.masses.append(float(elements[1]))
                return

        if self.currentsection == "Atoms":
            if len(elements) == 10:
                #if int(elements[0]) != len(self.atoms)+1:
                #    raise DataFileError("Atoms need to be ordered.")
                aid = int(elements[0])
                
                molid = int(elements[1])
                atomtype = int(elements[2])
                q = float(elements[3])
                x = float(elements[4])
                y = float(elements[5])
                z = float(elements[6])
                cx = int(elements[7])
                cy = int(elements[8])
                cz = int(elements[9])
            elif len(elements) == 7:
                aid = int(elements[0])
                
                molid = int(elements[1])
                atomtype = int(elements[2])
                q = float(elements[3])
                x = float(elements[4])
                y = float(elements[5])
                z = float(elements[6])
                cx = 0
                cy = 0
                cz = 0

            self.atoms[aid]=[molid, atomtype, q, x, y, z, cx, cy, cz]
            return

        if self.currentsection == "Bonds":
            if len(elements) == 4:
                #if int(elements[0]) != len(self.bonds)+1:
                #    raise DataFileError("Bonds need to be ordered.")
                bid = int(elements[0])
                bondtype = int(elements[1])
                n1 = int(elements[2])
                n2 = int(elements[3])

                self.bonds[bid]=[bondtype, n1, n2]
                return

        if self.currentsection == "Angles":
            if len(elements) == 5:
                #if int(elements[0]) != len(self.angles)+1:
                #    raise DataFileError("Angles need to be ordered.")
                cid = int(elements[0])
                angletype = int(elements[1])
                a1 = int(elements[2])
                a2 = int(elements[3])
                a3 = int(elements[4])

                self.angles[cid]=[angletype, a1, a2, a3]
                return

        if self.currentsection == "Dihedrals":
            if len(elements) == 6:
                #if int(elements[0]) != len(self.dihedrals)+1:
                #    raise DataFileError("Dihedrals need to be ordered.")
                did = int(elements[0])
                dihedraltype = int(elements[1])
                d1 = int(elements[2])
                d2 = int(elements[3])
                d3 = int(elements[4])
                d4 = int(elements[5])

                self.dihedrals[did]=[dihedraltype, d1, d2, d3, d4]
                return

        if self.currentsection == "Impropers":
            if len(elements) == 6:
                #if int(elements[0]) != len(self.impropers)+1:
                #    raise DataFileError("Impropers need to be ordered.")
                eid = int(elements[0])
                imptype = int(elements[1])
                i1 = int(elements[2])
                i2 = int(elements[3])
                i3 = int(elements[4])
                i4 = int(elements[5])

                self.impropers[eid]=[imptype, i1, i2, i3, i4]
                return
        if self.currentsection == "NA":
            # do nothing
            return

        raise DataFileError("Unable to handle line No. {}: {}".format(ln, line))

    def validateData(self):
        #if self.natoms != len(self.atoms):
        #    raise DataFileError("Number of atoms does not match declared number of atoms.")
        #if self.natomtypes != len(self.masses):
        #    raise DataFileError("Masses not specified for all atom types")
        #if self.nbonds != len(self.bonds):
        #    raise DataFileError("Number of bonds does not match declared number of bonds.")
        #if self.nangles != len(self.angles):
        #    raise DataFileError("Number of angles does not match declared number of angles.")
        #if self.ndihedrals != len(self.dihedrals):
        #    raise DataFileError("Number of dihedrals does not match declared number of dihedrals.")
        #if self.nimpropers != len(self.impropers):
        #    raise DataFileError("Number of impropers does not match declared number of impropers.")
        if self.reorder:
            self.natoms=len(self.atoms)
        else:
            self.natoms=np.max([k for i,(k,atom) in enumerate(self.atoms.items())])
        self.nbonds=len(self.bonds)
        self.nangles=len(self.angles)
        self.ndihedrals=len(self.dihedrals)
        self.nimpropers=len(self.impropers)
        self.natomtypes=np.max([atom[1] for i,(k,atom) in enumerate(self.atoms.items())])
        self.nbondtypes=np.max([atom[0] for i,(k,atom) in enumerate(self.bonds.items())])
        if len(self.angles):
            self.nangletypes=np.max([atom[0] for i,(k,atom) in enumerate(self.angles.items())])
        else:
            self.nangletypes=0
        if len(self.dihedrals):
            self.ndihedraltypes=np.max([atom[0] for i,(k,atom) in enumerate(self.dihedrals.items())])
        else:
            self.ndihedrals=0
        if len(self.impropers):
            self.nimpropertypes=np.max([atom[0] for i,(k,atom) in enumerate(self.impropers.items())])
        else:
            self.nimpropers=0

    def readFile(self, infile):
        self.currentSection = "Header"

        with open(infile, 'r') as df:
            for ln, line in enumerate(df):
                try:
                    self.analyseLine(ln, line)
                except DataFileError as err:
                    print("ERROR: " + str(err))
                
        #try:
        self.validateData()
        #except DataFileError as err:
        #    print("\nERROR!\n")
        #    print("The data file read from {} is inconsistent!".format(infile))
        #    print(err)

    def writeFile(self, outfile, reorder=False):
        self.reorder=reorder
        self.validateData()
        with open(outfile, 'w') as of:
            of.write("{}\n\n".format(self.comment))

            maxnum = max([self.natoms, self.nbonds, self.nangles,
                          self.ndihedrals, self.nimpropers])
            mw = int(math.log(maxnum, 10) + 1)

            of.write("{:{mw}d} atoms\n".format(self.natoms, mw=mw))
            of.write("{:{mw}d} bonds\n".format(self.nbonds, mw=mw))
            of.write("{:{mw}d} angles\n".format(self.nangles, mw=mw))
            of.write("{:{mw}d} dihedrals\n".format(self.ndihedrals, mw=mw))
            of.write("{:{mw}d} impropers\n\n".format(self.nimpropers, mw=mw))

            maxnum = max([self.natomtypes, self.nbondtypes, self.nangletypes,
                          self.ndihedraltypes, self.nimpropertypes])
            mw = int(math.log(maxnum, 10) + 1)

            of.write("{:{mw}d} atom types\n".format(self.natomtypes, mw=mw))
            of.write("{:{mw}d} bond types\n".format(self.nbondtypes, mw=mw))
            of.write("{:{mw}d} angle types\n".format(self.nangletypes, mw=mw))
            of.write("{:{mw}d} dihedral types\n".format(self.ndihedraltypes, mw=mw))
            of.write("{:{mw}d} improper types\n\n".format(self.nimpropertypes, mw=mw))

            of.write("{:.5f}  {:.5f} xlo xhi\n".format(self.box[0][0], self.box[1][0]))
            of.write("{:.5f}  {:.5f} ylo yhi\n".format(self.box[0][1], self.box[1][1]))
            of.write("{:.5f}  {:.5f} zlo zhi\n".format(self.box[0][2], self.box[1][2]))
            #of.write("{:.5f}  {:.5f}  {:.5f} xy xz yz\n\n".format(
            #         self.box[2][0], self.box[2][1], self.box[2][2]))

            of.write("\nMasses\n\n")

            mw = int(math.log(self.natomtypes, 10) + 1)

            for i, m in enumerate(self.masses):
                of.write("{:{mw}d} {:7.4f}\n".format(i+1,m,mw=mw))

            of.write("\nAtoms\n\n")

            wid = int(math.log(self.natoms, 10) + 1)
            wmol = int(math.log(max([a[0] for i,(k,a) in enumerate(self.atoms.items())]), 10) + 1)
            wtype = int(math.log(max([a[1] for i,(k,a) in enumerate(self.atoms.items())]), 10) + 1)

            #for i,(k,a) in enumerate(self.atoms.items()):
            ks=sorted(self.atoms)
            mapatom={}

            if reorder:
                for i, k in enumerate(ks):
                        a=self.atoms[k]
                        of.write("{:{wid}d} {:{wmol}d} {:{wtype}d} {:=9.6f} {:=9.6f} {:=9.6f} {:=9.6f}\n".format(i+1, *a,
                                 wid=wid, wmol=wmol, wtype=wtype))

                        mapatom[k]=i+1

            else:
                for i, k in enumerate(ks):
                    a=self.atoms[k]
                    of.write("{:{wid}d} {:{wmol}d} {:{wtype}d} {:=9.6f} {:=9.6f} {:=9.6f} {:=9.6f}\n".format(k, *a,
                             wid=wid, wmol=wmol, wtype=wtype))

            of.write("\nBonds\n\n")

            wid = int(math.log(self.nbonds, 10) + 1)
            wtype = int(math.log(self.nbondtypes, 10) + 1)
            wb = int(math.log(self.natoms, 10) + 1)

            #for i,(k,b) in enumerate(self.bonds.items()):
            ks=sorted(self.bonds)

            for i, k in enumerate(ks):
                b=self.bonds[k]

                if reorder:
                    b[1]=mapatom[b[1]]
                    b[2]=mapatom[b[2]]

                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d}\n".format(k, *b,
                         wid=wid, wtype=wtype, wb=wb))

            if len(self.angles):

                of.write("\nAngles\n\n")

                wid = int(math.log(self.nangles, 10) + 1)
                wtype = int(math.log(self.nangletypes, 10) + 1)

                #for i,(k,a) in enumerate(self.angles.items()):
                ks=sorted(self.angles)

                for i, k in enumerate(ks):
                    a=self.angles[k]

                    if reorder:
                        a[1]=mapatom[a[1]]
                        a[2]=mapatom[a[2]]
                        a[3]=mapatom[a[3]]

                    of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(k, *a,
                             wid=wid, wtype=wtype, wb=wb))
            if len(self.dihedrals):
                of.write("\nDihedrals\n\n")

                wid = int(math.log(self.ndihedrals, 10) + 1)
                wtype = int(math.log(self.ndihedraltypes, 10) + 1)

                #for i,(k,d) in enumerate(self.dihedrals.items()):
                ks=sorted(self.dihedrals)

                for i, k in enumerate(ks):
                    d=self.dihedrals[k]

                    if reorder:
                        d[1]=mapatom[d[1]]
                        d[2]=mapatom[d[2]]
                        d[3]=mapatom[d[3]]
                        d[4]=mapatom[d[4]]

                    of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(k, *d,
                             wid=wid, wtype=wtype, wb=wb))
            if len(self.impropers):
                of.write("\nImpropers\n\n")

                wid = int(math.log(self.nimpropers, 10) + 1)
                wtype = int(math.log(self.nimpropertypes, 10) + 1)

                #for i,(k,imp) in enumerate(self.impropers.items()):
                ks=sorted(self.impropers)

                for i, k in enumerate(ks):
                    imp=self.impropers[k]

                    if reorder:
                        imp[1]=mapatom[imp[1]]
                        imp[2]=mapatom[imp[2]]
                        imp[3]=mapatom[imp[3]]
                        imp[4]=mapatom[imp[4]]


                    of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(k, *imp, wid=wid, wtype=wtype, wb=wb))
