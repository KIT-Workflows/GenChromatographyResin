##
## Python class containing LAMMPS data
##

import math

class DataFileError(Exception):
    pass

class LammpsDataFF:
    def __init__(self):
        #self.comment = "Created with LammpsData"
        self.natomtypes = 0
        self.nbondtypes = 0
        self.nangletypes = 0
        self.ndihedraltypes = 0
        self.nimpropertypes = 0
        self.masses = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []

    def removeComments(self, line):
        if "#" in line:
            line = line.split('#')[1]
        return line.strip()

    def analyseLine(self, ln, line):
        line = line.strip()
        line = self.removeComments(line)
#        elements = line.split()
        if not line:
            return

        if True:
            if 'Masses' in line:
                self.currentsection = "Masses"
                return

            if 'Pair Coeffs' in line:
                self.currentsection = "Atoms"
                return

            if 'Bond Coeffs' in line:
                self.currentsection = "Bonds"
                return

            if 'Angle Coeffs' in line:
                self.currentsection = "Angles"
                return

            if 'Dihedral Coeffs' in line:
                self.currentsection = "Dihedrals"
                return

            if 'Improper Coeffs' in line:
                self.currentsection = "Impropers"
                return
        
        if self.currentsection == "Masses":
                self.masses.append(float(line.split()[1]))
                return

        if self.currentsection == "Atoms":
                self.atoms.append(line.split()[1:])
                return

        if self.currentsection == "Bonds":
                self.bonds.append(line.split()[1:])
                return

        if self.currentsection == "Angles":
                self.angles.append(line.split()[1:])
                return

        if self.currentsection == "Dihedrals":
                self.dihedrals.append(line.split()[1:])
                return

        if self.currentsection == "Impropers":
                self.impropers.append(line.split()[1:])
                return

        raise DataFileError("Unable to handle line No. {}: {}".format(ln, line))

    def validateData(self):
        self.natomtypes=len(self.atoms)
        self.nbondtypes=len(self.bonds)
        self.nangletypes=len(self.angles)
        self.ndihedraltypes=len(self.dihedrals)
        self.nimpropertypes=len(self.impropers)
        return

    def readFile(self, infile):
        self.currentSection = "Header"

        with open(infile, 'r') as df:
            for ln, line in enumerate(df):
                try:
                    self.analyseLine(ln, line)
                except DataFileError as err:
                    print("ERROR: " + str(err))
                
#        try:
        self.validateData()
#        except DataFileError as err:
#            print("\nERROR!\n")
#            print("The data file read from {} is inconsistent!".format(infile))
#            print(err)

    def writeFile(self, outfile, autoID=True):
        if autoID:
            self.validateData()

        with open(outfile, 'w') as of:
#            of.write("{}\n\n".format(self.comment))

#            maxnum = max([self.natoms, self.nbonds, self.nangles,
#                          self.ndihedrals, self.nimpropers])
#            mw = int(math.log(maxnum, 10) + 1)

#            of.write("{:{mw}d} atoms\n".format(self.natoms, mw=mw))
#            of.write("{:{mw}d} bonds\n".format(self.nbonds, mw=mw))
#            of.write("{:{mw}d} angles\n".format(self.nangles, mw=mw))
#            of.write("{:{mw}d} dihedrals\n".format(self.ndihedrals, mw=mw))
#            of.write("{:{mw}d} impropers\n\n".format(self.nimpropers, mw=mw))

#            maxnum = max([self.natomtypes, self.nbondtypes, self.nangletypes,
#                          self.ndihedraltypes, self.nimpropertypes])
#            mw = int(math.log(maxnum, 10) + 1)

#            of.write("{:{mw}d} atom types\n".format(self.natomtypes, mw=mw))
#            of.write("{:{mw}d} bond types\n".format(self.nbondtypes, mw=mw))
#            of.write("{:{mw}d} angle types\n".format(self.nangletypes, mw=mw))
#            of.write("{:{mw}d} dihedral types\n".format(self.ndihedraltypes, mw=mw))
#            of.write("{:{mw}d} improper types\n\n".format(self.nimpropertypes, mw=mw))

#            of.write("{:.5f}  {:.5f} xlo xhi\n".format(self.box[0][0], self.box[0][1]))
#            of.write("{:.5f}  {:.5f} ylo yhi\n".format(self.box[1][0], self.box[1][1]))
#            of.write("{:.5f}  {:.5f} zlo zhi\n".format(self.box[2][0], self.box[2][1]))
            #of.write("{:.5f}  {:.5f}  {:.5f} xy xz yz\n\n".format(
            #         self.box[2][0], self.box[2][1], self.box[2][2]))

#            of.write("\nMasses\n\n")

#            mw = int(math.log(self.natomtypes, 10) + 1)

#            for i, m in enumerate(self.masses):
#                of.write("{:{mw}d} {:7.4f}\n".format(i+1,m,mw=mw))

            of.write("Pair Coeffs\n\n")

#            wid = int(math.log(self.natoms, 10) + 1)
#            wmol = int(math.log(max([a[0] for a in self.atoms]), 10) + 1)
#            wtype = int(math.log(max([a[1] for a in self.atoms]), 10) + 1)

            for i, a in enumerate(self.atoms):
#                of.write("{:{wid}d} {:{wmol}d} {:{wtype}d} {:=9.6f} {:=9.6f} {:=9.6f} {:=9.6f}\n".format(i+1, *a,
#                         wid=wid, wmol=wmol, wtype=wtype))
                of.write("{} {}\n".format(i+1, " ".join(a)))

            of.write("\nBond Coeffs\n\n")

#            wid = int(math.log(self.nbonds, 10) + 1)
#            wtype = int(math.log(self.nbondtypes, 10) + 1)
#            wb = int(math.log(self.natoms, 10) + 1)

            for i, b in enumerate(self.bonds):
#                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d}\n".format(i+1, *b,
#                         wid=wid, wtype=wtype, wb=wb))
                of.write("{} {}\n".format(i+1, " ".join(b)))

            of.write("\nAngle Coeffs\n\n")

#            wid = int(math.log(self.nangles, 10) + 1)
#            wtype = int(math.log(self.nangletypes, 10) + 1)

            for i, a in enumerate(self.angles):
#                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(i+1, *a,
#                         wid=wid, wtype=wtype, wb=wb))
                of.write("{} {}\n".format(i+1, " ".join(a)))

            of.write("\nDihedral Coeffs\n\n")

#            wid = int(math.log(self.ndihedrals, 10) + 1)
#            wtype = int(math.log(self.ndihedraltypes, 10) + 1)

            for i, d in enumerate(self.dihedrals):
#                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(i+1, *d,
#                         wid=wid, wtype=wtype, wb=wb))
                of.write("{} {}\n".format(i+1, " ".join(d)))

            of.write("\nImproper Coeffs\n\n")

#                wid = int(math.log(self.nimpropers, 10) + 1)
#                wtype = int(math.log(self.nimpropertypes, 10) + 1)

            for i, imp in enumerate(self.impropers):
#                    of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(i+1, *imp, wid=wid, wtype=wtype, wb=wb))
                of.write("{} {}\n".format(i+1, " ".join(imp)))