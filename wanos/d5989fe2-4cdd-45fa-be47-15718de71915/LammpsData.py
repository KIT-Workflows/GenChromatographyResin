##
## Python class containing LAMMPS data
##

import math

class DataFileError(Exception):
    pass

class LammpsData:
    def __init__(self):
        self.comment = "Created with LammpsData"
        self.natoms = 1
        self.nbonds = 1
        self.nangles = 1
        self.ndihedrals = 1
        self.nimpropers = 1
        self.natomtypes = 1
        self.nbondtypes = 1
        self.nangletypes = 1
        self.ndihedraltypes = 1
        self.nimpropertypes = 1
        self.box = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        self.masses = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []

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
                self.box[0][1] = float(elements[1])
                return

            if ' '.join(elements[2:]) == "ylo yhi":
                self.box[1][0] = float(elements[0])
                self.box[1][1] = float(elements[1])
                return

            if ' '.join(elements[2:]) == "zlo zhi":
                self.box[2][0] = float(elements[0])
                self.box[2][1] = float(elements[1])
                return

        #if len(elements) == 6 and ' '.join(elements[3:]) == "xy xz yz":
        #    self.box[2][0] = float(elements[0])
        #    self.box[2][1] = float(elements[1])
        #    self.box[2][2] = float(elements[2])
        #    return

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
            if len(elements) == 7:
                if int(elements[0]) != len(self.atoms)+1:
                    raise DataFileError("Atoms need to be ordered.")

                
                molid = int(elements[1])
                atomtype = int(elements[2])
                q = float(elements[3])
                x = float(elements[4])
                y = float(elements[5])
                z = float(elements[6])

                self.atoms.append([molid, atomtype, q, x, y, z])
                return

        if self.currentsection == "Bonds":
            if len(elements) == 4:
                if int(elements[0]) != len(self.bonds)+1:
                    raise DataFileError("Bonds need to be ordered.")

                bondtype = int(elements[1])
                n1 = int(elements[2])
                n2 = int(elements[3])

                self.bonds.append([bondtype, n1, n2])
                return

        if self.currentsection == "Angles":
            if len(elements) == 5:
                if int(elements[0]) != len(self.angles)+1:
                    raise DataFileError("Angles need to be ordered.")

                angletype = int(elements[1])
                a1 = int(elements[2])
                a2 = int(elements[3])
                a3 = int(elements[4])

                self.angles.append([angletype, a1, a2, a3])
                return

        if self.currentsection == "Dihedrals":
            if len(elements) == 6:
                if int(elements[0]) != len(self.dihedrals)+1:
                    raise DataFileError("Dihedrals need to be ordered.")

                dihedraltype = int(elements[1])
                d1 = int(elements[2])
                d2 = int(elements[3])
                d3 = int(elements[4])
                d4 = int(elements[5])

                self.dihedrals.append([dihedraltype, d1, d2, d3, d4])
                return

        if self.currentsection == "Impropers":
            if len(elements) == 6:
                if int(elements[0]) != len(self.impropers)+1:
                    raise DataFileError("Impropers need to be ordered.")

                imptype = int(elements[1])
                i1 = int(elements[2])
                i2 = int(elements[3])
                i3 = int(elements[4])
                i4 = int(elements[5])

                self.impropers.append([imptype, i1, i2, i3, i4])
                return
        if self.currentsection == "NA":
            # do nothing
            return

        raise DataFileError("Unable to handle line No. {}: {}".format(ln, line))

    def validateData(self):
        if self.natoms != len(self.atoms):
            raise DataFileError("Number of atoms does not match declared number of atoms.")
        if self.natomtypes != len(self.masses):
            raise DataFileError("Masses not specified for all atom types")
        if self.nbonds != len(self.bonds):
            raise DataFileError("Number of bonds does not match declared number of bonds.")
        if self.nangles != len(self.angles):
            raise DataFileError("Number of angles does not match declared number of angles.")
        if self.ndihedrals != len(self.dihedrals):
            raise DataFileError("Number of dihedrals does not match declared number of dihedrals.")
        if self.nimpropers != len(self.impropers):
            raise DataFileError("Number of impropers does not match declared number of impropers.")

    def readFile(self, infile):
        self.currentSection = "Header"

        with open(infile, 'r') as df:
            for ln, line in enumerate(df):
                try:
                    self.analyseLine(ln, line)
                except DataFileError as err:
                    print("ERROR: " + str(err))
                
        try:
            self.validateData()
        except DataFileError as err:
            print("\nERROR!\n")
            print("The data file read from {} is inconsistent!".format(infile))
            print(err)

    def writeFile(self, outfile):
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

            of.write("{:.5f}  {:.5f} xlo xhi\n".format(self.box[0][0], self.box[0][1]))
            of.write("{:.5f}  {:.5f} ylo yhi\n".format(self.box[1][0], self.box[1][1]))
            of.write("{:.5f}  {:.5f} zlo zhi\n\n".format(self.box[2][0], self.box[2][1]))
            #of.write("{:.5f}  {:.5f}  {:.5f} xy xz yz\n\n".format(
            #         self.box[2][0], self.box[2][1], self.box[2][2]))

            of.write("Masses\n\n")

            mw = int(math.log(self.natomtypes, 10) + 1)

            for i, m in enumerate(self.masses):
                of.write("{:{mw}d} {:7.4f}\n".format(i+1,m,mw=mw))

            of.write("\nAtoms\n\n")

            wid = int(math.log(self.natoms, 10) + 1)
            wmol = int(math.log(max([a[0] for a in self.atoms]), 10) + 1)
            wtype = int(math.log(max([a[1] for a in self.atoms]), 10) + 1)

            for i, a in enumerate(self.atoms):
                of.write("{:{wid}d} {:{wmol}d} {:{wtype}d} {:=9.6f} {:=9.6f} {:=9.6f} {:=9.6f}\n".format(i+1, *a,
                         wid=wid, wmol=wmol, wtype=wtype))

            of.write("\nBonds\n\n")

            wid = int(math.log(self.nbonds, 10) + 1)
            wtype = int(math.log(self.nbondtypes, 10) + 1)
            wb = int(math.log(self.natoms, 10) + 1)

            for i, b in enumerate(self.bonds):
                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d}\n".format(i+1, *b,
                         wid=wid, wtype=wtype, wb=wb))

            of.write("\nAngles\n\n")

            wid = int(math.log(self.nangles, 10) + 1)
            wtype = int(math.log(self.nangletypes, 10) + 1)

            for i, a in enumerate(self.angles):
                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(i+1, *a,
                         wid=wid, wtype=wtype, wb=wb))

            of.write("\nDihedrals\n\n")

            wid = int(math.log(self.ndihedrals, 10) + 1)
            wtype = int(math.log(self.ndihedraltypes, 10) + 1)

            for i, d in enumerate(self.dihedrals):
                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(i+1, *d,
                         wid=wid, wtype=wtype, wb=wb))

            if(self.nimpropers!=0):
                of.write("\nImpropers\n\n")

            wid = int(math.log(self.nimpropers, 10) + 1)
            wtype = int(math.log(self.nimpropertypes, 10) + 1)

            for i, imp in enumerate(self.impropers):
                of.write("{:{wid}d} {:{wtype}d} {:{wb}d} {:{wb}d} {:{wb}d} {:{wb}d}\n".format(
                    i+1, *imp, wid=wid, wtype=wtype, wb=wb))